# =============================================================================
#  axi_solver.jl — the corrected finite Eshelby cell, axisymmetric declination.
#
#  The general 2 × N scheme of Adessina, Barthélémy, Lavergne & Ben Fraj
#  (*Int. J. Eng. Sci.* 119, 2017), one Fourier mode at a time.  For each mode
#  and each Kelvin basis tensor of that mode, two problems are solved on the
#  truncated cell:
#
#    E-problem   u|∂Ω = E·x                        → ⟨ε⟩ = 𝔸ᴱ:E, ⟨σ⟩ = 𝔹ᴱ:E
#    p-problem   u|∂Ω = +∇G(x):(V_D P)             → ⟨ε⟩ = 𝔸ᵖ:P, ⟨σ⟩ = 𝔹ᵖ:P
#
#  The polarization of the *infinite*-medium problem is its own consequence,
#
#      P = ⟨σ − ℂ₀:ε⟩_D = (𝔹ᴱ − ℂ₀:𝔸ᴱ):E + (𝔹ᵖ − ℂ₀:𝔸ᵖ):P ,
#
#  a linear fixed point solved in closed form by
#
#      X = [𝕀 − (𝔹ᵖ − ℂ₀:𝔸ᵖ)]⁻¹ : (𝔹ᴱ − ℂ₀:𝔸ᴱ) ,
#      𝔸 = 𝔸ᴱ + 𝔸ᵖ:X ,      𝔹 = 𝔹ᴱ + 𝔹ᵖ:X .
#
#  Transport is the same algebra with `ε ↦ ∇T`, `σ ↦ k∇T`, `ℂ₀ ↦ k₀` and the
#  scalar dipole field `T = M·x/(4π k₀ r³)`.
#
#  Cost: three assemblies (modes 0, 1, 2) and eight solves in *two* dimensions,
#  where the three-dimensional formulation of the paper needs six solves on a
#  mesh two orders of magnitude larger.
# =============================================================================

"One Fourier mode, discretized. Built once per inclusion and reused for every reference medium."
struct AxiModeSetup{D, C}
    m::Int
    ncomp::Int
    dh::D
    cv::C
    perm::Vector{Int}                    # component-major local dof ordering
    ch::Ferrite.ConstraintHandler        # outer boundary, data behind `bcref`
    bcref::Base.RefValue{Any}
    axis::Vector{Int}
    presc::Vector{Int}
    free::Vector{Int}
end

"Whole discretization of an axisymmetric inclusion: the grid and its modes."
struct AxiSetup{G}
    grid::G
    modes::Dict{Tuple{Symbol, Int}, Any}
end

_axi_setup(incl::MeanFieldHom.FEExcenteredSphere) =
    incl.cache.setup === nothing ?
    (incl.cache.setup = AxiSetup(_axi_grid(incl), Dict{Tuple{Symbol, Int}, Any}())) :
    incl.cache.setup

const _AXI_FIELDS = (:c1, :c2, :c3)

function _axi_mode_setup(incl::MeanFieldHom.FEExcenteredSphere, physics::Symbol, m::Int)
    s = _axi_setup(incl)
    return get!(s.modes, (physics, m)) do
        _build_axi_mode(s.grid, incl.mesh.order, physics, m)
    end
end

function _build_axi_mode(grid, order::Int, physics::Symbol, m::Int)
    ncomp = physics === :elasticity ? _axi_ncomp(m) : 1
    fields = _AXI_FIELDS[1:ncomp]

    ip = Ferrite.Lagrange{Ferrite.RefTriangle, order}()
    ip_geo = Ferrite.Lagrange{Ferrite.RefTriangle, 1}()
    qr = Ferrite.QuadratureRule{Ferrite.RefTriangle}(2 * order + 1)
    cv = Ferrite.CellValues(qr, ip, ip_geo)

    dh = Ferrite.DofHandler(grid)
    for f in fields
        Ferrite.add!(dh, f, ip)
    end
    Ferrite.close!(dh)
    perm = reduce(vcat, [collect(Ferrite.dof_range(dh, f)) for f in fields])

    # Outer boundary: values supplied per solve through `bcref`.
    bcref = Base.RefValue{Any}(_ -> ntuple(_ -> 0.0, ncomp))
    ch = Ferrite.ConstraintHandler(dh)
    for (k, f) in pairs(fields)
        Ferrite.add!(
            ch, Ferrite.Dirichlet(
                f, Ferrite.getfacetset(grid, AXI_SET_OUTER), (x, _t) -> bcref[](x)[k]
            )
        )
    end
    Ferrite.close!(ch)

    # Axis: the components regularity forces to vanish.
    zeros_ = physics === :elasticity ? _axi_axis_zeros(m) : (m == 0 ? () : (1,))
    axis = Int[]
    if !isempty(zeros_)
        cha = Ferrite.ConstraintHandler(dh)
        for k in zeros_
            Ferrite.add!(
                cha, Ferrite.Dirichlet(
                    fields[k], Ferrite.getfacetset(grid, AXI_SET_AXIS), (_x, _t) -> 0.0
                )
            )
        end
        Ferrite.close!(cha)
        axis = collect(cha.prescribed_dofs)
    end

    presc = sort!(union(collect(ch.prescribed_dofs), axis))
    free = setdiff(1:Ferrite.ndofs(dh), presc)
    return AxiModeSetup(m, ncomp, dh, cv, perm, ch, bcref, axis, presc, free)
end

# ─── Assembly ────────────────────────────────────────────────────────────────

"""
    _axi_assemble(ms, grid, Dmap, Bop) -> SparseMatrixCSC

Stiffness of one Fourier mode. `Dmap` maps a cell-set name to the material
matrix of that region, in the cylindrical `(ρ, θ, z)` basis; `Bop` builds the
generalized-strain operator of one scalar shape function. Both the elastic
(6 rows) and the transport (3 rows) problems go through this one routine.

The measure is `ρ dρ dz`: that single factor is what turns a plane problem
into a solid of revolution.
"""
function _axi_assemble(ms::AxiModeSetup, grid, Dmap, Bop)
    K = Ferrite.allocate_matrix(ms.dh)
    asm = Ferrite.start_assemble(K)
    nb = Ferrite.getnbasefunctions(ms.cv)
    nl = nb * ms.ncomp
    ke = zeros(nl, nl)
    Bq = zeros(size(first(values(Dmap)), 1), nl)

    for (setname, D) in Dmap
        for cell in Ferrite.CellIterator(ms.dh, Ferrite.getcellset(grid, setname))
            Ferrite.reinit!(ms.cv, cell)
            coords = Ferrite.getcoordinates(cell)
            fill!(ke, 0)
            for q in 1:Ferrite.getnquadpoints(ms.cv)
                ρ = Ferrite.spatial_coordinate(ms.cv, q, coords)[1]
                dΩ = Ferrite.getdetJdV(ms.cv, q) * ρ
                fill!(Bq, 0)
                for i in 1:nb
                    N = Ferrite.shape_value(ms.cv, q, i)
                    ∇N = Ferrite.shape_gradient(ms.cv, q, i)
                    Bi = Bop(ms.m, N, ∇N[1], ∇N[2], ρ)
                    for c in 1:ms.ncomp, r in axes(Bi, 1)
                        Bq[r, (c - 1) * nb + i] = Bi[r, c]
                    end
                end
                ke .+= (Bq' * D * Bq) .* dΩ
            end
            Ferrite.assemble!(asm, Ferrite.celldofs(cell)[ms.perm], ke)
        end
    end
    return K
end

# ─── Averages over the inclusion ─────────────────────────────────────────────

"""
    _axi_averages(ms, grid, Dmap, u, Bop, proj) -> (⟨primal⟩, ⟨dual⟩, V)

Volume averages over the inclusion (the `"core"` and `"shell"` regions) of the
generalized strain and of the associated generalized stress, projected on the
Kelvin basis of the mode, plus the inclusion volume `V = 2π ∫∫ ρ dρ dz`.

The azimuthal integration has already been performed analytically in `proj`;
what remains is the meridian quadrature with the `ρ dρ dz` measure.
"""
function _axi_averages(ms::AxiModeSetup, grid, Dmap, u, Bop, proj)
    nb = Ferrite.getnbasefunctions(ms.cv)
    nl = nb * ms.ncomp
    nrow = size(first(values(Dmap)), 1)
    Bq = zeros(nrow, nl)
    prim = zeros(length(proj(ms.m, zeros(nrow))))
    dual = zeros(length(prim))
    V = 0.0

    for setname in (AXI_SET_CORE, AXI_SET_SHELL)
        D = Dmap[setname]
        for cell in Ferrite.CellIterator(ms.dh, Ferrite.getcellset(grid, setname))
            Ferrite.reinit!(ms.cv, cell)
            coords = Ferrite.getcoordinates(cell)
            ue = u[Ferrite.celldofs(cell)[ms.perm]]
            for q in 1:Ferrite.getnquadpoints(ms.cv)
                ρ = Ferrite.spatial_coordinate(ms.cv, q, coords)[1]
                dΩ = Ferrite.getdetJdV(ms.cv, q) * ρ
                fill!(Bq, 0)
                for i in 1:nb
                    N = Ferrite.shape_value(ms.cv, q, i)
                    ∇N = Ferrite.shape_gradient(ms.cv, q, i)
                    Bi = Bop(ms.m, N, ∇N[1], ∇N[2], ρ)
                    for c in 1:ms.ncomp, r in axes(Bi, 1)
                        Bq[r, (c - 1) * nb + i] = Bi[r, c]
                    end
                end
                e = Bq * ue
                prim .+= proj(ms.m, e) .* dΩ
                dual .+= proj(ms.m, D * e) .* dΩ
                V += dΩ
            end
        end
    end
    return prim ./ V, dual ./ V, 2π * V
end

# ─── One mode, both problems ─────────────────────────────────────────────────

"""
    _axi_solve_mode(ms, grid, Dmap, Bop, proj, nload, bc_E, bc_p) -> (Aᴱ, 𝔹ᴱ, 𝔸ᵖ, 𝔹ᵖ, V)

Assemble one mode, factorize once, and run the `nload` affine problems followed
by the `nload` dipole problems. `bc_E(j)` and `bc_p(j, V)` return the boundary
data of load `j` as a function of the boundary point.
"""
function _axi_solve_mode(ms::AxiModeSetup, grid, Dmap, Bop, proj, nload, bc_E, bc_p)
    K = _axi_assemble(ms, grid, Dmap, Bop)
    F = LinearAlgebra.lu(K[ms.free, ms.free])
    Kfp = K[ms.free, ms.presc]
    nd = Ferrite.ndofs(ms.dh)

    function solve_with(f)
        ms.bcref[] = f
        Ferrite.update!(ms.ch, 0.0)
        u = zeros(nd)
        u[ms.ch.prescribed_dofs] .= ms.ch.inhomogeneities
        u[ms.axis] .= 0.0
        u[ms.free] .= F \ Vector(-Kfp * u[ms.presc])
        return u
    end

    A_E = zeros(nload, nload)
    B_E = zeros(nload, nload)
    V = 0.0
    for j in 1:nload
        e, s, V = _axi_averages(ms, grid, Dmap, solve_with(bc_E(j)), Bop, proj)
        A_E[:, j] .= e
        B_E[:, j] .= s
    end

    A_p = zeros(nload, nload)
    B_p = zeros(nload, nload)
    for j in 1:nload
        e, s, _ = _axi_averages(ms, grid, Dmap, solve_with(bc_p(j, V)), Bop, proj)
        A_p[:, j] .= e
        B_p[:, j] .= s
    end
    return A_E, B_E, A_p, B_p, V
end

"""
    _axi_correct(A_E, B_E, A_p, B_p, P0) -> (A, B, X)

The linear fixed point of the corrected boundary condition, on one modal
block: `X = (𝕀 − (𝔹ᵖ − P₀·𝔸ᵖ))⁻¹ (𝔹ᴱ − P₀·𝔸ᴱ)`, then `𝔸 = 𝔸ᴱ + 𝔸ᵖ X` and
`𝔹 = 𝔹ᴱ + 𝔹ᵖ X`.
"""
function _axi_correct(A_E, B_E, A_p, B_p, P0)
    Π_E = B_E - P0 * A_E
    Π_p = B_p - P0 * A_p
    X = (LinearAlgebra.I - Π_p) \ Π_E
    return A_E + A_p * X, B_E + B_p * X, X
end

# ─── Material matrices in the local frame ────────────────────────────────────

"Rotation matrix whose columns are the inclusion's local axes in global coordinates."
_axi_frame(incl) =
    hcat(MeanFieldHom.Core._frame_columns(MeanFieldHom.inclusion_basis(incl))...)

"Components of a 4th-order tensor in the local frame, as a 3×3×3×3 array."
function _to_local4(P::TensND.AbstractTens{4, 3}, R)
    G = MeanFieldHom.Core._C_array(P)
    L = zeros(Float64, 3, 3, 3, 3)
    @inbounds for p in 1:3, q in 1:3, r in 1:3, s in 1:3
        acc = 0.0
        for i in 1:3, j in 1:3, k in 1:3, l in 1:3
            acc += R[i, p] * R[j, q] * R[k, r] * R[l, s] * G[i, j, k, l]
        end
        L[p, q, r, s] = acc
    end
    return L
end

"Components of a 2nd-order tensor in the local frame, as a 3×3 matrix."
_to_local2(P::TensND.AbstractTens{2, 3}, R) = R' * TensND.components_canon(P) * R

"Rebuild a global-frame `Tens{4,3}` from local-frame Kelvin-Mandel components."
function _from_local66(M66::AbstractMatrix{Float64}, R)
    L = MeanFieldHom.Core.array_from_mandel66(M66)
    G = zeros(Float64, 3, 3, 3, 3)
    @inbounds for i in 1:3, j in 1:3, k in 1:3, l in 1:3
        acc = 0.0
        for p in 1:3, q in 1:3, r in 1:3, s in 1:3
            acc += R[i, p] * R[j, q] * R[k, r] * R[l, s] * L[p, q, r, s]
        end
        G[i, j, k, l] = acc
    end
    return TensND.TensCanonical(
        Tensors.Tensor{4, 3}((i, j, k, l) -> G[i, j, k, l])
    )
end

_from_local33(M::AbstractMatrix{Float64}, R) = TensND.TensCanonical(R * M * R')

"""
    _axi_check_ti(M66, what) -> M66

Refuse a constituent whose local-frame components are not transversely
isotropic about the symmetry axis: any other anisotropy couples the Fourier
modes, and the whole point of the axisymmetric formulation is that they do not.
"""
function _axi_check_ti(M66::AbstractMatrix{Float64}, what::AbstractString; rtol = 1.0e-8)
    ref = M66
    scale = maximum(abs, ref)
    # A TI tensor about e₃ has no (11,13), (11,23), (11,12), (33,·shear) coupling
    # and satisfies M₁₁ = M₂₂, M₁₃ = M₂₃, M₄₄ = M₅₅, M₆₆ = M₁₁ − M₁₂.
    bad = max(
        abs(ref[1, 1] - ref[2, 2]), abs(ref[1, 3] - ref[2, 3]),
        abs(ref[4, 4] - ref[5, 5]), abs(ref[6, 6] - (ref[1, 1] - ref[1, 2])),
        maximum(abs, ref[1:3, 4:6]), maximum(abs, ref[4:6, 1:3]),
        abs(ref[4, 5]), abs(ref[4, 6]), abs(ref[5, 6]),
    )
    bad ≤ rtol * scale || throw(
        ArgumentError(
            "$what must be transversely isotropic about the inclusion's " *
                "symmetry axis — an axisymmetric Fourier formulation cannot " *
                "represent any other anisotropy, because the modes would couple. " *
                "Deviation: $(round(bad / scale * 100, sigdigits = 3)) %."
        )
    )
    return ref
end

function _axi_check_ti2(M::AbstractMatrix{Float64}, what::AbstractString; rtol = 1.0e-8)
    scale = maximum(abs, M)
    bad = max(
        abs(M[1, 1] - M[2, 2]), abs(M[1, 2]), abs(M[2, 1]),
        abs(M[1, 3]), abs(M[3, 1]), abs(M[2, 3]), abs(M[3, 2]),
    )
    bad ≤ rtol * scale || throw(
        ArgumentError(
            "$what must be transversely isotropic about the inclusion's " *
                "symmetry axis; deviation $(round(bad / scale * 100, sigdigits = 3)) %."
        )
    )
    return M
end

"""
    _axi_iso_moduli(P₀) -> (μ, ν)

Shear modulus and Poisson ratio of the reference medium, refusing anything
that is not isotropic in *content* — the corrected boundary condition uses the
closed-form Kelvin dipole field, whose anisotropic counterpart (Pan-Chou,
Barnett-Willis) is not implemented.

`rtol` is deliberately loose (0.01 %). The tensors this solver *returns* are
isotropic only to the discretization error — a few parts in a million — so an
iterative scheme fed back its own estimate can never present a reference
isotropic to machine precision. A tighter guard would refuse `SelfConsistent`
and `AsymmetricSelfConsistent` on a perfectly legitimate isotropic problem.
Genuine anisotropy, the case the guard exists for, is orders of magnitude
larger; `IsoSymmetrize` remains the answer for it.
"""
function _axi_iso_moduli(C₀::TensND.AbstractTens{4, 3}; rtol = 1.0e-4)
    C_iso = MeanFieldHom.Core.isotropify(C₀)
    A, Aiso = MeanFieldHom.Core._C_array(C₀), MeanFieldHom.Core._C_array(C_iso)
    dev = maximum(abs, A .- Aiso)
    dev ≤ rtol * max(maximum(abs, Aiso), eps()) || throw(
        ArgumentError(
            "`FEExcenteredSphere` supports an isotropic reference medium only " *
                "(the corrected boundary condition uses the closed-form Kelvin " *
                "dipole field). The reference deviates from isotropy by " *
                "$(round(dev / maximum(abs, Aiso) * 100, sigdigits = 3)) %. Add " *
                "`symmetrize = IsoSymmetrize()` to the phase if this comes from an " *
                "iterative scheme."
        )
    )
    E, ν = MeanFieldHom.Core.extract_iso_moduli(C_iso)
    return Float64(E / (2 * (1 + ν))), Float64(ν)
end

function _axi_iso_scalar(K₀::TensND.AbstractTens{2, 3}; rtol = 1.0e-4)
    M = TensND.components_canon(K₀)
    k = (M[1, 1] + M[2, 2] + M[3, 3]) / 3
    dev = maximum(abs, M .- k .* LinearAlgebra.I(3))
    dev ≤ rtol * max(abs(k), eps()) || throw(
        ArgumentError(
            "`FEExcenteredSphere` supports an isotropic reference medium only; " *
                "the reference conductivity deviates by " *
                "$(round(dev / abs(k) * 100, sigdigits = 3)) %."
        )
    )
    return Float64(k)
end

# ─── The two physics ─────────────────────────────────────────────────────────

"""
    _axi_run_elastic(incl, C₀) -> NamedTuple

Full corrected solve in elasticity: the three modes, the fixed point, and the
two localization tensors reassembled in the global frame.
"""
function _axi_run_elastic(incl::MeanFieldHom.FEExcenteredSphere, C₀::TensND.AbstractTens{4, 3})
    μ, ν = _axi_iso_moduli(C₀)
    R = _axi_frame(incl)
    Dmap = Dict(
        AXI_SET_CORE => _axi_check_ti(
            MeanFieldHom.Core.mandel66_minor(_to_local4(incl.props[1], R)), "the core"
        ),
        AXI_SET_SHELL => _axi_check_ti(
            MeanFieldHom.Core.mandel66_minor(_to_local4(incl.props[2], R)), "the shell"
        ),
        AXI_SET_MATRIX => MeanFieldHom.Core.mandel66_minor(_to_local4(C₀, R)),
    )
    C0_66 = Dmap[AXI_SET_MATRIX]
    grid = _axi_setup(incl).grid

    blocks = map((0, 1, 2)) do m
        ms = _axi_mode_setup(incl, :elasticity, m)
        nload = m == 0 ? 2 : 1
        A_E, B_E, A_p, B_p, V = _axi_solve_mode(
            ms, grid, Dmap, _axi_B_elast, _axi_project, nload,
            j -> (x -> _axi_bc_affine(m, j, x[1], x[2])),
            (j, V) -> (x -> _axi_bc_dipole(m, j, x[1], x[2], μ, ν, V)),
        )
        cols = m == 0 ? (1:2) : m == 1 ? (3:3) : (5:5)
        Qm = _AXI_Q[:, cols]
        A, B, X = _axi_correct(A_E, B_E, A_p, B_p, Qm' * C0_66 * Qm)
        (; m, A_E, B_E, A_p, B_p, A, B, X, V)
    end

    A66 = _axi_assemble_66(blocks[1].A, blocks[2].A[1, 1], blocks[3].A[1, 1])
    B66 = _axi_assemble_66(blocks[1].B, blocks[2].B[1, 1], blocks[3].B[1, 1])
    A66_E = _axi_assemble_66(blocks[1].A_E, blocks[2].A_E[1, 1], blocks[3].A_E[1, 1])
    B66_E = _axi_assemble_66(blocks[1].B_E, blocks[2].B_E[1, 1], blocks[3].B_E[1, 1])
    incl.cache.assemblies += 1
    return (;
        A = _from_local66(A66, R), B = _from_local66(B66, R),
        A_uncorrected = _from_local66(A66_E, R),
        B_uncorrected = _from_local66(B66_E, R),
        blocks, volume = blocks[1].V,
    )
end

"""
    _axi_run_cond(incl, K₀) -> NamedTuple

Full corrected solve in transport: modes 0 and 1, the fixed point, and the two
localization tensors reassembled in the global frame.
"""
function _axi_run_cond(incl::MeanFieldHom.FEExcenteredSphere, K₀::TensND.AbstractTens{2, 3})
    k₀ = _axi_iso_scalar(K₀)
    R = _axi_frame(incl)
    Dmap = Dict(
        AXI_SET_CORE => _axi_check_ti2(_to_local2(incl.props[1], R), "the core"),
        AXI_SET_SHELL => _axi_check_ti2(_to_local2(incl.props[2], R), "the shell"),
        AXI_SET_MATRIX => _to_local2(K₀, R),
    )
    grid = _axi_setup(incl).grid

    blocks = map((0, 1)) do m
        ms = _axi_mode_setup(incl, :conduction, m)
        A_E, B_E, A_p, B_p, V = _axi_solve_mode(
            ms, grid, Dmap, _axi_B_cond, _axi_project_cond, 1,
            _j -> (x -> (_axi_bc_affine_cond(m, x[1], x[2]),)),
            (_j, V) -> (x -> (_axi_bc_dipole_cond(m, x[1], x[2], k₀, V),)),
        )
        A, B, X = _axi_correct(A_E, B_E, A_p, B_p, fill(k₀, 1, 1))
        (; m, A_E, B_E, A_p, B_p, A, B, X, V)
    end

    diag3(t, a) = [t 0.0 0.0; 0.0 t 0.0; 0.0 0.0 a]
    A33 = diag3(blocks[2].A[1, 1], blocks[1].A[1, 1])
    B33 = diag3(blocks[2].B[1, 1], blocks[1].B[1, 1])
    A33_E = diag3(blocks[2].A_E[1, 1], blocks[1].A_E[1, 1])
    B33_E = diag3(blocks[2].B_E[1, 1], blocks[1].B_E[1, 1])
    incl.cache.assemblies += 1
    return (;
        A = _from_local33(A33, R), B = _from_local33(B33, R),
        A_uncorrected = _from_local33(A33_E, R),
        B_uncorrected = _from_local33(B33_E, R),
        blocks, volume = blocks[1].V,
    )
end

# ─── The seam ────────────────────────────────────────────────────────────────

_axi_cache_key(P₀::TensND.AbstractTens{4, 3}) =
    map(x -> round(x, sigdigits = 12), Tuple(MeanFieldHom.Core._C_array(P₀)))
_axi_cache_key(P₀::TensND.AbstractTens{2, 3}) =
    map(x -> round(x, sigdigits = 12), Tuple(TensND.components_canon(P₀)))

_axi_run(incl, P₀::TensND.AbstractTens{4, 3}) = _axi_run_elastic(incl, P₀)
_axi_run(incl, P₀::TensND.AbstractTens{2, 3}) = _axi_run_cond(incl, P₀)

function _fe_axi_localization(
        incl::MeanFieldHom.FEExcenteredSphere, P₀::TensND.AbstractTens; kw...
    )
    res = get!(incl.cache.tensors, _axi_cache_key(P₀)) do
        _axi_run(incl, P₀)
    end
    return res.A, res.B
end

"""
    fe_axi_breakdown(incl, P₀) -> NamedTuple

Diagnostic view of the corrected axisymmetric solve. Besides the corrected
localization tensors `A`, `B` it returns their **uncorrected** counterparts —
those of the plain truncated cell, `u|∂Ω = E·x` — and the per-mode blocks
`(A_E, B_E, A_p, B_p, A, B, X)`, plus the measured inclusion volume.

`A_uncorrected` drifts with `radius_ratio` like `(a/R)³` while `A` does not:
that contrast is the practical proof that the correction is wired correctly.
Bypasses the cache.
"""
FE.fe_axi_breakdown(incl::MeanFieldHom.FEExcenteredSphere, P₀::TensND.AbstractTens) =
    _axi_run(incl, P₀)

"""
    fe_axi_mesh_report(incl) -> NamedTuple

Mesh diagnostics: cell and node counts per region, and the measured volumes of
the core, the shell and the whole cell against their exact values. Builds the
grid if it does not exist yet, and caches it.
"""
function FE.fe_axi_mesh_report(incl::MeanFieldHom.FEExcenteredSphere)
    grid = _axi_setup(incl).grid
    ip_geo = Ferrite.Lagrange{Ferrite.RefTriangle, 1}()
    qr = Ferrite.QuadratureRule{Ferrite.RefTriangle}(3)
    cv = Ferrite.CellValues(qr, ip_geo, ip_geo)
    vol(set) = begin
        V = 0.0
        for ci in Ferrite.getcellset(grid, set)
            coords = Ferrite.getcoordinates(grid, ci)
            Ferrite.reinit!(cv, coords)
            for q in 1:Ferrite.getnquadpoints(cv)
                V += Ferrite.getdetJdV(cv, q) *
                    Ferrite.spatial_coordinate(cv, q, coords)[1]
            end
        end
        2π * V
    end
    a = Float64(incl.a)
    ac = Float64(FE.core_radius(incl))
    R = incl.mesh.radius_ratio * a
    return (;
        ncells = Ferrite.getncells(grid),
        nnodes = Ferrite.getnnodes(grid),
        ncells_core = length(Ferrite.getcellset(grid, AXI_SET_CORE)),
        ncells_shell = length(Ferrite.getcellset(grid, AXI_SET_SHELL)),
        ncells_matrix = length(Ferrite.getcellset(grid, AXI_SET_MATRIX)),
        volume_core = vol(AXI_SET_CORE),
        volume_core_exact = 4π * ac^3 / 3,
        volume_shell = vol(AXI_SET_SHELL),
        volume_shell_exact = 4π * (a^3 - ac^3) / 3,
        volume_cell = vol(AXI_SET_CORE) + vol(AXI_SET_SHELL) + vol(AXI_SET_MATRIX),
        volume_cell_exact = 4π * R^3 / 3,
    )
end
