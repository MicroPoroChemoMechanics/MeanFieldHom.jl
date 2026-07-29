# =============================================================================
#  solver.jl — the corrected finite Eshelby cell for a flat crack.
#
#  Weak form: pure linear elasticity, `∫ σ(u):ε(v) dΩ = 0`, no body force.  The
#  crack is a **zero-thickness discontinuity** (duplicated nodes) whose lips are
#  traction-free *naturally* — there is no interface term, no multiplier, no
#  contact condition.  The only Dirichlet condition is on the outer sphere, and
#  its value is where the whole method lives.
#
#  Per evaluation (Adessina et al. 2017, crack declination — 3 + 3 instead of
#  the general 6 + 6):
#
#    1. one assembly, one Cholesky factorization of the free-free block;
#    2. three "traction" solves,  u|∂Ω = (𝕊₀:Σ⁽ⁱ⁾)·x  with Σ⁽ⁱ⁾·n̂ = eᵢ   → 𝐁ₛ
#    3. three "dipole"   solves,  u|∂Ω = -∇G(x):Πₘ,  Πₘ = b·S_f·ℂ₀:(eₘ⊗ˢn̂) → 𝐁ᵤ
#
#  The minus in step 3 is *not* the one in `axi_solver.jl`, which imposes
#  `+∇G:M`.  The universal rule is `u = +∇G:M` with `M` the polarization
#  moment; the moment of a displacement discontinuity is `M = -S_f ℂ₀:(⟨[[u]]⟩⊗ˢn̂)`,
#  and `Πₘ` above is its opposite, so the sign is carried here instead.  The
#  two files agree; they just put the minus in different places.
#    4. 𝐁∞ = (𝟏 - 𝐁ᵤ)⁻¹·𝐁ₛ
#
#  Step 3 is the correction: 𝐁ₛ is the COD of the *finite* cell, 𝐁ᵤ its
#  response to the crack's own far field, and step 4 solves the resulting
#  linear fixed point in closed form.
#
#  The COD itself is measured as a surface integral of the jump over each lip,
#  normalized by the semi-minor axis `b` — the convention of `cod_tensor`:
#
#      ⟨[[u]]⟩ / b = 𝐁 · (Σ·n̂).
# =============================================================================

"Discretization built once per crack geometry and reused for every `C₀`."
struct FESetup{G, D, C, F}
    grid::G
    dh::D
    cv::C
    fv::F
    lip_up::Set{Ferrite.FacetIndex}
    lip_dn::Set{Ferrite.FacetIndex}
    presc::Vector{Int}
    free::Vector{Int}
end

function _build_setup(crack::MeanFieldHom.FEEllipticCrack)
    a, b = Float64(crack.a), Float64(crack.b)
    opts = crack.mesh
    R = opts.radius_ratio * a

    gmsh.initialize()
    local grid0
    try
        gmsh.option.setNumber("General.Terminal", 0)
        _build_gmsh_crack_model(a, b, R, opts.htipdiv)
        grid0 = FerriteGmsh.togrid()
    finally
        gmsh.finalize()
    end

    grid, nweld = weld_crack_front(grid0, a, b)
    nweld == 0 && @warn "no crack-front node pair was welded — the crack front " *
        "may be split, which overestimates the opening" a b
    lip_up, lip_dn = split_crack_lips(grid)

    # Straight tetrahedra (P1 geometry) with a P1 or P2 displacement field.
    # Curving the geometry is *not* an option here: `setOrder(2)` runs after the
    # `Crack` plugin, and it curves the two lips' front edges differently, so
    # the welded front would come apart again at the mid-side nodes.
    ip_geo = Ferrite.Lagrange{Ferrite.RefTetrahedron, 1}()
    ip = Ferrite.Lagrange{Ferrite.RefTetrahedron, opts.order}()^3
    qr = Ferrite.QuadratureRule{Ferrite.RefTetrahedron}(2 * opts.order)
    fqr = Ferrite.FacetQuadratureRule{Ferrite.RefTetrahedron}(2 * opts.order)
    cv = Ferrite.CellValues(qr, ip, ip_geo)
    fv = Ferrite.FacetValues(fqr, ip, ip_geo)

    dh = Ferrite.DofHandler(grid)
    Ferrite.add!(dh, :u, ip)
    Ferrite.close!(dh)

    ch = _constraint_handler(dh, grid, _ -> Ferrite.Vec{3}((0.0, 0.0, 0.0)))
    presc = collect(ch.prescribed_dofs)
    free = setdiff(1:Ferrite.ndofs(dh), presc)

    return FESetup(grid, dh, cv, fv, lip_up, lip_dn, presc, free)
end

_setup(crack::MeanFieldHom.FEEllipticCrack) =
    crack.cache.setup === nothing ?
    (crack.cache.setup = _build_setup(crack)) : crack.cache.setup

function _constraint_handler(dh, grid, f)
    ch = Ferrite.ConstraintHandler(dh)
    Ferrite.add!(ch, Ferrite.Dirichlet(:u, Ferrite.getfacetset(grid, SET_OUTER), (x, _t) -> f(x)))
    Ferrite.close!(ch)
    Ferrite.update!(ch, 0.0)
    return ch
end

function _assemble(s::FESetup, C::Tensors.SymmetricTensor{4, 3, Float64})
    K = Ferrite.allocate_matrix(s.dh)
    asm = Ferrite.start_assemble(K)
    n = Ferrite.getnbasefunctions(s.cv)
    ke = zeros(n, n)
    for cell in Ferrite.CellIterator(s.dh)
        Ferrite.reinit!(s.cv, cell)
        fill!(ke, 0)
        for q in 1:Ferrite.getnquadpoints(s.cv)
            dΩ = Ferrite.getdetJdV(s.cv, q)
            for i in 1:n
                σi = C ⊡ Ferrite.shape_symmetric_gradient(s.cv, q, i)
                for j in i:n
                    ke[i, j] += (σi ⊡ Ferrite.shape_symmetric_gradient(s.cv, q, j)) * dΩ
                end
            end
        end
        for i in 1:n, j in 1:(i - 1)
            ke[i, j] = ke[j, i]
        end
        Ferrite.assemble!(asm, Ferrite.celldofs(cell), ke)
    end
    return K
end

"""
    _mean_jump(s, u, S_f, b) -> Vector{3}

`⟨[[u]]⟩ / b`, the crack opening averaged over the crack surface and
normalized by the semi-minor axis, from a surface integral of the trace of `u`
on each lip — no assumption on the opening profile.
"""
function _mean_jump(s::FESetup, u::Vector{Float64}, S_f::Float64, b::Float64)
    acc = zeros(3)
    for (lip, sgn) in ((s.lip_up, 1.0), (s.lip_dn, -1.0))
        for fi in lip
            Ferrite.reinit!(s.fv, Ferrite.getcoordinates(s.grid, fi[1]), fi[2])
            ue = u[Ferrite.celldofs(s.dh, fi[1])]
            for q in 1:Ferrite.getnquadpoints(s.fv)
                uq = Ferrite.function_value(s.fv, q, ue)
                dS = Ferrite.getdetJdV(s.fv, q)
                acc[1] += sgn * uq[1] * dS
                acc[2] += sgn * uq[2] * dS
                acc[3] += sgn * uq[3] * dS
            end
        end
    end
    return acc ./ (S_f * b)
end

# Σ⁽ⁱ⁾ with Σ⁽ⁱ⁾·e₃ = eᵢ, and eₘ ⊗ˢ e₃ — in the *local* crack frame, where the
# normal is e₃ by construction of the mesh.
const _SIGMA_T = (
    Tensors.SymmetricTensor{2, 3}([0.0 0.0 1.0; 0.0 0.0 0.0; 1.0 0.0 0.0]),
    Tensors.SymmetricTensor{2, 3}([0.0 0.0 0.0; 0.0 0.0 1.0; 0.0 1.0 0.0]),
    Tensors.SymmetricTensor{2, 3}([0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 1.0]),
)
const _EM_N = (
    Tensors.SymmetricTensor{2, 3}([0.0 0.0 0.5; 0.0 0.0 0.0; 0.5 0.0 0.0]),
    Tensors.SymmetricTensor{2, 3}([0.0 0.0 0.0; 0.0 0.0 0.5; 0.0 0.5 0.0]),
    Tensors.SymmetricTensor{2, 3}([0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 1.0]),
)

"""
    _solve_cod_local(crack, C_loc, μ, ν) -> (B_s, B_u, B_inf)

The 3 + 3 corrected scheme, entirely in the crack's local frame. Returns plain
`3×3` matrices.
"""
function _solve_cod_local(
        crack::MeanFieldHom.FEEllipticCrack,
        C::Tensors.SymmetricTensor{4, 3, Float64}, μ::Float64, ν::Float64
    )
    s = _setup(crack)
    a, b = Float64(crack.a), Float64(crack.b)
    S_f = π * a * b

    K = _assemble(s, C)
    crack.cache.assemblies += 1

    # One factorization of the free-free block serves all six right-hand sides.
    # (Assembling once and eliminating the prescribed dofs by hand, rather than
    # `apply!`-ing the constraints, which would consume the matrix.)
    F = LinearAlgebra.cholesky(LinearAlgebra.Symmetric(K[s.free, s.free]))
    Kfp = K[s.free, s.presc]
    nd = Ferrite.ndofs(s.dh)

    function solve_with(f)
        ch = _constraint_handler(s.dh, s.grid, f)
        g = ch.inhomogeneities
        u = zeros(nd)
        u[s.presc] .= g
        u[s.free] .= F \ Vector(-Kfp * g)
        return u
    end

    S₀ = inv(C)

    B_s = zeros(3, 3)
    for i in 1:3
        E = S₀ ⊡ _SIGMA_T[i]
        B_s[:, i] .= _mean_jump(s, solve_with(x -> E ⋅ x), S_f, b)
    end

    B_u = zeros(3, 3)
    for m in 1:3
        Π = Matrix((b * S_f) * (C ⊡ _EM_N[m]))
        B_u[:, m] .= _mean_jump(
            s, solve_with(
                x -> Ferrite.Vec{3}(
                    Tuple(-MeanFieldHom.Core._dipole_displacement_iso(μ, ν, x, Π))
                )
            ), S_f, b
        )
    end

    B_inf = (LinearAlgebra.I - B_u) \ B_s
    return B_s, B_u, B_inf
end
