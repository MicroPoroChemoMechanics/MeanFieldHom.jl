# =============================================================================
#  crack_backend.jl — the Ferrite implementation of the flat-crack contract.
#
#  Seven methods.  The mesh comes from the shared gmsh geometry, already welded
#  along the crack front; what is left here is the grid, the space, the dof
#  split, the Dirichlet lift, the assembly and the lip integral.
# =============================================================================

"Grid, space and the pieces the seam needs, built once per crack geometry."
struct FerriteCrackSpace{D, C, F}
    dh::D
    cv::C
    fv::F
    grid::Any
    lip_up::Set{Ferrite.FacetIndex}
    lip_dn::Set{Ferrite.FacetIndex}
    ch::Ferrite.ConstraintHandler
    bcref::Base.RefValue{Any}
    presc::Vector{Int}
    free::Vector{Int}
end

"""
    split_crack_lips(grid) -> (lip_up, lip_dn)

Split the crack facet set in two by the side of the plane its element sits on.
The element above the crack carries the trace `u(0⁺)`, hence the `+` lip.
"""
function split_crack_lips(grid::Ferrite.Grid)
    up, dn = Set{Ferrite.FacetIndex}(), Set{Ferrite.FacetIndex}()
    for fi in Ferrite.getfacetset(grid, FE.SET_CRACK)
        zc = sum(x -> x[3], Ferrite.getcoordinates(grid, fi[1])) / 4
        push!(zc > 0 ? up : dn, fi)
    end
    return up, dn
end

"Area of a facet set, on the geometric interpolation."
function facetset_area(grid::Ferrite.Grid, set, ip_geo, qr)
    fv = Ferrite.FacetValues(qr, ip_geo, ip_geo)
    A = 0.0
    for fi in set
        Ferrite.reinit!(fv, Ferrite.getcoordinates(grid, fi[1]), fi[2])
        for q in 1:Ferrite.getnquadpoints(fv)
            A += Ferrite.getdetJdV(fv, q)
        end
    end
    return A
end

# ─── Mesh ────────────────────────────────────────────────────────────────────

function FE.fe_crack_grid(::FE.FerriteBackend, crack::MeanFieldHom.FEEllipticCrack)
    a, b = Float64(crack.a), Float64(crack.b)
    opts = crack.mesh
    path = joinpath(mktempdir(), "mfh_crack.msh")
    gmsh.initialize()
    try
        gmsh.option.setNumber("General.Terminal", 0)
        FE._build_gmsh_crack_model(gmsh, a, b, opts.radius_ratio * a, opts.htipdiv)
        gmsh.write(path)
    finally
        gmsh.finalize()
    end
    nweld = FE._weld_msh_crack_front(path, a, b)
    nweld == 0 && @warn "no crack-front node pair was welded — the crack front " *
        "may be split, which overestimates the opening" a b
    return redirect_stdout(devnull) do
        FerriteGmsh.togrid(path)
    end
end

function FE.fe_crack_counts(::FE.FerriteBackend, grid)
    up, dn = split_crack_lips(grid)
    ip_geo = Ferrite.Lagrange{Ferrite.RefTetrahedron, 1}()
    qr = Ferrite.FacetQuadratureRule{Ferrite.RefTetrahedron}(2)
    return (;
        ncells = Ferrite.getncells(grid),
        nnodes = Ferrite.getnnodes(grid),
        nfacets_up = length(up),
        nfacets_dn = length(dn),
        area_up = facetset_area(grid, up, ip_geo, qr),
        area_dn = facetset_area(grid, dn, ip_geo, qr),
    )
end

# ─── Space ───────────────────────────────────────────────────────────────────

function FE.fe_crack_space(::FE.FerriteBackend, grid, order::Int)
    # Straight tetrahedra (P1 geometry) with a P1 or P2 displacement field.
    # Curving the geometry is *not* an option: `setOrder(2)` would run after the
    # `Crack` plugin and curve the two lips' front edges differently, pulling
    # the welded front apart again at the mid-side nodes.
    ip_geo = Ferrite.Lagrange{Ferrite.RefTetrahedron, 1}()
    ip = Ferrite.Lagrange{Ferrite.RefTetrahedron, order}()^3
    qr = Ferrite.QuadratureRule{Ferrite.RefTetrahedron}(2 * order)
    fqr = Ferrite.FacetQuadratureRule{Ferrite.RefTetrahedron}(2 * order)
    cv = Ferrite.CellValues(qr, ip, ip_geo)
    fv = Ferrite.FacetValues(fqr, ip, ip_geo)

    dh = Ferrite.DofHandler(grid)
    Ferrite.add!(dh, :u, ip)
    Ferrite.close!(dh)

    bcref = Base.RefValue{Any}(_ -> (0.0, 0.0, 0.0))
    ch = Ferrite.ConstraintHandler(dh)
    Ferrite.add!(
        ch, Ferrite.Dirichlet(
            :u, Ferrite.getfacetset(grid, FE.SET_OUTER),
            (x, _t) -> Ferrite.Vec{3}(Tuple(bcref[](x)))
        )
    )
    Ferrite.close!(ch)

    presc = collect(ch.prescribed_dofs)
    free = setdiff(1:Ferrite.ndofs(dh), presc)
    up, dn = split_crack_lips(grid)
    return FerriteCrackSpace(dh, cv, fv, grid, up, dn, ch, bcref, presc, free)
end

FE.fe_crack_dof_split(::FE.FerriteBackend, s::FerriteCrackSpace) =
    (Ferrite.ndofs(s.dh), s.free, s.presc)

function FE.fe_crack_set_dirichlet!(::FE.FerriteBackend, s::FerriteCrackSpace, u, f)
    s.bcref[] = f
    Ferrite.update!(s.ch, 0.0)
    u[s.ch.prescribed_dofs] .= s.ch.inhomogeneities
    return u
end

# ─── Operators ───────────────────────────────────────────────────────────────

function FE.fe_crack_stiffness(
        ::FE.FerriteBackend, s::FerriteCrackSpace,
        C::Tensors.SymmetricTensor{4, 3, Float64}
    )
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

function FE.fe_crack_mean_jump(
        ::FE.FerriteBackend, s::FerriteCrackSpace, u::Vector{Float64},
        S_f::Float64, b::Float64
    )
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
