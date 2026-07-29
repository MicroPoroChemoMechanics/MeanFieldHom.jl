"""
    MeanFieldHomFerriteExt

Finite-element resolution of the Eshelby problem, activated by
`import Ferrite, FerriteGmsh, Gmsh`.

Implements [`MeanFieldHom.FEEllipticCrack`](@ref): a flat elliptical crack
whose crack-opening-displacement tensor comes out of a finite-element solve
instead of a closed form, using the finite Eshelby cell with a first-order
corrected boundary condition of Adessina, Barthélémy, Lavergne & Ben Fraj,
*Int. J. Eng. Sci.* **119** (2017) 1-15.

The point of the exercise is that nothing downstream knows: because the type
subtypes `AbstractCrack` and declares a standard `shape_trait`, supplying
`cod_tensor` is enough for ℍ, ℕ, 𝐑, 𝐍_K, the bundled pair and the four
`delta_*` to be inherited, and the crack drops into every scheme —
`symmetrize` included.
"""
module MeanFieldHomFerriteExt

using MeanFieldHom
using TensND
using Ferrite
using FerriteGmsh
using Gmsh: gmsh
using Tensors
using LinearAlgebra

const FE = MeanFieldHom.FiniteElements

import MeanFieldHom.FiniteElements: _fe_cod_tensor, _fe_axi_localization

include("mesh.jl")
include("solver.jl")
include("axi_mesh.jl")
include("axi_fourier.jl")
include("axi_solver.jl")

# ─── Frame handling ──────────────────────────────────────────────────────────
#
#  The mesh is built once, in the crack's *local* frame (crack in the z = 0
#  plane, semi-axes along x and y, normal along z).  The reference medium is
#  therefore rotated into that frame before the solve and the resulting COD
#  tensor rotated back afterwards.  Besides being necessary, this is what makes
#  the memoization pay off under orientation averaging: a whole family of
#  identically-shaped cracks at different orientations, embedded in the same
#  (isotropic, or symmetrization-projected) matrix, share one single solve.

"Rotation matrix whose columns are the crack's local axes in global coordinates."
function _frame_matrix(crack)
    l̂, m̂, n̂ = MeanFieldHom.Core._frame_columns(MeanFieldHom.inclusion_basis(crack))
    return hcat(l̂, m̂, n̂)
end

"Components of `C₀` in the crack's local frame, as a `SymmetricTensor{4,3}`."
function _local_stiffness(crack, C₀::TensND.AbstractTens{4, 3})
    Cg = MeanFieldHom.Core._C_array(C₀)
    R = _frame_matrix(crack)
    Cl = zeros(Float64, 3, 3, 3, 3)
    @inbounds for p in 1:3, q in 1:3, r in 1:3, s in 1:3
        acc = 0.0
        for i in 1:3, j in 1:3, k in 1:3, l in 1:3
            acc += R[i, p] * R[j, q] * R[k, r] * R[l, s] * Cg[i, j, k, l]
        end
        Cl[p, q, r, s] = acc
    end
    return Tensors.SymmetricTensor{4, 3}((i, j, k, l) -> Cl[i, j, k, l])
end

"Cache key: the local-frame stiffness, rounded so that arithmetic noise between
two otherwise identical scheme iterations does not miss the cache."
_cache_key(C::Tensors.SymmetricTensor{4, 3, Float64}) =
    map(x -> round(x, sigdigits = 12), Tuple(C))

"""
    _isotropic_moduli(crack, C₀) -> (μ, ν)

Shear modulus and Poisson ratio of the reference medium, refusing anything
that is not isotropic.

The test is on the **content**, not on the TensND type: a self-consistent or
differential iterate arrives as a `TensCanonical` even when its content is
isotropic, and rejecting it on the type alone would rule out perfectly valid
uses (an isotropically-symmetrized crack family, for instance).
"""
function _isotropic_moduli(crack, C₀::TensND.AbstractTens{4, 3}; rtol = 1.0e-8)
    C_iso = MeanFieldHom.Core.isotropify(C₀)
    if !(C₀ isa TensND.TensISO{4, 3})
        A, Aiso = MeanFieldHom.Core._C_array(C₀), MeanFieldHom.Core._C_array(C_iso)
        dev = maximum(abs, A .- Aiso)
        dev ≤ rtol * max(maximum(abs, Aiso), eps()) || throw(
            ArgumentError(
                "`FEEllipticCrack` supports an isotropic reference medium only — " *
                    "the corrected boundary condition uses the closed-form Kelvin " *
                    "dipole field, whose anisotropic counterpart (Pan-Chou / " *
                    "Barnett-Willis) is not implemented. The reference medium " *
                    "deviates from isotropy by $(round(dev / maximum(abs, Aiso) * 100, sigdigits = 3)) %.\n" *
                    "This is the normal situation for an iterative scheme " *
                    "(`SelfConsistent`, `DifferentialScheme`) on a *parallel* crack " *
                    "family, whose effective medium is genuinely anisotropic. Add " *
                    "`symmetrize = IsoSymmetrize()` to the phase — the scheme then " *
                    "hands the kernel an isotropic reference — or use the " *
                    "closed-form `EllipticCrack`."
            )
        )
    end
    E, ν = MeanFieldHom.Core.extract_iso_moduli(C_iso)
    return Float64(E / (2 * (1 + ν))), Float64(ν)
end

# ─── The seam ────────────────────────────────────────────────────────────────

function _fe_cod_tensor(
        crack::MeanFieldHom.FEEllipticCrack, C₀::TensND.AbstractTens{4, 3}; kw...
    )
    μ, ν = _isotropic_moduli(crack, C₀)
    C_loc = _local_stiffness(crack, C₀)
    key = _cache_key(C_loc)

    B_loc = get!(crack.cache.tensors, key) do
        _, _, B = _solve_cod_local(crack, C_loc, μ, ν)
        B
    end

    R = _frame_matrix(crack)
    return TensND.TensCanonical(R * B_loc * R')
end

"""
    fe_cod_breakdown(crack, C₀) -> (; B_s, B_u, B_inf, B_s_glob, B_inf_glob)

Diagnostic view of the corrected solve: the COD tensor of the **finite** cell
`B_s`, the response `B_u` to the crack's own dipole far field, and the
infinite-medium result `B_inf = (1 - B_u)⁻¹ B_s`, all in the crack's local
frame, plus `B_s` and `B_inf` rotated back to the global frame.

`norm(B_u)` measures how much work the boundary correction is doing; it should
fall like `(a/R)³`, and `B_inf` — unlike `B_s` — should be insensitive to
`radius_ratio`. That contrast is the practical proof that the correction is
wired correctly.

Bypasses the cache.
"""
function FE.fe_cod_breakdown(
        crack::MeanFieldHom.FEEllipticCrack, C₀::TensND.AbstractTens{4, 3}
    )

    μ, ν = _isotropic_moduli(crack, C₀)
    C_loc = _local_stiffness(crack, C₀)
    B_s, B_u, B_inf = _solve_cod_local(crack, C_loc, μ, ν)
    R = _frame_matrix(crack)
    return (;
        B_s, B_u, B_inf,
        B_s_glob = TensND.TensCanonical(R * B_s * R'),
        B_inf_glob = TensND.TensCanonical(R * B_inf * R'),
    )
end

"""
    fe_mesh_report(crack) -> NamedTuple

Mesh diagnostics: cell and dof counts, number of welded crack-front node pairs,
the two lip facet counts and their measured areas against the exact `πab`.
Builds the discretization if it does not exist yet, and caches it.
"""
function FE.fe_mesh_report(crack::MeanFieldHom.FEEllipticCrack)
    s = _setup(crack)
    a, b = Float64(crack.a), Float64(crack.b)
    ip_geo = Ferrite.Lagrange{Ferrite.RefTetrahedron, 1}()
    qr = Ferrite.FacetQuadratureRule{Ferrite.RefTetrahedron}(2)
    area_up = facetset_area(s.grid, s.lip_up, ip_geo, qr)
    area_dn = facetset_area(s.grid, s.lip_dn, ip_geo, qr)
    return (;
        ncells = Ferrite.getncells(s.grid),
        nnodes = Ferrite.getnnodes(s.grid),
        ndofs = Ferrite.ndofs(s.dh),
        nfacets_up = length(s.lip_up),
        nfacets_dn = length(s.lip_dn),
        area_up, area_dn,
        area_exact = π * a * b,
    )
end

end # module
