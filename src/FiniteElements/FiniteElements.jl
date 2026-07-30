"""
    MeanFieldHom.FiniteElements

Inclusions whose response is obtained from a **finite-element resolution of the
Eshelby problem** rather than from a closed form — the package's own
demonstration that the [`CustomInclusions`](@ref MeanFieldHom.CustomInclusions)
contract is enough to reach every scheme.

| Type | Morphology | Discretization | Entry gate |
|---|---|---|---|
| [`FEEllipticCrack`](@ref) | flat elliptical crack | 3-D tetrahedra | COD tensor → crack algebra |
| [`FEExcenteredSphere`](@ref) | sphere with an off-centre spherical core | axisymmetric Fourier | B — the two localization tensors |

Both use the *finite Eshelby cell with a first-order corrected boundary
condition* of Adessina, Barthélémy, Lavergne & Ben Fraj, *Int. J. Eng. Sci.*
**119** (2017) 1-15: the infinite matrix is truncated to a ball of finite
radius and the truncation bias is removed by adding the inclusion's own dipole
far field to the imposed boundary displacement.

The types, the Fourier operators, the boundary data and the algebra of the
corrected boundary condition all live here. What a package extension supplies
is only the **discretization** — a mesh, scalar Lagrange spaces, an assembly
and a quadrature — through the nine generics of [`FEBackend`](@ref).

Two backends exist. [`FerriteBackend`](@ref) (`import Ferrite, FerriteGmsh,
Gmsh`) serves both inclusions; [`GridapBackend`](@ref) (`import Gridap,
GridapGmsh`) serves the axisymmetric one. An inclusion built without naming a
backend takes [`AutoBackend`](@ref) and picks one at its first solve; with
neither loaded, that solve errors informatively.

See `docs/src/manual/fe_inclusions.md` and
`docs/src/applications/recycled_aggregate.md`.
"""
module FiniteElements

using TensND

import LinearAlgebra
import Tensors

import ..Core
import ..Cracks
import ..Schemes

export FECache, fe_assembly_count, fe_reset!
export FEBackend, AutoBackend, FerriteBackend, GridapBackend
export FEMeshOptions, FEEllipticCrack, fe_cod_breakdown, fe_mesh_report
export FEAxiMeshOptions, FEExcenteredSphere
export fe_axi_breakdown, fe_axi_mesh_report, fe_axi_localization

include("common.jl")
include("backends.jl")
include("crack.jl")
include("excentered_sphere.jl")

# The axisymmetric solver: geometry, Fourier operators, algebra, driver.
include("axi_gmsh_geometry.jl")
include("axi_fourier.jl")
include("axi_algebra.jl")
include("axi_driver.jl")

# ─── Sensitivity is not available through a finite-element geometry ──────────
#
#  `Schemes._replace_geom_field` rebuilds a geometry by copying every
#  non-`<:Number` field by reference — which for these types includes the
#  `FECache`.  The perturbed inclusion would therefore share the original's
#  memoized tensors, whose key is the reference medium alone, and hand back the
#  *unperturbed* answer; the derivative would come out as exactly zero, quietly.
#  The grid stored in that same cache was built for the original geometry too.
#
#  Even without the cache the answer would be zero: the solve converts its
#  geometry to `Float64` on entry, so a `ForwardDiff.Dual` loses its
#  perturbation at the door.  Refusing is the only honest option.

const _FEGeometry = Union{FEEllipticCrack, FEExcenteredSphere}

_no_fe_sensitivity(geom, name) = error(
    "analytic sensitivity is not available through `$(nameof(typeof(geom)))`: " *
        "the finite-element solve runs in `Float64` and memoizes on the reference " *
        "medium alone, so differentiating `:$name` would silently return zero. " *
        "Use a finite difference over freshly constructed inclusions instead."
)

Schemes._replace_geom_field(geom::_FEGeometry, ::Val{name}, ::Nothing, value) where {name} =
    _no_fe_sensitivity(geom, name)
Schemes._replace_geom_field(geom::_FEGeometry, ::Val{name}, ::Int, value) where {name} =
    _no_fe_sensitivity(geom, name)

end # module
