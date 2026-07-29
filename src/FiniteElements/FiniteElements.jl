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

Only the **types** and their geometric interface live here; mesh generation and
the solves live in the package extension `MeanFieldHomFerriteExt`, loaded when
`Ferrite`, `FerriteGmsh` and `Gmsh` are available. Without them the seams below
error informatively.

See `docs/src/manual/fe_inclusions.md` and
`docs/src/applications/recycled_aggregate.md`.
"""
module FiniteElements

using TensND

import ..Core
import ..Cracks
import ..Schemes

export FECache, fe_assembly_count, fe_reset!
export FEMeshOptions, FEEllipticCrack, fe_cod_breakdown, fe_mesh_report
export FEAxiMeshOptions, FEExcenteredSphere
export fe_axi_breakdown, fe_axi_mesh_report, fe_axi_localization

include("common.jl")
include("crack.jl")
include("excentered_sphere.jl")

end # module
