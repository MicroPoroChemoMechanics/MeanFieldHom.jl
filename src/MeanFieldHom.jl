"""
    MeanFieldHom

Julia package for mean-field homogenisation of heterogeneous materials.

`MeanFieldHom` unifies the computation of Hill polarisation tensors for
ellipsoidal inhomogeneities, crack opening displacement (COD) tensors, stress
and displacement intensity factors, and — in the near future — full
homogenisation schemes, representative volume elements (RVEs), and
viscoelastic constitutive laws, sharing a common abstraction for
inclusions, algorithms, and material symmetry classes.

# Sub-modules

- `MeanFieldHom.Core`         — abstractions (`AbstractInclusion`,
  `AbstractAlgorithm`, `MaterialSymmetry`), shared numerics
  (Green / Newton kernels, Masson-style residue, DECUHR integrand), modulus
  extractors, and central dispatch.
- `MeanFieldHom.Elasticity`   — Hill polarisation for ellipsoidal inclusions
  (2D / 3D, isotropic and anisotropic matrix).
- `MeanFieldHom.Cracks`       — COD tensors, compliance contributions, SIF
  and DIF for elliptic and ribbon cracks.
- `MeanFieldHom.Conductivity` — 2nd-order Hill tensor for conductivity /
  diffusion problems.
- `MeanFieldHom.Schemes`      — placeholder for future mean-field
  homogenisation schemes (dilute, Mori–Tanaka, self-consistent, PCW, …).

# Shared generic interface

```julia
hill_tensor(ell::AbstractEllipsoidalInclusion, C₀; method=:auto, ...)
cod_tensor(crack::AbstractCrack, C₀; method=:auto, ...)
compliance_contribution(crack, C₀; method=:auto, ...)    # returns H (or R)
delta_compliance(crack, H, ε)                             # ΔS = factor · ε · H
delta_resistivity(crack, R, ε)                            # ΔR = factor · ε · R
sif(crack, C₀, Σ; method=:auto, ...)
dif(crack, C₀, Σ; method=:auto, ...)
```

All high-level entry points share the same algorithmic traits
(`Analytical`, `Residue`, `DECUHR`) and the same material-symmetry dispatch
rules. See the developer documentation (`docs/src/developer/`) for guidance
on extending the package with new inclusions, algorithms or schemes.
"""
module MeanFieldHom

using TensND

include("Elliptic/Elliptic.jl")
include("Core/Core.jl")
include("Elasticity/Elasticity.jl")
include("Cracks/Cracks.jl")
include("Conductivity/Conductivity.jl")
include("LayeredSpheres/LayeredSpheres.jl")
include("Schemes/Schemes.jl")

using .Elliptic
using .Core
using .Elasticity
using .Cracks
using .Conductivity
using .LayeredSpheres
using .Schemes

# ─── Localization + contribution (top-level: need all sub-module APIs) ──────
# Generics are declared in Core; Cracks already defines `compliance_contribution`,
# `delta_compliance`, `delta_resistivity`.  Extend both via qualified imports so
# that every method attaches to the same canonical function.
import .Core: strain_strain_loc, stress_strain_loc, strain_stress_loc,
    stress_stress_loc, gradient_gradient_loc, flux_gradient_loc,
    gradient_flux_loc, flux_flux_loc,
    stiffness_contribution, conductivity_contribution,
    resistivity_contribution, delta_stiffness, delta_conductivity
import .Cracks: compliance_contribution, delta_compliance, delta_resistivity

include("localization.jl")
include("contribution.jl")

# ── Abstractions ─────────────────────────────────────────────────────────────
export AbstractInclusion, AbstractEllipsoidalInclusion, AbstractCrack
export AbstractLayeredInclusion
export AbstractAlgorithm, Analytical, Residue, DECUHR, NestedQuadGK,
    CylinderQuadrature, Auto
export MaterialSymmetry, IsotropicSym, TransverselyIsotropicSym,
    OrthotropicSym, GeneralAnisotropicSym
export material_symmetry, dimension, inclusion_basis, shape_trait, shape_tensor
export eshelby_tensor

# ── Elasticity ───────────────────────────────────────────────────────────────
export Ellipsoid
export EllipsoidShape, Spherical, Prolate, Oblate, Triaxial, Circular, Elliptic
export Cylinder, CylindricalShape, CircularCylindrical, EllipticCylindrical
export newton_potential_3d_cylinder
export tens_IA, tens_UA, tens_VA
export hill_tensor

# ── Cracks ───────────────────────────────────────────────────────────────────
export CrackShape, Penny, EllipticShape, Ribbon
export EllipticCrack, RibbonCrack, PennyCrack
export crack_basis, aspect_ratio, semi_major, semi_minor, crack_normal
export cod_tensor, B_tensor
export cod_from_compliance, compliance_from_cod
export compliance_contribution
export delta_compliance, delta_resistivity
export sif, dif

# ── Localization & contribution (Eshelby dilute, Kachanov-Sevostianov) ───────
export strain_strain_loc, stress_strain_loc, strain_stress_loc, stress_stress_loc
export gradient_gradient_loc, flux_gradient_loc, gradient_flux_loc, flux_flux_loc
export stiffness_contribution, conductivity_contribution, resistivity_contribution
export delta_stiffness, delta_conductivity

# ── LayeredSphere (Hervé-Zaoui / Hervé-Luanco / Gurtin-Murdoch / Kapitza) ────
export LayeredSphere, AbstractInterface, PerfectInterface
export SpringInterface, MembraneInterface
export KapitzaInterface, SurfaceConductiveInterface
export layer_count, layer_radius, layer_modulus, layer_interface,
    layer_volume_fraction, outer_radius
export layer_strain_average, sphere_strain_average, cumulative_strain_average

# ── Elliptic integrals (type-generic) ────────────────────────────────────────
export ell_K, ell_E, ell_F, ell_RF, ell_RD

# ── Backwards-compat aliases ─────────────────────────────────────────────────
const HillAlgorithm = AbstractAlgorithm
const CrackAlgorithm = AbstractAlgorithm
export HillAlgorithm, CrackAlgorithm

end # module
