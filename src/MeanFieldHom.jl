"""
    MeanFieldHom

Julia package for mean-field homogenization of heterogeneous materials.

`MeanFieldHom` unifies the computation of Hill polarisation tensors for
ellipsoidal inhomogeneities, crack opening displacement (COD) tensors, stress
and displacement intensity factors, and — in the near future — full
homogenization schemes, representative volume elements (RVEs), and
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
  homogenization schemes (dilute, Mori–Tanaka, self-consistent, PCW, …).

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
include("Viscoelasticity/Viscoelasticity.jl")

using .Elliptic
using .Core
using .Elasticity
using .Cracks
using .Conductivity
using .LayeredSpheres
using .Schemes
using .Viscoelasticity

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
export Ellipsoid, Spheroid
export EllipsoidShape, Spherical, Prolate, Oblate, Triaxial, Circular, Elliptic
export Cylinder, CylindricalShape, CircularCylindrical, EllipticCylindrical
export newton_potential_3d_cylinder
export tens_IA, tens_UA, tens_VA
export hill_tensor
export k_mu, iso_stiffness, E_nu, iso_stiffness_E_nu
export hoenig_params, hoenig_stiffness

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
export is_homogeneous_inclusion
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

# ── Schemes : RVE + amounts + distribution shape + symmetrize ───────────────
export AbstractAmount, VolumeFraction, CrackDensity
export AbstractDistributionShape, UniformDistribution
export AbstractSymmetrize, NoSymmetrize, IsoSymmetrize, TISymmetrize
export isotropify, transverse_isotropify
export ti_average_mandel66, iso_average_mandel66
export best_fit_ti, best_fit_iso, best_fit_ortho
export polar_orientation_bins
export Phase, RVE
export add_matrix!, add_phase!
export matrix_phase, inclusion_phase_names
export phase_property, matrix_property
export volume_fraction, crack_density, matrix_volume_fraction
export phase_symmetrize
export validate_rve

# ── Schemes : scheme types + entry point ─────────────────────────────────────
export HomogenizationScheme
export Voigt, Reuss, Dilute, DiluteDual, MoriTanaka, Maxwell, PonteCastanedaWillis
export SelfConsistent, AsymmetricSelfConsistent
export AndersonDefault, NewtonDefault
export DifferentialTrajectory, Proportional, Sequential, CustomPath, Path, DifferentialScheme
export homogenize

# ── Schemes : sensitivities (autodiff via ForwardDiff strong dependency) ────
export AbstractParameter, AmountParameter, PropertyParameter,
    GeometryParameter, DistributionShapeParameter
export amount, property, geometry, shape_param
export get_param, set_param
export derivative, gradient, jacobian, sensitivity

# ── Viscoelasticity (ALV) ────────────────────────────────────────────────────
export AbstractViscoLaw, ViscoLaw, VALID_VISCO_MODES
export visco_mode, visco_eval
export maxwell_relaxation, kelvin_creep, maxwell_iso, kelvin_iso, heaviside_law
export trapezoidal_matrix
export volterra_inverse, volterra_product, volterra_divide, volterra_left_divide
export iso_params_from_blocks, iso_blocks_from_params
export ti_params_from_blocks, ti_blocks_from_params
export ortho_params_from_blocks, ortho_blocks_from_params
export AbstractALVKernel, ALVKernelISO, ALVKernelTI, ALVKernelOrtho
export hill_kernel
export dilute_concentration_alv, dilute_contribution_alv
export voigt_alv, reuss_alv, dilute_alv, dilute_dual_alv
export mori_tanaka_alv, maxwell_alv
export voigt_alv_iso, reuss_alv_iso, dilute_alv_iso, dilute_dual_alv_iso
export mori_tanaka_alv_iso, maxwell_alv_iso
export dilute_concentration_alv_iso, dilute_contribution_alv_iso
export voigt_alv_ti, reuss_alv_ti, dilute_alv_ti, dilute_dual_alv_ti
export mori_tanaka_alv_ti, maxwell_alv_ti
export dilute_concentration_alv_ti, dilute_contribution_alv_ti
export voigt_alv_ortho, reuss_alv_ortho, dilute_alv_ortho, dilute_dual_alv_ortho
export mori_tanaka_alv_ortho, maxwell_alv_ortho
export dilute_concentration_alv_ortho, dilute_contribution_alv_ortho
export self_consistent_alv, asymmetric_self_consistent_alv,
    pcw_alv, differential_alv
export bulk_localization_alv, bulk_state_seq_alv, shear_localization_alv
export strain_strain_loc_alv, stiffness_contribution_alv
export homogenize_alv, has_visco_property
export iso_order2_params_from_blocks, iso_order2_blocks_from_params
export hill_kernel_order2
export voigt_alv_order2, reuss_alv_order2, dilute_alv_order2,
    dilute_dual_alv_order2, mori_tanaka_alv_order2, maxwell_alv_order2
export dilute_concentration_alv_order2, dilute_contribution_alv_order2
export homogenize_alv_order2
export cod_kernel_alv, compliance_contribution_alv, delta_compliance_alv
export stiffness_contribution_alv, stiffness_contribution_alv_at, delta_stiffness_alv

# ── Backwards-compat aliases ─────────────────────────────────────────────────
const HillAlgorithm = AbstractAlgorithm
const CrackAlgorithm = AbstractAlgorithm
export HillAlgorithm, CrackAlgorithm

end # module
