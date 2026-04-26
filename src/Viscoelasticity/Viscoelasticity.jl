"""
    MeanFieldHom.Viscoelasticity

Ageing linear viscoelastic (ALV) homogenisation.  Provides:

  * [`ViscoLaw`](@ref) — relaxation `R(t,t')` or creep `J(t,t')` kernel,
    scalar- or 4-tensor-valued, with built-in Maxwell / Kelvin
    constructors.
  * [`trapezoidal_matrix`](@ref) — discretisation of the Stieltjes
    integral on a time grid into a lower-block-triangular matrix
    (`n×n` for scalar, `6n×6n` for 4-tensor in Mandel form).
  * [`volterra_inverse`](@ref) — block forward-substitution that takes
    a discrete relaxation kernel to the corresponding creep kernel
    (and vice versa).
  * [`visco_param`](@ref) / [`visco_assemble`](@ref) — conversions
    between symmetry-structured per-component scalar matrices and the
    full `6n×6n` block matrix.
  * `hill_kernel` — discrete ALV Hill polarisation tensor for an
    ellipsoidal inclusion, isotropic-matrix branch using the
    time-space decoupling formula
    [@barthelemyIJSS2016, App. *ALV Hill kernel*].
  * Time-domain viscoelastic homogenisation schemes (Voigt, Reuss,
    Dilute, DiluteDual, Mori-Tanaka, Maxwell, Self-Consistent),
    plugged into the existing [`MeanFieldHom.homogenize`](@ref)
    dispatcher whenever a phase carries a `ViscoLaw` property.

All ALV operators are stored as dense `Matrix{T}` of size `(B·n)×(B·n)`
(`B = 6` for 4-tensor, `B = 1` for scalar) with explicit zeros above
the block diagonal — this is the convention of
[@sanahuja2013] and the C++ ECHOES reference.
"""
module Viscoelasticity

using LinearAlgebra
using TensND

import ..Core
using ..Core
const MFH_Core = Core

import ..Elasticity
import ..Elasticity: tens_UA, tens_VA, Ellipsoid, Spheroid
import ..LayeredSpheres
using ..LayeredSpheres: LayeredSphere, layer_radius, layer_modulus,
                         layer_interface, AbstractInterface, PerfectInterface,
                         SpringInterface, MembraneInterface,
                         layer_count, layer_volume_fraction, outer_radius
import ..Schemes
using ..Schemes: RVE, HomogenizationScheme, Voigt, Reuss, Dilute, DiluteDual,
                  MoriTanaka, Maxwell, SelfConsistent, matrix_phase,
                  inclusion_phase_names, matrix_property, phase_property,
                  volume_fraction, matrix_volume_fraction

include("visco_law.jl")
include("trapezoidal.jl")
include("volterra_inverse.jl")
include("conversions.jl")
include("hill_alv.jl")
include("schemes_alv.jl")
include("iso_schemes_alv.jl")
include("schemes_alv_sc.jl")
include("layered_alv.jl")
include("homogenize_alv.jl")

# ── Exports ─────────────────────────────────────────────────────────────────
export AbstractViscoLaw, ViscoLaw, VALID_VISCO_MODES
export visco_mode, visco_eval
export maxwell_relaxation, kelvin_creep, maxwell_iso, kelvin_iso, heaviside_law
export trapezoidal_matrix
export volterra_inverse, volterra_product, volterra_divide, volterra_left_divide
export iso_params_from_blocks, iso_blocks_from_params
export hill_kernel
export dilute_concentration_alv, dilute_contribution_alv
export voigt_alv, reuss_alv, dilute_alv, dilute_dual_alv
export mori_tanaka_alv, maxwell_alv
export voigt_alv_iso, reuss_alv_iso, dilute_alv_iso, dilute_dual_alv_iso
export mori_tanaka_alv_iso, maxwell_alv_iso
export dilute_concentration_alv_iso, dilute_contribution_alv_iso
export self_consistent_alv
export bulk_localization_alv, bulk_state_seq_alv, shear_localization_alv
export strain_strain_loc_alv, stiffness_contribution_alv
export homogenize_alv, has_visco_property

end # module
