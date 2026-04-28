"""
    MeanFieldHom.Schemes

Mean-field homogenisation schemes.  Provides the [`RVE`](@ref) container
(matrix + named phases with their volume fractions or crack densities, plus
an optional distribution shape) and the suite of homogenisation
[`HomogenizationScheme`](@ref) types: bounds (`Voigt`, `Reuss`), one-shot
schemes with a matrix (`Dilute`, `DiluteDual`, `MoriTanaka`, `Maxwell`,
`PonteCastanedaWillis`), iterative self-consistent schemes
(`SelfConsistent`, `AsymmetricSelfConsistent`) and the differential
scheme (`Differential`) with user-selectable trajectory.

Public entry point: [`homogenize(rve, scheme; property=:C)`](@ref).  The
scheme can also be passed as a `Symbol` shortcut (`:MT`, `:SC`, …).
"""
module Schemes

using LinearAlgebra
using TensND
using ForwardDiff
using OrdinaryDiffEq

import ..Core
using ..Core
const MFH_Core = Core

# Forward declarations of inclusion types we touch from the other sub-modules
# (loaded earlier than Schemes by `MeanFieldHom.jl`).
import ..Elasticity: Ellipsoid, hill_tensor
import ..Cracks: compliance_contribution, delta_compliance, delta_resistivity

include("rve.jl")
include("symmetrize.jl")
include("scheme_types.jl")
include("homogenize.jl")
include("contribution_helpers.jl")
include("voigt.jl")
include("reuss.jl")
include("dilute.jl")
include("dilute_dual.jl")
include("mori_tanaka.jl")
include("maxwell.jl")
include("pcw.jl")
include("self_consistent.jl")
include("trajectory.jl")
include("differential.jl")
include("parameters.jl")
include("sensitivities.jl")

# ── Exports ────────────────────────────────────────────────────────────────
# Data model
export AbstractAmount, VolumeFraction, CrackDensity
export AbstractDistributionShape, UniformDistribution
export AbstractSymmetrize, NoSymmetrize, IsoSymmetrize, TISymmetrize
export Phase, RVE
export add_matrix!, add_phase!
export matrix_phase, inclusion_phase_names
export phase_property, matrix_property
export volume_fraction, crack_density, matrix_volume_fraction
export phase_symmetrize
export validate_rve

# Schemes
export HomogenizationScheme
export Voigt, Reuss, Dilute, DiluteDual, MoriTanaka, Maxwell, PonteCastanedaWillis
export SelfConsistent, AsymmetricSelfConsistent
export AndersonDefault, NewtonDefault
export DifferentialTrajectory, Proportional, Sequential, CustomPath, Path, DifferentialScheme

# Entry point
export homogenize

# Sensitivities — lentilles paramétriques + wrappers ForwardDiff (extension)
export AbstractParameter, AmountParameter, PropertyParameter,
    GeometryParameter, DistributionShapeParameter
export amount, property, geometry, shape_param
export get_param, set_param
export derivative, gradient, jacobian, sensitivity

end # module
