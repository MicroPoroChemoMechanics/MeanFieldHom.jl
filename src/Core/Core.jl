"""
    MeanFieldHom.Core

Shared abstractions and numerical kernels used throughout `MeanFieldHom`.

Contents
--------
- `abstractions.jl`       : inclusion hierarchy (`AbstractInclusion` …)
- `traits.jl`             : algorithm and material-symmetry traits
- `bases.jl`              : helpers around `TensND` bases
- `tensor_helpers.jl`     : low-level utilities (`_δ`, `_C_array`, Voigt)
- `rotational_average.jl` : exact SO(3) / azimuthal averages (ISO, TI) of
                             minor-symmetric tensors, incl. Mandel-block forms
- `moduli.jl`             : modulus extractors for the common symmetry classes
- `newton_potential.jl`   : Newton potentials (2D / 3D)
- `green_kernel.jl`       : acoustic tensor and its adjugate / determinant
- `green_residue.jl`      : Masson / Cauchy residue summation
- `green_helpers.jl`      : quadrature-agnostic Green-function helpers
- `quadrature.jl`         : uniform wrappers around `quadgk` and `hcubature`
- `dispatch.jl`           : central `_resolve_algo` mechanism
"""
module Core

using LinearAlgebra
using TensND
using QuadGK
using ..Elliptic
using Polynomials
using PolynomialRoots

include("abstractions.jl")
include("traits.jl")
include("bases.jl")
include("tensor_helpers.jl")
include("rotational_average.jl")
include("moduli.jl")
include("newton_potential.jl")
include("green_kernel.jl")
include("green_residue.jl")
include("green_helpers.jl")
include("quadrature.jl")
include("dispatch.jl")

# Abstractions
export AbstractInclusion, AbstractEllipsoidalInclusion,
    AbstractCrack, AbstractLayeredInclusion
export dimension, element_type, inclusion_basis, shape_trait, shape_tensor
export eshelby_tensor

# Traits — algorithms
export AbstractAlgorithm, Analytical, Residue, DECUHR, NestedQuadGK,
    CylinderQuadrature, Auto

# Traits — material symmetry
export MaterialSymmetry, IsotropicSym, TransverselyIsotropicSym,
    OrthotropicSym, GeneralAnisotropicSym, material_symmetry

# Modulus extractors (public — consumed by sub-modules and users)
export extract_iso_moduli, extract_ti_moduli, extract_iso_conductivity

# Exact rotation-group averages (public — used by Schemes, ALV and users)
export isotropify, transverse_isotropify
export ti_average_mandel66, iso_average_mandel66
export mandel66_minor, array_from_mandel66

# Newton potentials (public — used downstream and in tests)
export newton_potential_3d, newton_potential_2d, newton_potential_3d_cylinder

# Localization & contribution (generics; methods added at top level and in Cracks)
export strain_strain_loc, stress_strain_loc, strain_stress_loc, stress_stress_loc
export gradient_gradient_loc, flux_gradient_loc, gradient_flux_loc, flux_flux_loc
export stiffness_contribution, conductivity_contribution, resistivity_contribution
export delta_stiffness, delta_conductivity

end # module
