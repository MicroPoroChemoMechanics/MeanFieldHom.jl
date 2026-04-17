"""
    MeanFieldHom.Core

Shared abstractions and numerical kernels used throughout `MeanFieldHom`.

Contents
--------
- `abstractions.jl`       : inclusion hierarchy (`AbstractInclusion` …)
- `traits.jl`             : algorithm and material-symmetry traits
- `bases.jl`              : helpers around `TensND` bases
- `tensor_helpers.jl`     : low-level utilities (`_δ`, `_C_array`, Voigt)
- `moduli.jl`             : modulus extractors for the common symmetry classes
- `newton_potential.jl`   : Newton potentials (2D / 3D)
- `green_kernel.jl`       : acoustic tensor and its adjugate / determinant
- `green_residue.jl`      : Masson / Cauchy residue summation
- `green_decuhr.jl`       : shared 2D DECUHR / QuadGK integrand
- `quadrature.jl`         : uniform wrappers around `quadgk` and `hcubature`
- `dispatch.jl`           : central `_resolve_algo` mechanism
"""
module Core

using LinearAlgebra
using TensND
using DECUHR
using QuadGK
using GenericElliptic
using Polynomials
using PolynomialRoots

include("abstractions.jl")
include("traits.jl")
include("bases.jl")
include("tensor_helpers.jl")
include("moduli.jl")
include("newton_potential.jl")
include("green_kernel.jl")
include("green_residue.jl")
include("green_decuhr.jl")
include("quadrature.jl")
include("dispatch.jl")

# Abstractions
export AbstractInclusion, AbstractEllipsoidalInclusion,
    AbstractCrack, AbstractLayeredInclusion
export dimension, element_type, inclusion_basis, shape_trait

# Traits — algorithms
export AbstractAlgorithm, Analytical, Residue, DECUHR, Auto

# Traits — material symmetry
export MaterialSymmetry, IsotropicSym, TransverselyIsotropicSym,
    OrthotropicSym, GeneralAnisotropicSym, material_symmetry

# Modulus extractors (public — consumed by sub-modules and users)
export extract_iso_moduli, extract_ti_moduli, extract_iso_conductivity

# Newton potentials (public — used downstream and in tests)
export newton_potential_3d, newton_potential_2d

end # module
