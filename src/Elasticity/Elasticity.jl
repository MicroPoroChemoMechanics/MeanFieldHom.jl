"""
    MeanFieldHom.Elasticity

Hill polarisation tensors for ellipsoidal inclusions (2D / 3D,
isotropic / anisotropic matrix).  Public entry point:
[`hill_tensor`](@ref).
"""
module Elasticity

using LinearAlgebra
using TensND
using DECUHR
import Integrals
using QuadGK
using ..Elliptic
using Polynomials
using PolynomialRoots

import ..Core
using ..Core
const MFH_Core = Core

# Re-export the abstract inclusion supertype so that users can write
# `AbstractEllipsoidalInclusion` directly through `MeanFieldHom`.
import ..Core: AbstractEllipsoidalInclusion

include("ellipsoid.jl")
include("cylinder.jl")
include("auxiliary_tensors.jl")
include("hill_3d_iso.jl")
include("hill_3d_cylinder_iso.jl")
include("hill_3d_aniso_residue.jl")
include("hill_3d_aniso_nestedquadgk.jl")
include("hill_3d_aniso_decuhr.jl")
include("hill_3d_cylinder_aniso.jl")
include("hill_2d_iso.jl")
include("hill_2d_aniso.jl")
include("api.jl")

export Ellipsoid
export EllipsoidShape, Spherical, Prolate, Oblate, Triaxial, Circular, Elliptic
export Cylinder, CylindricalShape, CircularCylindrical, EllipticCylindrical
export tens_IA, tens_UA, tens_VA
export hill_tensor

end # module
