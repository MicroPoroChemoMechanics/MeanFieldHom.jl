"""
    MeanFieldHom.Cracks

COD tensors, compliance contributions, SIF and DIF for flat cracks
embedded in an elastic matrix of arbitrary anisotropy.  Public entry
points: [`cod_tensor`](@ref), [`compliance_contribution`](@ref),
[`sif`](@ref), [`dif`](@ref).
"""
module Cracks

using LinearAlgebra
using TensND
using Tensors
using QuadGK
using ..Elliptic
using Polynomials
using PolynomialRoots

import ..Core
using ..Core
const MFH_Core = Core

include("geometry.jl")
include("cod_H_bridge.jl")
include("cod_analytical.jl")
include("cod_analytical_thermal.jl")
include("green_residue.jl")
include("green_nestedquadgk.jl")
include("green_decuhr.jl")
include("cod_numerical.jl")
include("compliance.jl")
include("sif.jl")
include("api.jl")

# ── Geometry ─────────────────────────────────────────────────────────────────
export CrackShape, Penny, EllipticShape, Ribbon
export EllipticCrack, RibbonCrack
export PennyCrack
export crack_basis, aspect_ratio, semi_major, semi_minor, crack_normal

# ── COD / compliance ─────────────────────────────────────────────────────────
export cod_tensor, B_tensor
export cod_from_compliance, compliance_from_cod
export compliance_contribution
export delta_compliance, delta_resistivity
# (stiffness_contribution, conductivity_contribution, delta_stiffness,
# delta_conductivity are Core-level generics; methods defined above in
# compliance.jl are visible through Core's exports — nothing to re-export here.)

# ── SIF / DIF ────────────────────────────────────────────────────────────────
export sif, dif

end # module
