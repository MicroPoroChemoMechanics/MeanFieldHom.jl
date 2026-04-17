# =============================================================================
#  geometry.jl
#
#  Crack geometry types subtyping `Core.AbstractCrack{T}`.
# =============================================================================

"""
    CrackShape

Abstract supertype used as the second type parameter of
[`EllipticCrack`](@ref).  Concrete subtypes: [`Penny`](@ref) (a == b) and
[`EllipticShape`](@ref) (a > b).
"""
abstract type CrackShape end

"Circular (penny-shaped) flat crack: a == b."
struct Penny <: CrackShape end

"General flat elliptical crack: a > b."
struct EllipticShape <: CrackShape end

"""
    Ribbon

Placeholder type tag for ribbon-like cracks.  Mirrors the role of
[`CrackShape`](@ref) for elliptical cracks.
"""
struct Ribbon end

# =============================================================================
#  EllipticCrack
# =============================================================================

"""
    EllipticCrack{T, S<:CrackShape, B<:AbstractBasis}

Flat elliptical crack with semi-axes `a ≥ b` oriented by a local basis
`ℬ = (l̂, m̂, n̂)`.
"""
struct EllipticCrack{T<:Number, S<:CrackShape, B<:TensND.AbstractBasis} <:
        MFH_Core.AbstractCrack{T}
    a     :: T
    b     :: T
    basis :: B
end

# =============================================================================
#  RibbonCrack
# =============================================================================

"""
    RibbonCrack{T, B<:AbstractBasis}

Ribbon-like (tunnel) crack of half-width `b` along `m̂`, unbounded along
`l̂`, with normal `n̂`.
"""
struct RibbonCrack{T<:Number, B<:TensND.AbstractBasis} <: MFH_Core.AbstractCrack{T}
    b     :: T
    basis :: B
end

# =============================================================================
#  Internal helpers
# =============================================================================

function _classify_crack(::Type{T}, a, b) where {T}
    if T <: Real
        tol = max(a, b) * (1e-10 * one(T))
        return (a - b) ≤ tol ? Penny : EllipticShape
    else
        return isequal(a, b) ? Penny : EllipticShape
    end
end

# =============================================================================
#  Constructors
# =============================================================================

"""
    EllipticCrack(a, b; euler_angles=(0,0,0))

Elliptical flat crack with semi-axes `a` and `b` oriented by ZYZ Euler
angles.  When `T <: Real`, the semi-axes are sorted so that `a ≥ b`.
"""
function EllipticCrack(a::Ta, b::Tb;
                       euler_angles::NTuple{3,<:Real} = (0.0, 0.0, 0.0)) where {Ta,Tb}
    T = MFH_Core._floatlike(promote_type(Ta, Tb))
    a_, b_ = T <: Real && T(b) > T(a) ? (T(b), T(a)) : (T(a), T(b))
    basis  = MFH_Core._default_basis(T, euler_angles)
    S      = _classify_crack(T, a_, b_)
    return EllipticCrack{T, S, typeof(basis)}(a_, b_, basis)
end

"""
    EllipticCrack(a, b, R::AbstractMatrix)

Elliptical flat crack whose local frame columns are the columns of `R`.
"""
function EllipticCrack(a::Ta, b::Tb, R::AbstractMatrix) where {Ta,Tb}
    T = MFH_Core._floatlike(promote_type(Ta, Tb))
    basis = TensND.RotatedBasis(Matrix{Float64}(R))
    S = _classify_crack(T, T(a), T(b))
    return EllipticCrack{T, S, typeof(basis)}(T(a), T(b), basis)
end

"""
    EllipticCrack(a, b, basis::AbstractBasis)

Elliptical flat crack with an already-constructed TensND basis.
"""
function EllipticCrack(a::Ta, b::Tb, basis::TensND.AbstractBasis) where {Ta,Tb}
    T = MFH_Core._floatlike(promote_type(Ta, Tb))
    S = _classify_crack(T, T(a), T(b))
    return EllipticCrack{T, S, typeof(basis)}(T(a), T(b), basis)
end

"""
    PennyCrack(a; euler_angles=(0,0,0))

Convenience constructor for a circular (penny-shaped) flat crack.
"""
PennyCrack(a; euler_angles::NTuple{3,<:Real}=(0.0, 0.0, 0.0)) =
    EllipticCrack(a, a; euler_angles=euler_angles)

"""
    RibbonCrack(b; euler_angles=(0,0,0))

Ribbon-like (tunnel) crack of half-width `b`.
"""
function RibbonCrack(b::Tb;
                     euler_angles::NTuple{3,<:Real} = (0.0, 0.0, 0.0)) where {Tb}
    T = MFH_Core._floatlike(Tb)
    basis = MFH_Core._default_basis(T, euler_angles)
    return RibbonCrack{T, typeof(basis)}(T(b), basis)
end

"""
    RibbonCrack(b, R::AbstractMatrix)

Ribbon-like crack with local frame columns `R = [l̂ | m̂ | n̂]`.
"""
function RibbonCrack(b::Tb, R::AbstractMatrix) where {Tb}
    T = MFH_Core._floatlike(Tb)
    basis = TensND.RotatedBasis(Matrix{Float64}(R))
    return RibbonCrack{T, typeof(basis)}(T(b), basis)
end

"""
    RibbonCrack(b, basis::AbstractBasis)

Ribbon-like crack with an already-constructed TensND basis.
"""
function RibbonCrack(b::Tb, basis::TensND.AbstractBasis) where {Tb}
    T = MFH_Core._floatlike(Tb)
    return RibbonCrack{T, typeof(basis)}(T(b), basis)
end

# =============================================================================
#  Interface implementations (Core abstractions)
# =============================================================================

MFH_Core.dimension(::MFH_Core.AbstractCrack) = 3
MFH_Core.inclusion_basis(c::MFH_Core.AbstractCrack) = c.basis
MFH_Core.shape_trait(::EllipticCrack{T,S}) where {T,S} = S
MFH_Core.shape_trait(::RibbonCrack) = Ribbon

# =============================================================================
#  Accessors
# =============================================================================

"Return the local basis ``(\\hat{\\mathbf l}, \\hat{\\mathbf m}, \\hat{\\mathbf n})`` of the crack."
crack_basis(c::MFH_Core.AbstractCrack) = c.basis

"Return the unit normal ``\\hat{\\mathbf n}`` of the crack plane."
crack_normal(c::MFH_Core.AbstractCrack) = TensND.tensbasis(c.basis, 3)

"""
    aspect_ratio(c)

Aspect ratio ``η = b/a`` for an elliptical crack, or `zero(T)` for a
[`RibbonCrack`](@ref) (limit case).
"""
aspect_ratio(c::EllipticCrack) = c.b / c.a
aspect_ratio(c::RibbonCrack{T}) where {T} = zero(T)

"Semi-major axis of an elliptical crack."
semi_major(c::EllipticCrack) = c.a
semi_major(c::RibbonCrack) = error("`semi_major` is undefined for a ribbon crack (infinite).  Use `semi_minor` for the half-width.")

"Semi-minor axis of an elliptical crack, or half-width of a ribbon crack."
semi_minor(c::EllipticCrack) = c.b
semi_minor(c::RibbonCrack)   = c.b

"""
    crack_chi(c)

Dimensionless coefficient ``χ`` linking the average crack opening to the
maximum opening:  ``χ^{\\mathcal E} = 2/3`` for an [`EllipticCrack`](@ref),
``χ^{\\mathcal R} = π/4`` for a [`RibbonCrack`](@ref).
"""
crack_chi(::EllipticCrack{T}) where {T<:Number} = T(2) / T(3)
crack_chi(::RibbonCrack{T})   where {T<:Number} = T(π) / T(4)

# =============================================================================
#  Pretty printing
# =============================================================================

function Base.show(io::IO, c::EllipticCrack{T,S}) where {T,S}
    shape = S === Penny ? "penny-shaped" : "elliptic"
    print(io, "EllipticCrack{", T, "} (", shape, ", a=", c.a, ", b=", c.b, ")")
end

function Base.show(io::IO, c::RibbonCrack{T}) where {T}
    print(io, "RibbonCrack{", T, "} (half-width b=", c.b, ")")
end
