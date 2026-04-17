# =============================================================================
#  bases.jl
#
#  Low-level helpers operating on TensND `AbstractBasis` objects or on
#  the Julia scalar types themselves — shared by every sub-module.
# =============================================================================

"""
    _floatlike(::Type{T})

Promote integer element types to their floating-point counterpart (so that
`Integer`-valued semi-axes become `Float64` or `BigFloat`) while leaving
every `AbstractFloat`, `ForwardDiff.Dual` or symbolic type untouched.
"""
_floatlike(::Type{T}) where {T<:Integer} = float(T)
_floatlike(::Type{T}) where {T<:Number}  = T

"""
    _basis_eltype(::Type{T})

Element type of the default `CanonicalBasis` associated with a given
scalar element type `T`.  Real types always use `Float64` (to match the
historical TensND default); non-real scalars (SymPy, Symbolics, …) keep
their own element type so that the `CanonicalBasis` and the tensor data
share the same element type.
"""
_basis_eltype(::Type{T}) where {T<:Real}  = Float64
_basis_eltype(::Type{T}) where {T<:Number} = T

"""
    _basis_col(basis, m::Int) -> NTuple

Return the `m`-th column of the TensND basis as an `NTuple{d}` — the
coordinates of the `m`-th basis vector in the canonical frame.
"""
function _basis_col(basis::TensND.AbstractBasis, m::Int)
    E = TensND.vecbasis(basis, :cov)
    d = size(E, 1)
    return ntuple(i -> E[i, m], d)
end

"""
    _frame_columns(basis) -> (l̂, m̂, n̂)

Return the three axis vectors of a 3D basis as plain `Vector`s — a
convenience form used by the numerical (residue / DECUHR) algorithms.
The returned element type follows the basis element type.
"""
@inline function _frame_columns(basis::TensND.AbstractBasis)
    l̂ = Vector(TensND.components_canon(TensND.tensbasis(basis, 1)))
    m̂ = Vector(TensND.components_canon(TensND.tensbasis(basis, 2)))
    n̂ = Vector(TensND.components_canon(TensND.tensbasis(basis, 3)))
    return l̂, m̂, n̂
end

"""
    _default_basis(::Type{T}, euler_angles::NTuple{3,<:Real})

Build the default TensND basis associated with a set of ZYZ Euler
angles: `CanonicalBasis{3,_basis_eltype(T)}` when all angles are zero,
`RotatedBasis(euler_angles...)` otherwise.
"""
function _default_basis(::Type{T}, euler_angles::NTuple{3,<:Real}) where {T}
    Tbasis = _basis_eltype(T)
    return all(iszero, euler_angles) ? TensND.CanonicalBasis{3, Tbasis}() :
                                       TensND.RotatedBasis(euler_angles...)
end

# ─── Shape classification helpers ─────────────────────────────────────────────

"""
    _classify_shape_3d(::Type{T}, a, b, c)

Classify a 3D ellipsoid from its (already sorted) semi-axes into one of
the shape traits `Spherical`, `Prolate`, `Oblate`, `Triaxial`.  The
actual shape types are defined in the `Elasticity` sub-module; the
classification logic is delegated to a caller-supplied tuple of type
tags — this keeps `Core` dependency-free from `Elasticity`.

Returns an integer code instead of a type to avoid circular references:
    1 → equivalent of `Spherical`   (a == b == c)
    2 → equivalent of `Prolate`     (a > b == c)
    3 → equivalent of `Oblate`      (a == b > c)
    4 → equivalent of `Triaxial`    (a > b > c)
"""
function _classify_shape_3d(::Type{T}, a, b, c) where {T}
    if T <: Real
        tol  = max(a, b, c) * (1e-10 * one(T))
        AeqB = (a - b) ≤ tol
        BeqC = (b - c) ≤ tol
    else
        AeqB = isequal(a, b)
        BeqC = isequal(b, c)
    end
    if     AeqB &&  BeqC; return 1   # Spherical
    elseif !AeqB &&  BeqC; return 2   # Prolate
    elseif  AeqB && !BeqC; return 3   # Oblate
    else;                   return 4  # Triaxial
    end
end

"""
    _classify_shape_2d(::Type{T}, a, b)

Classify a 2D ellipse from its (already sorted) semi-axes:
    1 → equivalent of `Circular`  (a == b)
    2 → equivalent of `Elliptic`  (a > b)
"""
function _classify_shape_2d(::Type{T}, a, b) where {T}
    is_equal = T <: Real ? (a - b) ≤ max(a, b) * (1e-10 * one(T)) : isequal(a, b)
    return is_equal ? 1 : 2
end
