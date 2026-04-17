# =============================================================================
#  ellipsoid.jl
#
#  Ellipsoidal-inclusion geometry type subtyping
#  `Core.AbstractEllipsoidalInclusion{dim,T}`. Shape classification
#  happens at construction time and is encoded in the type parameter `S`
#  — no runtime branching in downstream `_kernel` methods.
# =============================================================================

"""
    EllipsoidShape

Abstract supertype for the shape classification of an [`Ellipsoid`](@ref).

Concrete subtypes (3-D): [`Spherical`](@ref), [`Prolate`](@ref),
[`Oblate`](@ref), [`Triaxial`](@ref).

Concrete subtypes (2-D): [`Circular`](@ref), [`Elliptic`](@ref).
"""
abstract type EllipsoidShape end

"3-D sphere: a = b = c."
struct Spherical <: EllipsoidShape end

"3-D prolate spheroid: a > b = c (axis of revolution = e₁)."
struct Prolate <: EllipsoidShape end

"3-D oblate spheroid: a = b > c (axis of revolution = e₃)."
struct Oblate <: EllipsoidShape end

"3-D triaxial ellipsoid: a > b > c."
struct Triaxial <: EllipsoidShape end

"2-D circle: a = b."
struct Circular <: EllipsoidShape end

"2-D ellipse: a > b."
struct Elliptic <: EllipsoidShape end

# Map `_classify_shape_3d` / `_classify_shape_2d` integer codes to shape
# types (`Core` does not depend on `EllipsoidShape`, so we do the mapping
# here, in `Elasticity`).
const _SHAPE_3D = (Spherical, Prolate, Oblate, Triaxial)
const _SHAPE_2D = (Circular, Elliptic)

# =============================================================================
#  Ellipsoid struct
# =============================================================================

"""
    Ellipsoid{dim, S<:EllipsoidShape, T<:Number, B<:AbstractBasis}

Ellipsoidal inclusion with semi-axes `semi_axes` (a₁ ≥ a₂ [≥ a₃] when `T<:Real`,
or in the caller-provided order when `T` is symbolic) and an orientation
`basis` describing the principal frame relative to the global (canonical)
frame.

The shape `S` is determined at construction time:
- 3-D: `Spherical`, `Prolate`, `Oblate`, or `Triaxial`
- 2-D: `Circular` or `Elliptic`

`T` can be any `Number` subtype (`Float64`, `ForwardDiff.Dual`, `SymPy.Sym`,
`Symbolics.Num`, …).
"""
struct Ellipsoid{dim, S<:EllipsoidShape, T<:Number, B<:TensND.AbstractBasis} <:
        MFH_Core.AbstractEllipsoidalInclusion{dim,T}
    semi_axes :: NTuple{dim, T}
    basis     :: B
end

# ── 3-D constructors ──────────────────────────────────────────────────────────

"""
    Ellipsoid(a, b, c; euler_angles=(θ,ϕ,ψ))

3-D ellipsoid with semi-axes `a`, `b`, `c` (sorted descending when `T<:Real`)
oriented by ZYZ Euler angles `(θ,ϕ,ψ)`.  Angles default to `(0,0,0)`.
"""
function Ellipsoid(a::Ta, b::Tb, c::Tc;
                   euler_angles::NTuple{3,<:Real} = (0.0, 0.0, 0.0)) where {Ta,Tb,Tc}
    T = MFH_Core._floatlike(promote_type(Ta, Tb, Tc))
    if T <: Real
        axes_sorted = NTuple{3,T}(sort([T(a), T(b), T(c)]; rev=true))
    else
        axes_sorted = (T(a), T(b), T(c))
    end
    Tbasis = MFH_Core._basis_eltype(T)
    basis  = all(iszero, euler_angles) ? TensND.CanonicalBasis{3, Tbasis}() :
                                         TensND.RotatedBasis(euler_angles...)
    code = MFH_Core._classify_shape_3d(T, axes_sorted...)
    S = _SHAPE_3D[code]
    return Ellipsoid{3, S, T, typeof(basis)}(axes_sorted, basis)
end

"""
    Ellipsoid(a, b, c, R::AbstractMatrix)

3-D ellipsoid whose principal axes are the columns of the rotation matrix `R`.
"""
function Ellipsoid(a::Ta, b::Tb, c::Tc, R::AbstractMatrix) where {Ta,Tb,Tc}
    T = MFH_Core._floatlike(promote_type(Ta, Tb, Tc))
    if T <: Real
        axes_sorted = NTuple{3,T}(sort([T(a), T(b), T(c)]; rev=true))
    else
        axes_sorted = (T(a), T(b), T(c))
    end
    basis = TensND.RotatedBasis(Matrix{Float64}(R))
    code = MFH_Core._classify_shape_3d(T, axes_sorted...)
    S = _SHAPE_3D[code]
    return Ellipsoid{3, S, T, typeof(basis)}(axes_sorted, basis)
end

# ── 2-D constructors ──────────────────────────────────────────────────────────

"""
    Ellipsoid(a, b; angle=0.0)

2-D ellipse with semi-axes `a`, `b` (sorted descending when `T<:Real`) and
orientation angle `θ` (radians) of the major axis w.r.t. the first global axis.
"""
function Ellipsoid(a::Ta, b::Tb; angle::Real = 0.0) where {Ta,Tb}
    T = MFH_Core._floatlike(promote_type(Ta, Tb))
    if T <: Real
        a1, a2 = T(a) >= T(b) ? (T(a), T(b)) : (T(b), T(a))
    else
        a1, a2 = T(a), T(b)
    end
    Tbasis = MFH_Core._basis_eltype(T)
    basis  = iszero(angle) ? TensND.CanonicalBasis{2, Tbasis}() : TensND.RotatedBasis(float(angle))
    code = MFH_Core._classify_shape_2d(T, a1, a2)
    S = _SHAPE_2D[code]
    return Ellipsoid{2, S, T, typeof(basis)}((a1, a2), basis)
end

# ── Sphere / circle ───────────────────────────────────────────────────────────

"""
    Ellipsoid(r; dim=3)

Sphere (3-D) or circle (2-D) of radius `r`.
"""
function Ellipsoid(r::T; dim::Int = 3) where {T<:Number}
    Tf   = MFH_Core._floatlike(T)
    axes = ntuple(_ -> Tf(r), dim)
    basis = TensND.CanonicalBasis{dim, MFH_Core._basis_eltype(Tf)}()
    S = dim == 3 ? Spherical : Circular
    return Ellipsoid{dim, S, Tf, typeof(basis)}(axes, basis)
end

# ── Interface implementations ────────────────────────────────────────────────

MFH_Core.dimension(::Ellipsoid{dim}) where {dim} = dim
MFH_Core.inclusion_basis(ell::Ellipsoid) = ell.basis
MFH_Core.shape_trait(::Ellipsoid{dim,S}) where {dim,S} = S

# ── Shape helpers used by scripts and downstream code ──────────────────────

"Return the spatial dimension of the ellipsoid."
getdim(::Ellipsoid{dim}) where {dim} = dim

"Return the i-th semi-axis (1-indexed)."
semi_axis(ell::Ellipsoid, i::Int) = ell.semi_axes[i]

"Aspect ratio η = a₂/a₁  (≤ 1)."
aspect_ratio_eta(ell::Ellipsoid{3}) = ell.semi_axes[2] / ell.semi_axes[1]

"Aspect ratio ω = a₃/a₁  (≤ η ≤ 1)."
aspect_ratio_omega(ell::Ellipsoid{3}) = ell.semi_axes[3] / ell.semi_axes[1]

"Aspect ratio ρ = a₂/a₁  (≤ 1)."
aspect_ratio_rho(ell::Ellipsoid{2}) = ell.semi_axes[2] / ell.semi_axes[1]

# ── Ellipsoid convenience wrappers around the Core Newton potentials ────────

"""
    newton_potential_3d(ell::Ellipsoid{3})

Ellipsoid-level convenience wrapper that forwards to
[`Core.newton_potential_3d`](@ref).
"""
MFH_Core.newton_potential_3d(ell::Ellipsoid{3}) =
    MFH_Core.newton_potential_3d(ell.semi_axes[1], ell.semi_axes[2], ell.semi_axes[3])

"""
    newton_potential_2d(ell::Ellipsoid{2})

Ellipsoid-level convenience wrapper that forwards to
[`Core.newton_potential_2d`](@ref).
"""
MFH_Core.newton_potential_2d(ell::Ellipsoid{2}) =
    MFH_Core.newton_potential_2d(ell.semi_axes[1], ell.semi_axes[2])
