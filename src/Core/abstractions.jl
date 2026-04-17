# =============================================================================
#  abstractions.jl
#
#  Unified inclusion-type hierarchy used throughout `MeanFieldHom`.
#
#  Every inclusion geometry (ellipsoid, crack, multi-layer, user-defined, …)
#  subtypes `AbstractInclusion{T}`, where `T` is the element type of the
#  geometric scalars (semi-axes, half-widths, radii, …).  Sub-hierarchies
#  organise the dispatch tables at the next level:
#
#     AbstractInclusion{T}
#       ├── AbstractEllipsoidalInclusion{dim,T}   — ellipsoids (2D / 3D)
#       ├── AbstractCrack{T}                      — flat cracks (3D)
#       └── AbstractLayeredInclusion{dim,T}       — multi-layer (scaffold)
#
#  The interface below (`dimension`, `element_type`, `inclusion_basis`,
#  `shape_trait`) is the minimal contract every inclusion is expected to
#  implement.  It is deliberately declared here as *stub* `function`
#  definitions (no methods) so that sub-modules can add their own methods
#  without ambiguities — the sub-module always does
#
#      import ..Core: dimension, inclusion_basis, shape_trait
#      Core.dimension(::MyInclusion) = …
#
#  (`element_type` has a single generic fallback — see below.)
# =============================================================================

"""
    AbstractInclusion{T<:Number}

Root abstract supertype for every inclusion geometry recognised by
`MeanFieldHom`.  The type parameter `T` is the element type of the
geometric scalars (semi-axes, half-widths, …) and propagates through
every tensor produced by the package, supporting `Float64`,
`ForwardDiff.Dual`, `SymPy.Sym`, `Symbolics.Num`, …
"""
abstract type AbstractInclusion{T<:Number} end

"""
    AbstractEllipsoidalInclusion{dim,T} <: AbstractInclusion{T}

Supertype for ellipsoidal inclusions — solid ellipsoids (and their
degenerate limits: spheres, cylinders, discs …).  The first type
parameter `dim` encodes the spatial dimension (2 or 3).
"""
abstract type AbstractEllipsoidalInclusion{dim,T} <: AbstractInclusion{T} end

"""
    AbstractCrack{T} <: AbstractInclusion{T}

Supertype for flat-crack geometries (elliptic / ribbon / penny).  Cracks
always live in 3D physical space — the spatial dimension of the geometry
(= 2 in the crack plane) is not exposed in the type, as downstream
algorithms uniformly operate on the 3D stiffness tensor.
"""
abstract type AbstractCrack{T} <: AbstractInclusion{T} end

"""
    AbstractLayeredInclusion{dim,T} <: AbstractInclusion{T}

Scaffold supertype for multi-layer inclusions (spheres with concentric
shells, cylinders with coatings, …).  No concrete subtype is shipped
yet; see `docs/src/developer/roadmap.md`.
"""
abstract type AbstractLayeredInclusion{dim,T} <: AbstractInclusion{T} end

# ─── Minimal interface ───────────────────────────────────────────────────────

"""
    dimension(incl::AbstractInclusion) -> Int

Spatial dimension of the inclusion's ambient space (2 or 3 for the
concrete inclusions shipped with the package).
"""
function dimension end

"""
    element_type(incl::AbstractInclusion{T}) -> Type{T}

Element type of the geometric scalars stored in the inclusion
(`Float64`, `ForwardDiff.Dual`, `SymPy.Sym`, …).
"""
element_type(::AbstractInclusion{T}) where {T} = T

"""
    inclusion_basis(incl::AbstractInclusion) -> TensND.AbstractBasis

Local principal basis of the inclusion (principal frame for an
ellipsoid, ``(\\hat l, \\hat m, \\hat n)`` for a crack, …).  Used by
downstream algorithms to rotate the matrix stiffness / conductivity
into the inclusion frame.
"""
function inclusion_basis end

"""
    shape_trait(incl::AbstractInclusion) -> Type

Concrete shape classification of the inclusion, used as a type
parameter for Holy-style dispatch in the downstream kernels.  Typical
values: `Spherical`, `Prolate`, `Oblate`, `Triaxial`, `Circular`,
`Elliptic` (ellipsoids), `Penny`, `EllipticShape`, `Ribbon` (cracks).
"""
function shape_trait end
