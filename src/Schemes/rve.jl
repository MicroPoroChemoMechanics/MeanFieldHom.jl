# =============================================================================
#  rve.jl — Representative Volume Element (RVE) data model.
#
#  An `RVE` aggregates a *matrix* phase with one or more *inclusion* phases,
#  storing for each phase its geometry (an `AbstractInclusion`), its material
#  properties (a `Symbol => AbstractTens` map: `:C` for stiffness, `:K` for
#  conductivity, …) and an *amount* — a `VolumeFraction` for ellipsoidal
#  inclusions, a `CrackDensity` for cracks.  Volume fractions are stored
#  **at the RVE level**, not on the inclusions: a single inclusion is still
#  usable for localization-tensor calculations (`hill_tensor`,
#  `strain_strain_loc`, …) without any fraction-related machinery.
#
#  The matrix amount is implicit (`f_matrix = 1 - Σ f_inc`) and crack
#  densities are excluded from that sum (their volume contribution → 0
#  in the penny limit).
#
#  PCW / Maxwell additionally need a *distribution shape* describing the
#  outer envelope of the phase distribution; this is stored in the
#  `distribution_shape` field through an `AbstractDistributionShape`
#  hierarchy that allows future extension to pairwise distributions
#  ([Willis 1982](@cite willis1982)) without breaking the public API.
# =============================================================================

# =============================================================================
#  Amounts: volume fraction (ellipsoidal phases) and crack density (cracks)
# =============================================================================

"""
    AbstractAmount{T<:Number}

Supertype for the *quantity* attached to a phase in a `RVE`. Two concrete
subtypes:

- [`VolumeFraction`](@ref) — for solid (ellipsoidal) inclusions and the
  matrix; obeys the unit-sum constraint
  `f_matrix = 1 - Σ_other VolumeFraction.value`.
- [`CrackDensity`](@ref) — for flat cracks (Budiansky-O'Connell density);
  does **not** participate in the unit-sum constraint, since the volume
  contribution of a flat crack vanishes in the penny limit while the
  density remains finite.

The type parameter `T` is the element type of the stored value
(`Float64`, `ForwardDiff.Dual`, `Complex{Float64}`, …) and propagates
through every scheme that consumes the amount.
"""
abstract type AbstractAmount{T <: Number} end

"""
    VolumeFraction(f) <: AbstractAmount

Volume fraction of a solid inclusion (or of the matrix).
"""
struct VolumeFraction{T <: Number} <: AbstractAmount{T}
    value::T
end

"""
    CrackDensity(ε) <: AbstractAmount

Budiansky-O'Connell crack density of a population of flat cracks.
"""
struct CrackDensity{T <: Number} <: AbstractAmount{T}
    value::T
end

"""
    amount_value(a::AbstractAmount) -> Number

Return the scalar value carried by an `AbstractAmount`.
"""
amount_value(a::AbstractAmount) = a.value

Base.eltype(::Type{<:AbstractAmount{T}}) where {T} = T
Base.eltype(a::AbstractAmount) = eltype(typeof(a))

"""
    _sums_to_unit(a::AbstractAmount) -> Bool

Whether the amount counts towards the matrix-fraction complement
`f_matrix = 1 - Σ_phase _sums_to_unit·value`. `true` for
[`VolumeFraction`](@ref), `false` for [`CrackDensity`](@ref).
"""
_sums_to_unit(::VolumeFraction) = true
_sums_to_unit(::CrackDensity) = false

# =============================================================================
#  Symmetrize: orientation-distribution projection of a phase's contribution
# =============================================================================

"""
    AbstractSymmetrize

Specifies how a phase's *localization tensor* (and the derived stiffness /
compliance / conductivity / resistivity contributions) is averaged over an
orientation distribution before being used in the homogenization formula.

Three concrete subtypes are shipped :

- [`NoSymmetrize`](@ref) (default) — keep the contribution as computed for
  the single-orientation inclusion stored in the phase.
- [`IsoSymmetrize`](@ref) — average over **all** rotations (uniform spatial
  distribution of orientations) ; produces an isotropic projection.
- [`TISymmetrize`](@ref) — average over rotations around a specified axis
  (uniaxial uniform distribution) ; produces a transversely-isotropic
  projection.

This mirrors C++ ECHOES's `symmetrize=[ISO]` / `symmetrize=[TI]` keyword on
`ellipsoid()`, but moved to the *RVE* side (just like volume fractions) :
the same inclusion type can be re-used in different RVEs with different
distribution assumptions, and a single inclusion remains usable for
localisation-tensor calculations without any RVE.
"""
abstract type AbstractSymmetrize end

"""
    NoSymmetrize() <: AbstractSymmetrize

Default. The localization tensor is used as computed for the single
orientation defined by the inclusion's basis.
"""
struct NoSymmetrize <: AbstractSymmetrize end

"""
    IsoSymmetrize() <: AbstractSymmetrize

The localization tensor is averaged over a *uniform spatial distribution*
of orientations, equivalent to projecting onto the isotropic basis
`(J, K_proj)` for 4th-order tensors and onto the spherical part for
2nd-order tensors. Produces an isotropic phase contribution regardless of
the inclusion's actual shape.
"""
struct IsoSymmetrize <: AbstractSymmetrize end

"""
    TISymmetrize(axis = (0, 0, 1)) <: AbstractSymmetrize

The localization tensor is averaged over rotations about `axis` (uniaxial
uniform distribution). Produces a transversely-isotropic phase
contribution with that symmetry axis.
"""
struct TISymmetrize{T <: Number} <: AbstractSymmetrize
    axis::NTuple{3, T}
end
TISymmetrize() = TISymmetrize((0.0, 0.0, 1.0))
TISymmetrize(axis::AbstractVector) = TISymmetrize(NTuple{3}(Tuple(axis)))

# Coercer for kwargs : accept a Symbol shortcut, an `AbstractSymmetrize`, or
# nothing (no projection).
_to_symmetrize(::Nothing) = NoSymmetrize()
_to_symmetrize(s::AbstractSymmetrize) = s
function _to_symmetrize(s::Symbol)
    s === :none && return NoSymmetrize()
    s === :iso  && return IsoSymmetrize()
    s === :ISO  && return IsoSymmetrize()
    s === :ti   && return TISymmetrize()
    s === :TI   && return TISymmetrize()
    throw(ArgumentError("unknown symmetrize Symbol :$(s); expected :none, :iso or :ti (or pass an AbstractSymmetrize instance for non-default axis)"))
end

# =============================================================================
#  Distribution shape: PCW / Maxwell outer-envelope descriptor
# =============================================================================

"""
    AbstractDistributionShape

Supertype for the *outer envelope* of the phase distribution used by the
[`Maxwell`](@ref) and [`PonteCastanedaWillis`](@ref) schemes.

Currently a single concrete subtype is shipped:

- [`UniformDistribution`](@ref) — a single shape applied to every
  inclusion phase (Maxwell 1873 ; Ponte-Castañeda & Willis 1995).

Future extension (placeholder, *not* implemented in this PR): a
`PairwiseDistribution` carrying a per-pair `(i, j) ↦ shape` mapping
([Willis 1982](@cite willis1982)).  Adding it will only require a new
concrete subtype + matching `_evaluate(rve, ::Maxwell|::PonteCastanedaWillis, …)`
methods — no public-API change.
"""
abstract type AbstractDistributionShape end

"""
    UniformDistribution(shape::AbstractInclusion) <: AbstractDistributionShape

Single distribution shape applied to every inclusion phase. The default
constructor `UniformDistribution()` returns a unit sphere (isotropic
distribution, recovers Mori-Tanaka in the limit `P_d = P_inc`).
"""
struct UniformDistribution{S <: AbstractInclusion} <: AbstractDistributionShape
    shape::S
end
UniformDistribution() = UniformDistribution(Ellipsoid(1.0))

"""
    distribution_shape_of(d::UniformDistribution) -> AbstractInclusion

Return the inclusion describing the (single) distribution envelope.
"""
distribution_shape_of(d::UniformDistribution) = d.shape

"""
    _to_distribution_shape(x) -> AbstractDistributionShape

Coerce `x` to a concrete `AbstractDistributionShape`. Accepts:

- `nothing` → `UniformDistribution(Ellipsoid(1.0))` (default sphere),
- an `AbstractInclusion` → wrapped as `UniformDistribution(x)`,
- an `AbstractDistributionShape` → passed through.
"""
_to_distribution_shape(::Nothing) = UniformDistribution()
_to_distribution_shape(s::AbstractInclusion) = UniformDistribution(s)
_to_distribution_shape(s::AbstractDistributionShape) = s

# =============================================================================
#  Phase: geometry + material properties
# =============================================================================

"""
    Phase(geometry::AbstractInclusion, properties::Dict{Symbol,<:AbstractTens})

A single phase of a [`RVE`](@ref): one inclusion *geometry* (ellipsoid,
crack, …) together with one or several material *property tensors*
indexed by symbol (`:C` for stiffness, `:K` for conductivity, …).

The geometry is field-typed `AbstractInclusion` (rather than parametric
`Phase{I}`) so that a heterogeneous RVE mixing ellipsoids and cracks can
be stored in a single `Dict{Symbol,Phase}` without losing information at
construction time. Specialization happens at the dispatch site
(`hill_tensor(phase.geometry, …)`, `cod_tensor(phase.geometry, …)`).
"""
mutable struct Phase
    geometry::AbstractInclusion
    properties::Dict{Symbol, Any}
end

Phase(geometry::AbstractInclusion, properties::AbstractDict) =
    Phase(geometry, Dict{Symbol, Any}(properties...))

# =============================================================================
#  RVE: ordered collection of phases + matrix tag + distribution shape
# =============================================================================

"""
    RVE{T<:Number, S<:Union{Nothing,AbstractDistributionShape}}

Multi-phase representative volume element. Fields:

- `matrix_name::Symbol` — name of the phase that plays the role of the
  matrix (its amount is implicit, computed as `1 - Σ_inc f_inc`).
- `phase_names::Vector{Symbol}` — phases in insertion order (the matrix
  is the first entry).
- `phases::Dict{Symbol,Phase}` — geometry + properties of each phase.
- `amounts::Dict{Symbol,AbstractAmount{T}}` — volume fraction or crack
  density of each non-matrix phase. The matrix entry, if present, is
  ignored when computing `matrix_volume_fraction`.
- `distribution_shape::S` — outer envelope used by Maxwell / PCW;
  defaults to a unit sphere wrapped in [`UniformDistribution`](@ref).

`T` is the element type of every amount in the RVE — it drives the
propagation of `ForwardDiff.Dual` / `Complex{Float64}` through fractions
independently of the moduli, which can carry their own element type on
each phase.

Construction is two-step:

```julia
rve = RVE(:M; T = Float64, distribution_shape = nothing)
add_matrix!(rve, ellipsoid_matrix, Dict(:C => C0))
add_phase!(rve, :I1, ellipsoid_inc, Dict(:C => C1); fraction = 0.2)
add_phase!(rve, :CRACK, penny_crack, Dict(:C => C0); density = 0.05)
```

See also [`add_matrix!`](@ref), [`add_phase!`](@ref),
[`matrix_volume_fraction`](@ref), [`validate_rve`](@ref).
"""
mutable struct RVE{T <: Number, S <: Union{Nothing, AbstractDistributionShape}}
    matrix_name::Symbol
    phase_names::Vector{Symbol}
    phases::Dict{Symbol, Phase}
    amounts::Dict{Symbol, AbstractAmount{T}}
    symmetrize::Dict{Symbol, AbstractSymmetrize}
    distribution_shape::S
end

"""
    RVE(matrix_name::Symbol; T = Float64, distribution_shape = nothing)

Construct an empty RVE. The matrix phase is referenced by `matrix_name`
but **not** added — call [`add_matrix!`](@ref) next. Element type of
the amounts is fixed by the `T` keyword (default `Float64`); use
`T = ForwardDiff.Dual{...}` or `T = Complex{Float64}` for AD or
frequency-domain workflows.
"""
function RVE(
        matrix_name::Symbol;
        T::Type{<:Number} = Float64,
        distribution_shape = nothing
    )
    ds = _to_distribution_shape(distribution_shape)
    return RVE{T, typeof(ds)}(
        matrix_name,
        Symbol[],
        Dict{Symbol, Phase}(),
        Dict{Symbol, AbstractAmount{T}}(),
        Dict{Symbol, AbstractSymmetrize}(),
        ds,
    )
end

# =============================================================================
#  Mutators
# =============================================================================

"""
    add_matrix!(rve, geometry, properties::AbstractDict; symmetrize = nothing)

Register the matrix phase. Must be called before any [`add_phase!`](@ref).
The matrix has no explicit amount (its volume fraction is implicit).

Pass a `symmetrize = :iso | :ti | TISymmetrize(axis) | NoSymmetrize()` kwarg
to declare an orientation-distribution projection of the matrix's
localization tensor (see [`AbstractSymmetrize`](@ref)).
"""
function add_matrix!(
        rve::RVE, geometry::AbstractInclusion, properties::AbstractDict;
        symmetrize = nothing
    )
    name = rve.matrix_name
    haskey(rve.phases, name) &&
        throw(ArgumentError("matrix phase :$(name) already registered"))
    rve.phases[name] = Phase(geometry, properties)
    pushfirst!(rve.phase_names, name)
    sym = _to_symmetrize(symmetrize)
    if !(sym isa NoSymmetrize)
        rve.symmetrize[name] = sym
    end
    return rve
end

"""
    add_phase!(rve, name::Symbol, geometry, properties::AbstractDict;
               fraction = nothing, density = nothing, symmetrize = nothing)

Register an inclusion phase with the given `geometry` and material
`properties`. Exactly one of `fraction` (for ellipsoidal inclusions and
solid inhomogeneities) or `density` (for cracks) must be supplied;
`fraction` produces a [`VolumeFraction`](@ref), `density` a
[`CrackDensity`](@ref).

Both `fraction` and `density` are converted to the RVE's amount eltype
`T` at insertion.

The optional `symmetrize` kwarg declares an orientation-distribution
projection of this phase's localization tensor : `:iso` (uniform spatial
distribution → isotropic projection), `:ti` (uniaxial uniform around
z-axis), `TISymmetrize(axis)` (around an arbitrary axis), or pass an
explicit [`AbstractSymmetrize`](@ref) instance. The default
[`NoSymmetrize`](@ref) keeps the inclusion's actual single-orientation
tensor.
"""
function add_phase!(
        rve::RVE{T}, name::Symbol, geometry::AbstractInclusion,
        properties::AbstractDict;
        fraction = nothing, density = nothing,
        symmetrize = nothing
    ) where {T}
    name === rve.matrix_name &&
        throw(ArgumentError("name :$(name) is reserved for the matrix phase"))
    haskey(rve.phases, name) &&
        throw(ArgumentError("phase :$(name) is already registered"))
    (fraction === nothing) == (density === nothing) &&
        throw(ArgumentError("specify exactly one of `fraction=…` or `density=…`"))

    rve.phases[name] = Phase(geometry, properties)
    push!(rve.phase_names, name)
    rve.amounts[name] = if fraction !== nothing
        VolumeFraction{T}(convert(T, fraction))
    else
        CrackDensity{T}(convert(T, density))
    end
    sym = _to_symmetrize(symmetrize)
    if !(sym isa NoSymmetrize)
        rve.symmetrize[name] = sym
    end
    return rve
end

# =============================================================================
#  Accessors
# =============================================================================

"""
    matrix_phase(rve::RVE) -> Phase

Return the matrix `Phase`. Errors if the matrix has not been registered
(call [`add_matrix!`](@ref) first).
"""
function matrix_phase(rve::RVE)
    haskey(rve.phases, rve.matrix_name) ||
        throw(ArgumentError("matrix phase :$(rve.matrix_name) is not yet registered — call add_matrix! first"))
    return rve.phases[rve.matrix_name]
end

"""
    inclusion_phase_names(rve::RVE) -> Vector{Symbol}

Names of the non-matrix phases in insertion order.
"""
inclusion_phase_names(rve::RVE) =
    Symbol[n for n in rve.phase_names if n != rve.matrix_name]

"""
    phase_property(rve, name::Symbol, key::Symbol) -> AbstractTens

Return the property tensor `key` (e.g. `:C`, `:K`) of the phase named
`name`.
"""
function phase_property(rve::RVE, name::Symbol, key::Symbol)
    haskey(rve.phases, name) ||
        throw(ArgumentError("no phase named :$(name) in RVE"))
    p = rve.phases[name]
    haskey(p.properties, key) ||
        throw(ArgumentError("phase :$(name) does not carry property :$(key)"))
    return p.properties[key]
end

"""
    matrix_property(rve, key::Symbol) -> AbstractTens

Shortcut for `phase_property(rve, rve.matrix_name, key)`.
"""
matrix_property(rve::RVE, key::Symbol) = phase_property(rve, rve.matrix_name, key)

"""
    volume_fraction(rve, name::Symbol) -> Number

Volume fraction of phase `name`. Returns `zero(T)` if the phase carries a
[`CrackDensity`](@ref) instead of a [`VolumeFraction`](@ref).
"""
function volume_fraction(rve::RVE{T}, name::Symbol) where {T}
    name === rve.matrix_name && return matrix_volume_fraction(rve)
    a = rve.amounts[name]
    return a isa VolumeFraction ? amount_value(a) : zero(T)
end

"""
    crack_density(rve, name::Symbol) -> Number

Crack density of phase `name`. Returns `zero(T)` if the phase carries a
[`VolumeFraction`](@ref) instead of a [`CrackDensity`](@ref).
"""
function crack_density(rve::RVE{T}, name::Symbol) where {T}
    haskey(rve.amounts, name) || return zero(T)
    a = rve.amounts[name]
    return a isa CrackDensity ? amount_value(a) : zero(T)
end

"""
    phase_symmetrize(rve, name::Symbol) -> AbstractSymmetrize

Return the orientation-distribution projection declared for phase `name`.
Defaults to [`NoSymmetrize`](@ref) if none was set.
"""
phase_symmetrize(rve::RVE, name::Symbol) =
    get(rve.symmetrize, name, NoSymmetrize())

"""
    matrix_volume_fraction(rve::RVE) -> Number

Implicit matrix volume fraction `1 - Σ_inc f_inc` (only
[`VolumeFraction`](@ref) entries contribute; [`CrackDensity`](@ref)
entries are ignored).
"""
function matrix_volume_fraction(rve::RVE{T}) where {T}
    f_inc = zero(T)
    for (_, a) in rve.amounts
        if _sums_to_unit(a)
            f_inc += amount_value(a)
        end
    end
    return one(T) - f_inc
end

# =============================================================================
#  Validation
# =============================================================================

"""
    validate_rve(rve::RVE)

Sanity-check the RVE: matrix registered, all amounts non-negative, sum
of `VolumeFraction` entries ≤ 1.  Throws `ArgumentError` on
hard failures; emits `@warn` if `f_inc > 1` (non-physical RVE — useful
for symbolic / Dual exploration but flagged).
"""
function validate_rve(rve::RVE)
    haskey(rve.phases, rve.matrix_name) ||
        throw(ArgumentError("RVE has no matrix phase :$(rve.matrix_name); call add_matrix! first"))
    for (name, a) in rve.amounts
        v = amount_value(a)
        if v isa Real && v < 0
            throw(ArgumentError("phase :$(name) has negative amount $(v)"))
        end
    end
    fm = matrix_volume_fraction(rve)
    if fm isa Real && fm < 0
        @warn "RVE has matrix volume fraction $(fm) < 0 — total inclusion volume fraction exceeds 1"
    end
    return rve
end

# =============================================================================
#  Pretty printing
# =============================================================================

function Base.show(io::IO, ::MIME"text/plain", rve::RVE{T, S}) where {T, S}
    println(io, "RVE{$T} with ", length(rve.phase_names), " phase(s)")
    println(io, "  matrix : :$(rve.matrix_name)")
    for name in rve.phase_names
        name === rve.matrix_name && continue
        a = rve.amounts[name]
        kind = a isa VolumeFraction ? "f" : "ε"
        println(io, "  inclusion : :$(name)   $kind = $(amount_value(a))")
    end
    print(io, "  distribution_shape : ", rve.distribution_shape)
    return
end

Base.show(io::IO, rve::RVE{T}) where {T} =
    print(io, "RVE{$T}(:$(rve.matrix_name), $(length(rve.phase_names)) phases)")
