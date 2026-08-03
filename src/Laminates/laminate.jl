# =============================================================================
#  laminate.jl — the `Laminate` cell: a periodic stack of parallel layers.
#
#  Deliberately built to mirror `RVE` (two-step construction, per-member
#  property dicts, per-member element types, immutable rebuild) so that the
#  parameter lenses, the declarative nesting and the sensitivity wrappers
#  serve both cells through the same `AbstractHomogenizationCell` contract.
#
#  What it is NOT: an `AbstractInclusion`. A laminate is a *cell* — a unit of
#  homogenization — not something embeddable in a matrix. Embedding a
#  laminated inclusion would require its Hill tensor, which is a different
#  (and open) problem.
# =============================================================================

"""
    Layer(properties)

One layer of a [`Laminate`](@ref). Carries only material properties: unlike
a [`Phase`](@ref), a layer has no inclusion geometry of its own — its
geometry is the stacking direction, shared by the whole cell.

A property value may itself be a [`Homogenized`](@ref) cell + scheme,
resolved lazily at `homogenize` time (declarative multiscale nesting).
"""
mutable struct Layer
    properties::Dict{Symbol, Any}
end

Layer(properties::AbstractDict) = Layer(Dict{Symbol, Any}(properties...))

"""
    Laminate{T<:Number, B<:TensND.AbstractBasis{3}} <: AbstractHomogenizationCell

Periodic unit cell of parallel layers normal to `n` — the deterministic,
matrix-free counterpart of an [`RVE`](@ref), solved exactly by
[`Laminated`](@ref).

Fields:

- `layer_names::Vector{Symbol}` — layers in **stacking order**;
- `layers::Dict{Symbol,Layer}` — per-layer property dicts;
- `thicknesses::Dict{Symbol,Any}` — `h_i`, each keeping its own element type
  (`Float64`, `Dual`, `Sym`, …), exactly as `RVE.amounts` does;
- `interfaces::Vector{AbstractInterface}` — one per layer; entry `k` is the
  interface **on top of** layer `k`, entry `N` closing the cell back onto
  layer 1 by periodicity;
- `basis::B` — orthonormal `(ℓ, m, n)`, whose **third** axis is the layer
  normal. Storing a basis rather than a bare `n` fixes `(ℓ, m)`
  deterministically, which matters for anisotropic layers;
- `period` — cached `L = Σ h_i`.

Thicknesses rather than fractions are the stored primitive: the physics needs
both `f_i = h_i/L` and `L`, and one dict gives both. `L` is what carries the
**size effect** of imperfect interfaces, which enter with weight `1/L` — an
interface *density*. With perfect interfaces the result depends on the `f_i`
alone and `L` is irrelevant.

`T` is the declared element-type **floor**, exactly as in [`RVE`](@ref): a
`Dual` or symbolic thickness lives happily in a plain `Laminate()`.

Construction is two-step, mirroring `RVE`:

```julia
lam = Laminate(; normal = (0, 0, 1))
add_layer!(lam, :A, Dict(:C => C_A); fraction = 0.4)
add_layer!(lam, :B, Dict(:C => C_B); fraction = 0.6,
           interface = SpringInterface(1.0e-3, 2.0e-3))
homogenize(lam, Laminated(), :C)
```

See also [`add_layer!`](@ref), [`layer_property`](@ref),
[`laminate_period`](@ref), [`validate_laminate`](@ref).
"""
mutable struct Laminate{T <: Number, B <: TensND.AbstractBasis{3}} <:
    MFH_Core.AbstractHomogenizationCell
    layer_names::Vector{Symbol}
    layers::Dict{Symbol, Layer}
    thicknesses::Dict{Symbol, Any}
    interfaces::Vector{AbstractInterface}
    basis::B
    # Cached `L = Σ h_i`, refreshed by `add_layer!` (the only writer of
    # `thicknesses` on a live cell), so that reading it is a field access.
    period::Any
end

"""
    Laminate(; normal = (0, 0, 1), euler_angles = nothing, basis = nothing,
             T = Float64)

Construct an empty laminate. Layers are added next with [`add_layer!`](@ref).

The frame is given in exactly one of three ways (or defaulted):

- `normal = (nx, ny, nz)` — completed into `(ℓ, m, n̂)` by
  `Core._frame_from_normal`; the default `(0, 0, 1)` yields a
  `CanonicalBasis`, and the kernel then skips the frame rotation entirely;
- `euler_angles = (θ, ϕ, ψ)` — ZYZ angles, as everywhere else in the package;
- `basis = …` — an explicit `TensND` basis whose third axis is the normal
  (the only route available for a symbolic frame).

`T` declares the element-type floor of the thicknesses; like an `RVE`'s, it is
a floor and not a constraint — a wider thickness is stored as such.
"""
function Laminate(;
        normal = nothing,
        euler_angles = nothing,
        basis = nothing,
        T::Type{<:Number} = Float64
    )
    n_given = count(!isnothing, (normal, euler_angles, basis))
    n_given ≤ 1 || throw(
        ArgumentError(
            "Laminate: give at most one of `normal`, `euler_angles`, `basis` " *
                "to define the frame; got $(n_given)"
        )
    )
    b = if basis !== nothing
        basis isa TensND.AbstractBasis{3} ||
            throw(ArgumentError("Laminate: `basis` must be a 3-D TensND basis"))
        basis
    elseif euler_angles !== nothing
        MFH_Core._default_basis(T, euler_angles isa Tuple ? euler_angles : (euler_angles,))
    elseif normal !== nothing
        MFH_Core._frame_from_normal(normal)
    else
        TensND.CanonicalBasis{3, MFH_Core._basis_eltype(T)}()
    end
    return Laminate{T, typeof(b)}(
        Symbol[], Dict{Symbol, Layer}(), Dict{Symbol, Any}(),
        AbstractInterface[], b, zero(T),
    )
end

Laminate{T}(; kw...) where {T <: Number} = Laminate(; T = T, kw...)

# =============================================================================
#  Mutators
# =============================================================================

"""
    add_layer!(lam, name::Symbol, properties::AbstractDict;
               thickness = nothing, fraction = nothing,
               interface = PerfectInterface())

Append a layer to the top of the stack.

Give **exactly one** of `thickness` (an absolute height, so the period grows
with the stack) or `fraction` (a share of the total, i.e. of the *final*
period — [`validate_laminate`](@ref) then checks `Σ f ≈ 1` rather than
silently rescaling). Mixing the two forms across layers is rejected: with
imperfect interfaces the absolute period is physically meaningful, so a
half-specified stack would be ambiguous.

`interface` is the condition **on top of** this layer; the last layer's
interface closes the cell onto the first by periodicity. The five types of
`LayeredSpheres` are reused unchanged: [`PerfectInterface`](@ref),
[`SpringInterface`](@ref), [`MembraneInterface`](@ref),
[`KapitzaInterface`](@ref), [`SurfaceConductiveInterface`](@ref).
"""
function add_layer!(
        lam::Laminate{T}, name::Symbol, properties::AbstractDict;
        thickness = nothing, fraction = nothing,
        interface::AbstractInterface = PerfectInterface()
    ) where {T}
    haskey(lam.layers, name) &&
        throw(ArgumentError("layer :$(name) is already registered in this Laminate"))
    (thickness === nothing) == (fraction === nothing) && throw(
        ArgumentError(
            "add_layer!(:$(name)): give exactly one of `thickness = …` or `fraction = …`"
        )
    )
    h = thickness !== nothing ? thickness : fraction
    if h isa Real && h < 0
        throw(ArgumentError("layer :$(name) has negative thickness $(h)"))
    end
    lam.layers[name] = Layer(properties)
    push!(lam.layer_names, name)
    lam.thicknesses[name] = MFH_Core._cell_promote(T, h)
    push!(lam.interfaces, interface)
    lam.period = _compute_period(lam)
    return lam
end

# Total period `L = Σ h_i`, seeded with `zero(T)` so that the result carries
# the effective element type (a `Dual` as soon as one thickness is a `Dual`)
# and, for a symbolic laminate, stays `h₁ + h₂` rather than `0.0 + h₁ + h₂`.
function _compute_period(lam::Laminate{T}) where {T}
    L = zero(T)
    for name in lam.layer_names
        L = L + lam.thicknesses[name]
    end
    return L
end

# =============================================================================
#  Accessors
# =============================================================================

"""
    layer_names(lam) -> Vector{Symbol}

Layer names in stacking order.
"""
layer_names(lam::Laminate) = lam.layer_names

"""
    layer_count(lam) -> Int

Number of layers in one period.
"""
layer_count(lam::Laminate) = length(lam.layer_names)

"""
    layer_property(lam, name::Symbol, key::Symbol)

Property `key` (`:C`, `:K`, …) of layer `name`.

A stored [`Homogenized`](@ref) is resolved here — memoized for the duration
of the enclosing `homogenize` call — so the kernel always sees a plain
tensor. Use [`layer_property_raw`](@ref) to inspect the stored value.
"""
function layer_property(lam::Laminate, name::Symbol, key::Symbol)
    haskey(lam.layers, name) ||
        throw(ArgumentError("no layer named :$(name) in Laminate"))
    l = lam.layers[name]
    haskey(l.properties, key) ||
        throw(ArgumentError("layer :$(name) does not carry property :$(key)"))
    return MFH_Core.resolve_property(l.properties[key], key)
end

"""
    layer_property_raw(lam, name::Symbol, key::Symbol)

The value stored under `key` on layer `name`, **without** resolving a
[`Homogenized`](@ref).
"""
function layer_property_raw(lam::Laminate, name::Symbol, key::Symbol)
    haskey(lam.layers, name) ||
        throw(ArgumentError("no layer named :$(name) in Laminate"))
    l = lam.layers[name]
    haskey(l.properties, key) ||
        throw(ArgumentError("layer :$(name) does not carry property :$(key)"))
    return l.properties[key]
end

"""
    layer_thickness(lam, name::Symbol) -> Number

Thickness `h_i` of layer `name`.
"""
function layer_thickness(lam::Laminate, name::Symbol)
    haskey(lam.thicknesses, name) ||
        throw(ArgumentError("no layer named :$(name) in Laminate"))
    return lam.thicknesses[name]
end

"""
    layer_volume_fraction(lam, name::Symbol) -> Number

Volume fraction `f_i = h_i / L` of layer `name`.
"""
layer_volume_fraction(lam::Laminate, name::Symbol) =
    layer_thickness(lam, name) / laminate_period(lam)

"""
    laminate_period(lam) -> Number

Period `L = Σ h_i` of the cell. Imperfect interfaces enter the effective
property with weight `1/L`, so this is what sets their size effect; with
perfect interfaces the result is independent of it.
"""
laminate_period(lam::Laminate) = lam.period

"""
    laminate_basis(lam) -> AbstractBasis

The orthonormal frame `(ℓ, m, n)` of the cell; its third axis is the layer
normal.
"""
laminate_basis(lam::Laminate) = lam.basis

"""
    laminate_normal(lam) -> NTuple{3}

The layer normal `n`, i.e. the third axis of [`laminate_basis`](@ref),
in canonical components.
"""
laminate_normal(lam::Laminate) = MFH_Core._basis_col(lam.basis, 3)

"""
    layer_interface(lam, k::Integer) -> AbstractInterface

The interface on top of layer `k` (the `k`-th in stacking order); `k = N`
closes the cell onto layer 1 by periodicity.
"""
function layer_interface(lam::Laminate, k::Integer)
    1 ≤ k ≤ length(lam.interfaces) ||
        throw(BoundsError("laminate has $(length(lam.interfaces)) interfaces; asked for $(k)"))
    return lam.interfaces[k]
end

"""
    eltype(lam::Laminate) -> Type

Effective thickness element type: the promotion of the declared floor with
every stored thickness. Use `eltype(typeof(lam))` for the declared floor.
"""
Base.eltype(lam::Laminate{T}) where {T} =
    mapfoldl(typeof, promote_type, values(lam.thicknesses); init = T)

Base.eltype(::Type{<:Laminate{T}}) where {T} = T

# =============================================================================
#  Validation
# =============================================================================

"""
    validate_laminate(lam)

Sanity-check the cell: at least one layer, one interface per layer,
non-negative thicknesses, a strictly positive period, and volume fractions
summing to one.

The requirement that `validate_rve` imposes and this one cannot — a
registered matrix phase — is exactly why [`validate_cell`](@ref) exists.
Non-`Real` thicknesses (symbolic, `Dual`) skip the inequality checks, as
elsewhere in the package.
"""
function validate_laminate(lam::Laminate)
    isempty(lam.layer_names) &&
        throw(ArgumentError("Laminate has no layer; call add_layer! first"))
    length(lam.interfaces) == length(lam.layer_names) || throw(
        ArgumentError(
            "Laminate has $(length(lam.interfaces)) interfaces for " *
                "$(length(lam.layer_names)) layers (one per layer, periodic)"
        )
    )
    for name in lam.layer_names
        h = lam.thicknesses[name]
        if h isa Real && h < 0
            throw(ArgumentError("layer :$(name) has negative thickness $(h)"))
        end
    end
    L = laminate_period(lam)
    if L isa Real && L ≤ 0
        throw(ArgumentError("Laminate period must be > 0; got $(L)"))
    end
    return lam
end

MFH_Core.validate_cell(lam::Laminate) = validate_laminate(lam)

# =============================================================================
#  The `AbstractHomogenizationCell` contract
# =============================================================================

MFH_Core.cell_member_names(lam::Laminate) = lam.layer_names

MFH_Core.cell_container_property(lam::Laminate, name::Symbol, key::Symbol) =
    layer_property_raw(lam, name, key)

function MFH_Core.cell_set_property(lam::Laminate, name::Symbol, key::Symbol, value)
    haskey(lam.layers, name) ||
        throw(ArgumentError("no layer named :$(name) in Laminate"))
    l = lam.layers[name]
    new_props = Dict{Symbol, Any}(l.properties)
    new_props[key] = value
    new_layers = Dict{Symbol, Layer}(lam.layers)
    new_layers[name] = Layer(new_props)
    return _rebuild_laminate(lam; layers = new_layers)
end

"""
    _rebuild_laminate(lam; layers, thicknesses, interfaces, T) -> Laminate

Immutable rebuild used by every `set_param` path: a fresh cell sharing the
untouched fields with `lam`, never a mutation. Mirrors `_rebuild_rve`.
"""
function _rebuild_laminate(
        lam::Laminate;
        layers = lam.layers,
        thicknesses = lam.thicknesses,
        interfaces = lam.interfaces,
        T = eltype(typeof(lam))
    )
    new_lam = Laminate{T, typeof(lam.basis)}(
        copy(lam.layer_names),
        layers isa Dict{Symbol, Layer} ? layers : Dict{Symbol, Layer}(layers),
        thicknesses isa Dict{Symbol, Any} ? thicknesses : Dict{Symbol, Any}(thicknesses),
        interfaces isa Vector{AbstractInterface} ? interfaces :
            Vector{AbstractInterface}(interfaces),
        lam.basis,
        zero(T),
    )
    # Property / interface lenses leave the thicknesses alone, which is the
    # common case under autodiff — carry the cached period over instead of
    # summing the dict again.
    new_lam.period = thicknesses === lam.thicknesses && T === eltype(typeof(lam)) ?
        lam.period : _compute_period(new_lam)
    return new_lam
end

# =============================================================================
#  Pretty printing
# =============================================================================

function Base.show(io::IO, ::MIME"text/plain", lam::Laminate{T}) where {T}
    Te = eltype(lam)
    tag = Te === T ? "$T" : "$T → $Te"
    n = laminate_normal(lam)
    println(io, "Laminate{$tag} with ", layer_count(lam), " layer(s)")
    println(io, "  normal : (", join(map(x -> round(float(x); digits = 4), n), ", "), ")")
    println(io, "  period : ", laminate_period(lam))
    for (k, name) in enumerate(lam.layer_names)
        h = lam.thicknesses[name]
        println(io, "  layer  : :$(name)   h = $(h)")
        itf = lam.interfaces[k]
        itf isa PerfectInterface || println(io, "      interface above : ", itf)
    end
    return
end

Base.show(io::IO, lam::Laminate{T}) where {T} =
    print(io, "Laminate{$T}($(layer_count(lam)) layers, L = $(laminate_period(lam)))")
