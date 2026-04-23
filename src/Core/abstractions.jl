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
abstract type AbstractInclusion{T <: Number} end

"""
    AbstractEllipsoidalInclusion{dim,T} <: AbstractInclusion{T}

Supertype for ellipsoidal inclusions — solid ellipsoids (and their
degenerate limits: spheres, cylinders, discs …).  The first type
parameter `dim` encodes the spatial dimension (2 or 3).
"""
abstract type AbstractEllipsoidalInclusion{dim, T} <: AbstractInclusion{T} end

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
abstract type AbstractLayeredInclusion{dim, T} <: AbstractInclusion{T} end

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

"""
    shape_tensor(incl::AbstractInclusion) -> AbstractTens{2}

Symmetric 2nd-order tensor encoding both the semi-axes and the
orientation of the inclusion in the global (canonical) frame:

```math
\\mathbf A = \\mathbf R \\; \\mathrm{diag}(a_1, a_2, \\dots) \\; \\mathbf R^{\\!T}
```

where ``\\mathbf R`` is the rotation matrix mapping the canonical frame
onto the inclusion's local basis and the diagonal entries are the
semi-axes in the order dictated by the local basis.

Conventions for degenerate cases:

| Inclusion               | Diagonal (principal frame)      |
| ----------------------- | ------------------------------- |
| `Ellipsoid{3}`          | `(a₁, a₂, a₃)`                  |
| `Ellipsoid{2}`          | `(a₁, a₂)`                      |
| `Cylinder`              | `(Inf, b, c)` — axis ``e_1``    |
| `EllipticCrack`         | `(a, b, 0)`  — normal ``e_3``   |
| `RibbonCrack`           | `(Inf, b, 0)`                   |
"""
function shape_tensor end

"""
    eshelby_tensor(incl, C₀; method=:auto, abstol, reltol, maxiters) -> AbstractTens

Eshelby tensor of the inclusion `incl` embedded in a matrix of
stiffness / conductivity `C₀`, derived from the Hill polarisation
tensor ``\\mathbb P`` (or ``\\mathbf P``) by the relations

```math
\\mathbb S = \\mathbb P : \\mathbb C_0
\\qquad\\text{(order 4, elasticity)}
```

```math
\\mathbf s = \\mathbf P \\cdot \\mathbf K_0
\\qquad\\text{(order 2, conductivity / diffusion)}
```

The appropriate method is selected by dispatch on the order of `C₀`:
an `AbstractTens{4, 3}` (elasticity) triggers the double contraction
``\\mathbb P \\;\\underset{s}{:}\\; \\mathbb C_0``, while an
`AbstractTens{2, 3}` (conductivity) triggers the simple contraction
``\\mathbf P \\cdot \\mathbf K_0``.

All keyword arguments (`method`, `abstol`, `reltol`, `maxiters`) are
forwarded verbatim to [`hill_tensor`](@ref); see its docstring for the
set of admissible algorithm traits.

See also [`hill_tensor`](@ref).
"""
function eshelby_tensor end

# =============================================================================
#  Localization & contribution public API stubs
#
#  Declared here so that every sub-module (Elasticity, Cracks,
#  Conductivity, LayeredSphere, user extensions) can add methods via
#  `import ..Core: stiffness_contribution` + method definition, all
#  attaching to a single generic function.  Definitions live in
#  `src/localization.jl` and `src/contribution.jl` (loaded at
#  MeanFieldHom top level after all sub-modules).
# =============================================================================

"""
    strain_strain_loc(incl, C₁, C₀; kw...)  -> Tens{4,3}

Dilute strain-strain localization tensor `A_εε` (Eshelby).
"""
function strain_strain_loc end

"""
    stress_strain_loc(incl, C₁, C₀; kw...)  -> Tens{4,3}
"""
function stress_strain_loc end

"""
    strain_stress_loc(incl, C₁, C₀; kw...)  -> Tens{4,3}
"""
function strain_stress_loc end

"""
    stress_stress_loc(incl, C₁, C₀; kw...)  -> Tens{4,3}
"""
function stress_stress_loc end

"""
    gradient_gradient_loc(incl, K₁, K₀; kw...) -> Tens{2,3}
"""
function gradient_gradient_loc end

"""
    flux_gradient_loc(incl, K₁, K₀; kw...)    -> Tens{2,3}
"""
function flux_gradient_loc end

"""
    gradient_flux_loc(incl, K₁, K₀; kw...)    -> Tens{2,3}
"""
function gradient_flux_loc end

"""
    flux_flux_loc(incl, K₁, K₀; kw...)        -> Tens{2,3}
"""
function flux_flux_loc end

"""
    stiffness_contribution(incl, C₁, C₀; kw...) -> Tens{4,3}
    stiffness_contribution(crack, C₀; kw...)   -> Tens{4,3}

Size-independent stiffness contribution tensor `N` of an inclusion in
a matrix `C₀`.  For a dilute family of volume fraction `f`:
`ΔC_eff = f · N` (see [`delta_stiffness`](@ref)).
"""
function stiffness_contribution end

"""
    conductivity_contribution(incl, K₁, K₀; kw...) -> Tens{2,3}
    conductivity_contribution(crack, K₀; kw...)     -> Tens{2,3}

Size-independent conductivity contribution tensor for the 2nd-order
transport problem.  Analogue of [`stiffness_contribution`](@ref).
"""
function conductivity_contribution end

"""
    resistivity_contribution(incl, K₁, K₀; kw...) -> Tens{2,3}

Size-independent resistivity contribution tensor of an inclusion
(2nd-order analogue of [`compliance_contribution`](@ref) for solid
ellipsoids).
"""
function resistivity_contribution end

"""
    delta_stiffness(N, f) -> Tens{4,3}

Dilute effective-stiffness correction `ΔC = f · N` from the size-
independent contribution tensor `N` and the volume fraction `f`.
"""
function delta_stiffness end

"""
    delta_conductivity(N_K, f) -> Tens{2,3}

Dilute effective-conductivity correction `ΔK = f · N_K`.
"""
function delta_conductivity end
