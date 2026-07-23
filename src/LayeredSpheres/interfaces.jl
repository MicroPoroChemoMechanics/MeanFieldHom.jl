# =============================================================================
#  interfaces.jl — imperfect interface models for `LayeredSphere`.
#
#  Four physically-motivated interface types are provided, organized as
#  a "primal / dual" pair per physics:
#
#   Elasticity
#   ----------
#   - `SpringInterface(kn, kt)` — displacement jump (primal):
#     `[u_n] = kn · t_n`, `[u_t] = kt · t_t`.
#   - `MembraneInterface(κs, μs)` — traction jump (dual, surface
#     elasticity): surface stiffness introduces a jump in the normal
#     component of the traction proportional to the surface strain.
#
#   Conductivity (thermal / electric / Darcy)
#   -----------------------------------------
#   - `KapitzaInterface(ρ)` — temperature jump (primal, interfacial
#     thermal resistance): `[T] = ρ · q_n`.
#   - `SurfaceConductiveInterface(ks)` — flux jump (dual, highly
#     conductive 2D layer): `[q_n] = -divₛ(ks ∇ₛ T)`.
#
#  `PerfectInterface` is the trivial limit of any of them (k→0 for the
#  primal types, ks→0 / κs=μs=0 for the dual types).
# =============================================================================

"""
    AbstractInterface{T}

Root supertype for interface conditions in a `LayeredSphere`.  Concrete
subtypes determine the jump matrix applied to the state vector
`(u_r, σ_rr)` (bulk), `(U, V, σ_rr, σ_rθ)` (shear), or `(T, q_n)`
(conductivity).
"""
abstract type AbstractInterface{T <: Number} end

"""
    PerfectInterface{T}()

Perfect (continuous) interface: all state-vector components are
continuous.
"""
struct PerfectInterface{T <: Number} <: AbstractInterface{T} end

PerfectInterface() = PerfectInterface{Float64}()

"""
    SpringInterface{T}(kn::T, kt::T)

Imperfect interface of "spring" type with two compliances — normal
`kn` and tangential `kt`:

```
[u_n] = kn · t_n,        [u_t] = kt · t_t,       t_n, t_t continuous.
```

The `kn = kt = 0` limit recovers [`PerfectInterface`](@ref); the
`kn, kt → ∞` limit is a free-surface (fully decoupled layer boundary).
"""
struct SpringInterface{T <: Number} <: AbstractInterface{T}
    kn::T
    kt::T
end

# Convenience: normal-only compliance (kt = 0 ≡ no tangential jump).
SpringInterface(k::Number) = SpringInterface(k, zero(k))

"""
    MembraneInterface{T}(κs::T, μs::T)

Imperfect interface of surface-elastic (Gurtin–Murdoch "membrane") type —
the dual analog of [`SpringInterface`](@ref) and the elastic counterpart of
Echoes' `DUALDISC`.  The interface behaves as a 2D elastic shell with
surface moduli `κs = λs + μs` (surface dilatation, matching Echoes' `ks`)
and surface shear `μs`.  Displacement is continuous across the interface
and the surface strain generates a traction jump (`[σ·n] = −divₛσˢ`).  On a
spherical interface of radius `r`, the bulk (`Y₀`) mode jump is

```
[σ_rr] = (4 κs / r²) · u_r,
```

and the shear (`Y₂`-harmonic) mode jump, with `u_r = U P₂`,
`u_θ = W dP₂/dθ`, is

```
[σ_rr] = ( 4κs U − 12κs W) / r²,
[σ_rθ] = (−2κs U + (6κs + 4μs) W) / r².
```

The `κs = μs = 0` limit recovers [`PerfectInterface`](@ref).  These jumps
reproduce Echoes' `DUALDISC` concentration tensors and effective moduli to
machine precision.
"""
struct MembraneInterface{T <: Number} <: AbstractInterface{T}
    κs::T
    μs::T
end

"""
    KapitzaInterface{T}(resistance::T)

Thermal imperfect interface with scalar thermal resistance:
`[T] = resistance · q_n`, with `q_n` continuous.  Primal analog of
[`SpringInterface`](@ref).
"""
struct KapitzaInterface{T <: Number} <: AbstractInterface{T}
    resistance::T
end

"""
    SurfaceConductiveInterface{T}(conductance::T)

Highly-conductive 2D surface layer (dual analog of
[`MembraneInterface`](@ref)).  Introduces a flux jump driven by the
surface Laplacian of the temperature; for the spherical harmonic `Y_n`
on a spherical interface of radius `r`,

```
[q_n] = -n(n+1) · conductance · T / r².
```

`conductance = 0` recovers [`PerfectInterface`](@ref).
"""
struct SurfaceConductiveInterface{T <: Number} <: AbstractInterface{T}
    conductance::T
end

Base.eltype(::AbstractInterface{T}) where {T} = T
Base.eltype(::Type{<:AbstractInterface{T}}) where {T} = T
