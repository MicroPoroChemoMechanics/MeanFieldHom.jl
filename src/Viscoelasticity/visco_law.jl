# =============================================================================
#  visco_law.jl — viscoelastic kernel `R(t,t')` or `J(t,t')`.
#
#  A `ViscoLaw` wraps an `eval` callable taking two times and returning
#  either a scalar (`Real` or `Complex`), a 4-tensor (`TensND.AbstractTens{4,3}`),
#  or a `6×6` Mandel matrix.  The `mode` field selects whether the kernel
#  is a relaxation (`R(t,t') = σ / ε` for a unit strain step at `t'`) or
#  a creep (`J(t,t')` for a unit stress step) law.
# =============================================================================

"""
    AbstractViscoLaw

Root supertype for viscoelastic kernels.  Concrete subtypes provide an
`eval(law, t, t_p)` method returning a scalar or a 4-tensor.
"""
abstract type AbstractViscoLaw end

"""
    VALID_VISCO_MODES :: Tuple{Symbol, Symbol}

Allowed values for `ViscoLaw.mode` :

  * `:relaxation` — `R(t, t')` maps strain history to stress.
  * `:creep`      — `J(t, t')` maps stress history to strain.
"""
const VALID_VISCO_MODES = (:relaxation, :creep)

"""
    ViscoLaw(eval_fun, mode::Symbol = :relaxation)

Concrete viscoelastic law.  Wraps an `eval_fun::F` callable
`(t, t_p) -> X` together with its `mode` Symbol (`:relaxation` or
`:creep`).  The output type `X` may be:

  * a `Real` or `Complex` scalar (scalar Volterra kernel) ;
  * a `TensND.AbstractTens{4,3}` (4-tensor relaxation / creep tensor) ;
  * an `AbstractMatrix` of size `6×6` already in Mandel form.

Construct via [`ViscoLaw`](@ref) directly or use the convenience
constructors [`maxwell_relaxation`](@ref), [`kelvin_creep`](@ref),
[`maxwell_iso`](@ref), [`kelvin_iso`](@ref).
"""
struct ViscoLaw{F} <: AbstractViscoLaw
    eval_fun::F
    mode::Symbol

    function ViscoLaw(eval_fun::F, mode::Symbol = :relaxation) where {F}
        mode in VALID_VISCO_MODES ||
            throw(ArgumentError(
                "ViscoLaw mode must be one of $(VALID_VISCO_MODES); got :$(mode)"
            ))
        return new{F}(eval_fun, mode)
    end
end

# Direct call syntax: `law(t, t_p)` evaluates the kernel.
(law::ViscoLaw)(t, t_p) = law.eval_fun(t, t_p)

"""
    visco_mode(law) -> Symbol

Return the mode (`:relaxation` or `:creep`) of `law`.
"""
visco_mode(law::ViscoLaw) = law.mode

"""
    visco_eval(law, t, t_p)

Evaluate the kernel at `(t, t_p)`.
"""
visco_eval(law::ViscoLaw, t, t_p) = law.eval_fun(t, t_p)

# ── Convenience constructors ────────────────────────────────────────────────

"""
    heaviside_law(C; mode = :relaxation)

Return a `ViscoLaw` corresponding to a purely elastic kernel
`R(t,t') = C · H(t-t')`.  Useful for testing the elastic limit of the
ALV pipeline.
"""
function heaviside_law(C; mode::Symbol = :relaxation)
    eval_fun = (t, t_p) -> (t ≥ t_p) ? C : zero(C)
    return ViscoLaw(eval_fun, mode)
end

"""
    maxwell_relaxation(C_inf, C_branches, taus; mode = :relaxation)

Build a generalised Maxwell relaxation kernel

```
R(t, t') = C_inf + Σ_i C_branches[i] · exp(-(t - t')/taus[i])
```

for `t ≥ t'`, `0` otherwise.  `C_inf` and `C_branches[i]` may be
scalars or 4-tensors (`TensND.AbstractTens{4,3}`); `taus[i]` is a
positive relaxation time of the `i`-th branch.
"""
function maxwell_relaxation(C_inf, C_branches::AbstractVector, taus::AbstractVector;
                            mode::Symbol = :relaxation)
    length(C_branches) == length(taus) ||
        throw(ArgumentError("maxwell_relaxation: C_branches and taus must have the same length"))
    eval_fun = function (t, t_p)
        if t < t_p
            return zero(C_inf)
        end
        result = deepcopy(C_inf)
        Δt = t - t_p
        @inbounds for i in eachindex(C_branches)
            result = result + C_branches[i] * exp(-Δt / taus[i])
        end
        return result
    end
    return ViscoLaw(eval_fun, mode)
end

"""
    kelvin_creep(J_0, J_branches, taus; mode = :creep)

Build a Kelvin (or Kelvin-Voigt-Generalised) creep kernel

```
J(t, t') = J_0 + Σ_i J_branches[i] · (1 - exp(-(t - t')/taus[i]))
```

for `t ≥ t'`, `0` otherwise.  `J_0` is the instantaneous compliance,
`J_branches[i]` the `i`-th branch compliance, `taus[i]` its retardation
time.
"""
function kelvin_creep(J_0, J_branches::AbstractVector, taus::AbstractVector;
                      mode::Symbol = :creep)
    length(J_branches) == length(taus) ||
        throw(ArgumentError("kelvin_creep: J_branches and taus must have the same length"))
    eval_fun = function (t, t_p)
        if t < t_p
            return zero(J_0)
        end
        result = deepcopy(J_0)
        Δt = t - t_p
        @inbounds for i in eachindex(J_branches)
            result = result + J_branches[i] * (one(eltype_data(J_branches[i])) - exp(-Δt / taus[i]))
        end
        return result
    end
    return ViscoLaw(eval_fun, mode)
end

# Helper for kelvin_creep — extract a scalar element type to write `1 - exp(...)`
# in a way compatible with both scalar and tensor branches.
@inline eltype_data(x::Number) = typeof(x)
@inline eltype_data(x::TensND.AbstractTens) = eltype(x)
@inline eltype_data(x::AbstractArray) = eltype(x)

"""
    maxwell_iso(k, mu, eta_k, eta_mu) -> ViscoLaw

Convenience: build an isotropic Maxwell **relaxation** 4-tensor kernel

```
R(t, t') = 3 k · exp(-(t - t')/eta_k) · 𝕁
         + 2 mu · exp(-(t - t')/eta_mu) · 𝕂
```

with `𝕁` the spherical projector `(1/3) 𝟙 ⊗ 𝟙` and `𝕂 = 𝕀 - 𝕁` the
deviatoric projector.  The output is a `TensND.TensISO{4,3}` at every
`(t, t')` with `t ≥ t'`.
"""
function maxwell_iso(k, mu, eta_k, eta_mu)
    eval_fun = function (t, t_p)
        if t < t_p
            T = promote_type(typeof(k), typeof(mu), typeof(eta_k), typeof(eta_mu),
                             typeof(t), typeof(t_p))
            return TensND.TensISO{3}(zero(T), zero(T))
        end
        Δt = t - t_p
        α = 3 * k * exp(-Δt / eta_k)
        β = 2 * mu * exp(-Δt / eta_mu)
        return TensISO{3}(α, β)
    end
    return ViscoLaw(eval_fun, :relaxation)
end

"""
    kelvin_iso(k_0, mu_0, k_branches, mu_branches, taus_k, taus_mu) -> ViscoLaw

Convenience: build an isotropic Kelvin **creep** 4-tensor kernel

```
J(t, t') = (1/(3 k_0)) 𝕁 + (1/(2 mu_0)) 𝕂
         + Σ_i (1/(3 k_branches[i])) (1 - exp(-(t-t')/taus_k[i])) 𝕁
         + Σ_i (1/(2 mu_branches[i])) (1 - exp(-(t-t')/taus_mu[i])) 𝕂
```

`k_branches`, `mu_branches`, `taus_k`, `taus_mu` may be empty if no
Kelvin branches are required (instantaneous-only compliance).
"""
function kelvin_iso(k_0, mu_0,
                    k_branches::AbstractVector = Float64[],
                    mu_branches::AbstractVector = Float64[],
                    taus_k::AbstractVector = Float64[],
                    taus_mu::AbstractVector = Float64[])
    length(k_branches) == length(taus_k) ||
        throw(ArgumentError("kelvin_iso: k_branches and taus_k length mismatch"))
    length(mu_branches) == length(taus_mu) ||
        throw(ArgumentError("kelvin_iso: mu_branches and taus_mu length mismatch"))
    eval_fun = function (t, t_p)
        T = promote_type(typeof(k_0), typeof(mu_0), typeof(t), typeof(t_p))
        if t < t_p
            return TensND.TensISO{3}(zero(T), zero(T))
        end
        α = T(1) / k_0    # 3 (1/(3 k_0)) for the J-axis (TensISO data is (3K, 2μ))
        β = T(1) / mu_0   # 2 (1/(2 mu_0)) for the K-axis
        Δt = t - t_p
        @inbounds for i in eachindex(k_branches)
            α += (T(1) / k_branches[i]) * (one(T) - exp(-Δt / taus_k[i]))
        end
        @inbounds for i in eachindex(mu_branches)
            β += (T(1) / mu_branches[i]) * (one(T) - exp(-Δt / taus_mu[i]))
        end
        return TensISO{3}(α, β)
    end
    return ViscoLaw(eval_fun, :creep)
end

# ── Pretty printing ──────────────────────────────────────────────────────────

function Base.show(io::IO, law::ViscoLaw)
    return print(io, "ViscoLaw(:", law.mode, ")")
end
