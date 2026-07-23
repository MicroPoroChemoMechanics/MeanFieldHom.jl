# =============================================================================
#  localization.jl — generic localization tensors of the dilute Eshelby
#  problem (Kachanov-Sevostianov / Barthélémy et al. 2021).
#
#  Four tensors relate the fields within an inclusion of stiffness `C₁`
#  (or conductivity `K₁`) embedded in an infinite matrix `C₀` (`K₀`) to
#  the remote far-field `ε∞` / `σ∞`:
#
#      strain_strain_loc : ε_inc = A_εε : ε∞
#      stress_strain_loc : σ_inc = A_σε : ε∞
#      strain_stress_loc : ε_inc = A_εσ : σ∞
#      stress_stress_loc : σ_inc = A_σσ : σ∞
#
#  The pivot formula (dilute) is
#
#      A_εε = [𝕀 + ℙ(incl, C₀) : (C₁ - C₀)]⁻¹,
#
#  where ℙ = `hill_tensor(incl, C₀)`.  The three remaining tensors are
#  derived algebraically:
#
#      A_σε = C₁ : A_εε
#      A_εσ = A_εε : S₀    (S₀ = C₀⁻¹)
#      A_σσ = C₁ : A_εε : S₀ = A_σε : S₀
#
#  Conductivity analogs (2-tensor fields, 2-tensor Hill/moduli) use
#  the same formulas with `·` in place of `:`.
#
#  Type-genericity: the implementation works for Float64, BigFloat,
#  ForwardDiff.Dual, SymPy.Sym and Symbolics.Num so long as `hill_tensor`
#  does; it relies only on TensND algebra (`+`, `-`, `⊡`, `inv`).
# =============================================================================

"""
    _identity_4sym(::Type{T}) -> TensISO{4,3}

Symmetric 4-tensor identity `𝕀_{ijkl} = ½(δ_{ik}δ_{jl} + δ_{il}δ_{jk})`
in its most compact `TensISO` form (3D).  `𝕀 ⊡ X = X` for any symmetric
`Tens{4,3}`.
"""
_identity_4sym(::Type{T}) where {T <: Number} = TensISO{3}(one(T), one(T))

"""
    _identity_2(::Type{T}) -> TensISO{2,3}

Identity 2-tensor `δ_{ij}` in 3D (`TensISO{2,3}`).  `𝟙 · x = x` for any
`Tens{2,3}` or 3-vector.
"""
_identity_2(::Type{T}) where {T <: Number} = TensISO{3}(one(T))

# =============================================================================
#  Elastic localization (4-tensor fields)
# =============================================================================

"""
    strain_strain_loc(incl, C₁, C₀; kw...) -> Tens{4,3}

Dilute **strain-strain localization tensor** `A_εε`: connects the
average strain in an `AbstractInclusion` of stiffness `C₁` to the remote
strain `ε∞`:

```
ε_inc = A_εε : ε∞,
A_εε  = [𝕀 + ℙ(incl, C₀) : (C₁ - C₀)]⁻¹.
```

Keyword arguments are forwarded to [`hill_tensor`](@ref).

See also [`stress_strain_loc`](@ref), [`strain_stress_loc`](@ref),
[`stress_stress_loc`](@ref).
"""
function strain_strain_loc(
        incl::AbstractInclusion,
        C₁::TensND.AbstractTens{4, 3},
        C₀::TensND.AbstractTens{4, 3};
        kw...
    )
    T = promote_type(eltype(C₁), eltype(C₀))
    P = hill_tensor(incl, C₀; kw...)
    δC = C₁ - C₀
    return inv(_identity_4sym(T) + (P ⊡ δC))
end

"""
    stress_strain_loc(incl, C₁, C₀; kw...) -> Tens{4,3}

Dilute **stress-strain localization tensor** `A_σε = C₁ : A_εε`:
`σ_inc = A_σε : ε∞`.
"""
function stress_strain_loc(
        incl::AbstractInclusion,
        C₁::TensND.AbstractTens{4, 3},
        C₀::TensND.AbstractTens{4, 3};
        kw...
    )
    return C₁ ⊡ strain_strain_loc(incl, C₁, C₀; kw...)
end

"""
    strain_stress_loc(incl, C₁, C₀; kw...) -> Tens{4,3}

Dilute **strain-stress localization tensor** `A_εσ = A_εε : S₀`:
`ε_inc = A_εσ : σ∞`.  `S₀ = C₀⁻¹` is built internally.
"""
function strain_stress_loc(
        incl::AbstractInclusion,
        C₁::TensND.AbstractTens{4, 3},
        C₀::TensND.AbstractTens{4, 3};
        kw...
    )
    return strain_strain_loc(incl, C₁, C₀; kw...) ⊡ inv(C₀)
end

"""
    stress_stress_loc(incl, C₁, C₀; kw...) -> Tens{4,3}

Dilute **stress-stress localization tensor** `A_σσ = C₁ : A_εε : S₀`:
`σ_inc = A_σσ : σ∞`.
"""
function stress_stress_loc(
        incl::AbstractInclusion,
        C₁::TensND.AbstractTens{4, 3},
        C₀::TensND.AbstractTens{4, 3};
        kw...
    )
    return C₁ ⊡ strain_strain_loc(incl, C₁, C₀; kw...) ⊡ inv(C₀)
end

# =============================================================================
#  Conductivity localization (2-tensor fields)
# =============================================================================

"""
    gradient_gradient_loc(incl, K₁, K₀; kw...) -> Tens{2,3}

Dilute **gradient-gradient localization tensor** `A_∇∇` for the 2nd
order transport problem:

```
∇T_inc = A_∇∇ · ∇T∞,
A_∇∇   = [𝟙 + ℙ(incl, K₀) · (K₁ - K₀)]⁻¹.
```

Conductivity analog of [`strain_strain_loc`](@ref).  Keyword arguments
are forwarded to [`hill_tensor`](@ref).
"""
function gradient_gradient_loc(
        incl::AbstractInclusion,
        K₁::TensND.AbstractTens{2, 3},
        K₀::TensND.AbstractTens{2, 3};
        kw...
    )
    T = promote_type(eltype(K₁), eltype(K₀))
    P = hill_tensor(incl, K₀; kw...)
    δK = K₁ - K₀
    return inv(_identity_2(T) + (P ⋅ δK))
end

"""
    flux_gradient_loc(incl, K₁, K₀; kw...) -> Tens{2,3}

Dilute **flux-gradient localization tensor** `A_q∇ = K₁ · A_∇∇`:
`q_inc = A_q∇ · ∇T∞`.
"""
function flux_gradient_loc(
        incl::AbstractInclusion,
        K₁::TensND.AbstractTens{2, 3},
        K₀::TensND.AbstractTens{2, 3};
        kw...
    )
    return K₁ ⋅ gradient_gradient_loc(incl, K₁, K₀; kw...)
end

"""
    gradient_flux_loc(incl, K₁, K₀; kw...) -> Tens{2,3}

Dilute **gradient-flux localization tensor** `A_∇q = A_∇∇ · R₀`
(with `R₀ = K₀⁻¹`): `∇T_inc = A_∇q · q∞`.
"""
function gradient_flux_loc(
        incl::AbstractInclusion,
        K₁::TensND.AbstractTens{2, 3},
        K₀::TensND.AbstractTens{2, 3};
        kw...
    )
    return gradient_gradient_loc(incl, K₁, K₀; kw...) ⋅ inv(K₀)
end

"""
    flux_flux_loc(incl, K₁, K₀; kw...) -> Tens{2,3}

Dilute **flux-flux localization tensor** `A_qq = K₁ · A_∇∇ · R₀`:
`q_inc = A_qq · q∞`.
"""
function flux_flux_loc(
        incl::AbstractInclusion,
        K₁::TensND.AbstractTens{2, 3},
        K₀::TensND.AbstractTens{2, 3};
        kw...
    )
    return K₁ ⋅ gradient_gradient_loc(incl, K₁, K₀; kw...) ⋅ inv(K₀)
end
