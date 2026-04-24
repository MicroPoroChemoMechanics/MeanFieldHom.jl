# =============================================================================
#  sif.jl — stress / displacement intensity factors.
# =============================================================================

"""
    sif(crack, C₀, Σ; y₀=nothing, method=:auto, kw...) -> (𝐊, (Kᴵ, Kᴵᴵ, Kᴵᴵᴵ))

Stress intensity factor vector ``\\hat{\\mathbf K}`` at a point of the
crack front, together with its ``(K_{I},K_{II},K_{III})`` decomposition
on ``(\\hat{\\mathbf n},\\hat{\\boldsymbol\\nu},\\hat{\\boldsymbol\\tau})``
([Irwin 1957](@cite irwin1957),
 [Kassir & Sih 1968](@cite kassir1968),
 [Willis 1968](@cite willis1968);
 energy release rate identity ``G = \\hat{\\mathbf K}\\cdot\\hat{\\mathbf N}``
 in [Barnett & Asaro 1972](@cite barnett1972),
 [Rice 1989](@cite rice1989)).

For a ribbon crack (``\\hat{\\boldsymbol\\nu}=\\pm\\hat{\\mathbf m}``)
``\\hat{\\mathbf K}^{\\mathcal R} = \\sqrt{\\pi b}\\,\\boldsymbol\\Sigma\\cdot\\hat{\\mathbf n}``
(independent of the matrix stiffness).
For an elliptic crack, ``\\hat{\\mathbf K}`` is obtained from the COD
tensor ``\\mathbf B^{\\mathcal E}`` of the actual crack and the COD
tensor ``\\mathbf B^{\\mathcal R}`` of the tangent ribbon crack at the
observation point:

```
K̂ = (3/8) π^{3/2} √b √(b ‖S† · ŷ₀★‖)
    · (B^𝓡(ν̂, n̂))⁻¹ · B^𝓔(m̂, n̂, η) · Σ·n̂ .
```

The central identity
``\\hat{\\mathbf K} = \\pi\\,(\\mathbf B^{\\mathcal R})^{-1}\\cdot\\hat{\\mathbf N}``
is purely local
([Kanaun 1981](@cite kanaun1981), [Kunin 1983](@cite kunin1983),
 [Kanaun & Levin 2009](@cite kanaun2009)).
"""
function sif end

# Ribbon crack (paper eq. 737)
function sif(
        crack::RibbonCrack{T}, C₀, Σ;
        y₀ = nothing, method::Symbol = :auto, kw...
    ) where {T}
    b = crack.b
    l̂, m̂, n̂ = (TensND.tens_basis(crack_basis(crack), i) for i in 1:3)
    𝐊 = sqrt(T(π) * b) * (Σ ⋅ n̂)
    Kᴵ = 𝐊 ⋅ n̂
    Kᴵᴵ = 𝐊 ⋅ m̂
    Kᴵᴵᴵ = 𝐊 ⋅ l̂
    return 𝐊, (Kᴵ, Kᴵᴵ, Kᴵᴵᴵ)
end

# Elliptical crack (paper eq. 719)
function sif(
        crack::EllipticCrack{T}, C₀, Σ;
        y₀ = nothing, method::Symbol = :auto, kw...
    ) where {T}
    a = crack.a
    b = crack.b
    ℬ = crack_basis(crack)
    l̂, m̂, n̂ = (TensND.tens_basis(ℬ, i) for i in 1:3)

    𝐒_inv = inv(a) * (l̂ ⊗ l̂) + inv(b) * (m̂ ⊗ m̂)

    y0 = y₀ === nothing ? m̂ : y₀
    S⁻¹_y0 = 𝐒_inv ⋅ y0
    n_Sy = norm(S⁻¹_y0)
    ν̂ = TensND.change_tens(S⁻¹_y0 / n_Sy, ℬ)
    τ̂ = TensND.Tens(TensND.change_tens(n̂, ℬ) × TensND.change_tens(ν̂, ℬ), ℬ)

    ℬ_ν = TensND.Basis(
        hcat(
            TensND.components_canon(τ̂),
            TensND.components_canon(ν̂),
            TensND.components_canon(n̂)
        )
    )

    B_ℰ = cod_tensor(crack, C₀; method = method, kw...)
    ribbon_ref = RibbonCrack(b, ℬ_ν)
    B_ℛ = cod_tensor(ribbon_ref, C₀; method = method, kw...)

    𝐊 = (3 * T(π)^(T(3) / 2) * b / 8) * sqrt(b * n_Sy) *
        inv(B_ℛ) ⋅ B_ℰ ⋅ Σ ⋅ n̂

    Kᴵ = 𝐊 ⋅ n̂
    Kᴵᴵ = 𝐊 ⋅ ν̂
    Kᴵᴵᴵ = 𝐊 ⋅ τ̂
    return 𝐊, (Kᴵ, Kᴵᴵ, Kᴵᴵᴵ)
end

# =============================================================================
#  Displacement intensity factor
# =============================================================================

"""
    dif(crack, C₀, Σ; method=:auto, kw...) -> Tens{1,3}
"""
function dif(crack::MFH_Core.AbstractCrack, C₀, Σ; method::Symbol = :auto, kw...)
    B = cod_tensor(crack, C₀; method = method, kw...)
    n̂ = TensND.tens_basis(crack_basis(crack), 3)
    return B ⋅ Σ ⋅ n̂
end

# =============================================================================
#  Thermal (2nd-order) — heat-flux / temperature intensity factors.
# =============================================================================

# Heat-flux intensity factor `K_T` — thermal analogue of the elasticity SIF.
# For a flat crack driven by remote heat flux ``\mathbf q^∞``, the crack-tip
# singular field scales as ``\sim K_T/\sqrt{r}``:
#   - Ribbon:    ``K_T = \sqrt{\pi b}\,(\hat n\cdot \mathbf q^∞)``
#   - Elliptic:  ``K_T = (3\pi^{3/2} b/8)\sqrt{b\,n_S}\,(b^{\mathcal E}/b^{\mathcal R})
#                      \,(\hat n\cdot \mathbf q^∞)``
# Only the mode I analogue exists in the scalar-temperature case (no shear mode).

"""
    sif(crack::RibbonCrack, K₀::AbstractTens{2,3}, q∞; kw...) -> Real

Thermal SIF (heat-flux intensity factor) of a ribbon crack:
``K_T = \\sqrt{\\pi b}\\;\\hat{\\mathbf n}\\cdot\\mathbf q^{\\infty}``.
"""
function sif(
        crack::RibbonCrack{T},
        K₀::TensND.AbstractTens{2, 3},
        q∞;
        y₀ = nothing, method::Symbol = :auto, kw...
    ) where {T}
    b = crack.b
    n̂ = TensND.tens_basis(crack_basis(crack), 3)
    return sqrt(T(π) * b) * (n̂ ⋅ q∞)
end

# Elliptic crack — thermal
function sif(
        crack::EllipticCrack{T},
        K₀::TensND.AbstractTens{2, 3},
        q∞;
        y₀ = nothing, method::Symbol = :auto, kw...
    ) where {T}
    a = crack.a
    b = crack.b
    ℬ = crack_basis(crack)
    l̂, m̂, n̂ = (TensND.tens_basis(ℬ, i) for i in 1:3)

    𝐒_inv = inv(a) * (l̂ ⊗ l̂) + inv(b) * (m̂ ⊗ m̂)

    y0 = y₀ === nothing ? m̂ : y₀
    S⁻¹_y0 = 𝐒_inv ⋅ y0
    n_Sy = norm(S⁻¹_y0)
    ν̂ = TensND.change_tens(S⁻¹_y0 / n_Sy, ℬ)
    τ̂ = TensND.Tens(TensND.change_tens(n̂, ℬ) × TensND.change_tens(ν̂, ℬ), ℬ)

    ℬ_ν = TensND.Basis(
        hcat(
            TensND.components_canon(τ̂),
            TensND.components_canon(ν̂),
            TensND.components_canon(n̂)
        )
    )

    b_ℰ = cod_tensor(crack, K₀; method = method, kw...)
    ribbon_ref = RibbonCrack(b, ℬ_ν)
    b_ℛ = cod_tensor(ribbon_ref, K₀; method = method, kw...)

    return (3 * T(π)^(T(3) / 2) * b / 8) * sqrt(b * n_Sy) *
        (b_ℰ / b_ℛ) * (n̂ ⋅ q∞)
end

"""
    dif(crack, K₀::AbstractTens{2,3}, q∞; method=:auto, kw...) -> Real

Temperature intensity factor (analogue of displacement intensity
factor) for a flat crack driven by a remote heat-flux vector
``\\mathbf q^{\\infty}``:

```
[[T]]_avg = b · (\\hat{\\mathbf n}·\\mathbf q^{\\infty}) .
```

Returns a scalar (vs the `Tens{1,3}` returned by the elasticity `dif`)
since the temperature field is scalar.
"""
function dif(
        crack::MFH_Core.AbstractCrack,
        K₀::TensND.AbstractTens{2, 3},
        q∞;
        method::Symbol = :auto, kw...
    )
    b = cod_tensor(crack, K₀; method = method, kw...)
    n̂ = TensND.tens_basis(crack_basis(crack), 3)
    return b * (n̂ ⋅ q∞)
end
