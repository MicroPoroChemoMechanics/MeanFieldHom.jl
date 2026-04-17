# =============================================================================
#  sif.jl — stress / displacement intensity factors.
# =============================================================================

"""
    sif(crack, C₀, Σ; y₀=nothing, method=:auto, kw...) -> (𝐊, (Kᴵ, Kᴵᴵ, Kᴵᴵᴵ))
"""
function sif end

# Ribbon crack (paper eq. 737)
function sif(crack::RibbonCrack{T}, C₀, Σ;
             y₀=nothing, method::Symbol=:auto, kw...) where {T}
    b  = crack.b
    l̂, m̂, n̂ = (TensND.tensbasis(crack_basis(crack), i) for i in 1:3)
    𝐊 = sqrt(T(π) * b) * (Σ ⋅ n̂)
    Kᴵ   = 𝐊 ⋅ n̂
    Kᴵᴵ  = 𝐊 ⋅ m̂
    Kᴵᴵᴵ = 𝐊 ⋅ l̂
    return 𝐊, (Kᴵ, Kᴵᴵ, Kᴵᴵᴵ)
end

# Elliptical crack (paper eq. 719)
function sif(crack::EllipticCrack{T}, C₀, Σ;
             y₀=nothing, method::Symbol=:auto, kw...) where {T}
    a = crack.a
    b = crack.b
    ℬ = crack_basis(crack)
    l̂, m̂, n̂ = (TensND.tensbasis(ℬ, i) for i in 1:3)

    𝐒_inv = inv(a) * (l̂ ⊗ l̂) + inv(b) * (m̂ ⊗ m̂)

    y0 = y₀ === nothing ? m̂ : y₀
    S⁻¹_y0 = 𝐒_inv ⋅ y0
    n_Sy = norm(S⁻¹_y0)
    ν̂ = TensND.change_tens(S⁻¹_y0 / n_Sy, ℬ)
    τ̂ = TensND.Tens(TensND.change_tens(n̂, ℬ) × TensND.change_tens(ν̂, ℬ), ℬ)

    ℬ_ν = TensND.Basis(hcat(TensND.components_canon(τ̂),
                            TensND.components_canon(ν̂),
                            TensND.components_canon(n̂)))

    B_ℰ = cod_tensor(crack, C₀; method=method, kw...)
    ribbon_ref = RibbonCrack(b, ℬ_ν)
    B_ℛ = cod_tensor(ribbon_ref, C₀; method=method, kw...)

    𝐊 = (3 * T(π)^(T(3)/2) * b / 8) * sqrt(b * n_Sy) *
        inv(B_ℛ) ⋅ B_ℰ ⋅ Σ ⋅ n̂

    Kᴵ   = 𝐊 ⋅ n̂
    Kᴵᴵ  = 𝐊 ⋅ ν̂
    Kᴵᴵᴵ = 𝐊 ⋅ τ̂
    return 𝐊, (Kᴵ, Kᴵᴵ, Kᴵᴵᴵ)
end

# =============================================================================
#  Displacement intensity factor
# =============================================================================

"""
    dif(crack, C₀, Σ; method=:auto, kw...) -> Tens{1,3}
"""
function dif(crack::MFH_Core.AbstractCrack, C₀, Σ; method::Symbol=:auto, kw...)
    B = cod_tensor(crack, C₀; method=method, kw...)
    n̂ = TensND.tensbasis(crack_basis(crack), 3)
    return B ⋅ Σ ⋅ n̂
end
