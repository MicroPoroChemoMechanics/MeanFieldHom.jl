# =============================================================================
#  cod_analytical.jl — closed-form COD tensors.
#  Moduli extractors live in Core.moduli.jl; this file only contains the
#  COD formulas.
# =============================================================================

"""
    _elliptic_CS(η) -> (𝒞, 𝒮, ℰ)

Angular integrals (paper eq. 1532).
"""
@inline function _elliptic_CS(η::T) where {T <: Number}
    η² = η^2
    k² = one(T) - η²
    if T <: Real && iszero(k²)
        q = T(π) / T(4)
        return q, q, T(π) / T(2)
    else
        𝒦 = GenericElliptic.ell_K(k²)
        ℰ = GenericElliptic.ell_E(k²)
        𝒞 = (ℰ - η² * 𝒦) / k²
        𝒮 = (𝒦 - ℰ) / k²
        return 𝒞, 𝒮, ℰ
    end
end

# =============================================================================
#  ISO matrix
# =============================================================================

"""
    _cod_iso_ellipse(c::EllipticCrack, E, ν) -> Tens{2,3}
Paper eq. (1514).
"""
function _cod_iso_ellipse(c::EllipticCrack{T}, E, ν) where {T <: Number}
    η = aspect_ratio(c)
    𝒞, 𝒮, ℰ = _elliptic_CS(η)
    χ = 8 * (one(T) - ν^2) / (3 * E)
    Bll = χ / ((one(T) - ν) * 𝒞 + η^2 * 𝒮)
    Bmm = χ / ((one(T) - ν) * η^2 * 𝒮 + 𝒞)
    Bnn = χ / ℰ
    return TensND.Tens(Diagonal([Bll, Bmm, Bnn]), crack_basis(c))
end

"""
    _cod_iso_ribbon(c::RibbonCrack, E, ν) -> Tens{2,3}
Paper eq. (1576).
"""
function _cod_iso_ribbon(c::RibbonCrack{T}, E, ν) where {T <: Number}
    χ = T(π) * (one(T) - ν^2) / E
    Bll = χ / (one(T) - ν)
    return TensND.Tens(Diagonal([Bll, χ, χ]), crack_basis(c))
end

# =============================================================================
#  TI matrix (axis ≡ n̂)
# =============================================================================

@inline function _ti_sigma_gamma(E::T, H::T, ν₁::T, ν₂::T, Γ::T) where {T <: Number}
    return sqrt(T(2)) * sqrt(
        (one(T) - Γ * ν₂) / (Γ * (one(T) - ν₁)) +
            sqrt((one(T) - H * ν₂^2) / (H * (one(T) - ν₁^2))),
    )
end

"""
    _cod_ti_ellipse(c, E, H, ν₁, ν₂, Γ) -> Tens{2,3}
Paper eqs. 1596-1634.
"""
function _cod_ti_ellipse(c::EllipticCrack{T}, E::T, H::T, ν₁::T, ν₂::T, Γ::T) where {T <: Number}
    η = aspect_ratio(c)
    σᵞ = _ti_sigma_gamma(E, H, ν₁, ν₂, Γ)
    𝒞, 𝒮, ℰ = _elliptic_CS(η)
    χ = 4 * σᵞ * (one(T) - ν₁^2) / (3 * E)
    β = σᵞ / 2 * sqrt(Γ) * (one(T) - ν₁)
    Bll = χ / (β * 𝒞 + η^2 * 𝒮)
    Bmm = χ / (β * η^2 * 𝒮 + 𝒞)
    Bnn = χ * sqrt((one(T) - H * ν₂^2) / (H * (one(T) - ν₁^2))) / ℰ
    return TensND.Tens(Diagonal([Bll, Bmm, Bnn]), crack_basis(c))
end

# Penny specialisation (η = 1)
function _cod_ti_ellipse(c::EllipticCrack{T, Penny}, E::T, H::T, ν₁::T, ν₂::T, Γ::T) where {T <: Number}
    σᵞ = _ti_sigma_gamma(E, H, ν₁, ν₂, Γ)
    χ = T(π) * (one(T) - ν₁^2) / E
    Bnn = χ * sqrt((one(T) - H * ν₂^2) / (H * (one(T) - ν₁^2)))
    Bmm = χ * σᵞ / 2
    Bll = χ / (sqrt(Γ) * (one(T) - ν₁))
    return TensND.Tens(Diagonal([Bll, Bmm, Bnn]), crack_basis(c))
end

"""
    _cod_ti_ribbon(c, E, H, ν₁, ν₂, Γ) -> Tens{2,3}
Paper eqs. 1669-1681.
"""
function _cod_ti_ribbon(c::RibbonCrack{T}, E::T, H::T, ν₁::T, ν₂::T, Γ::T) where {T <: Number}
    σᵞ = _ti_sigma_gamma(E, H, ν₁, ν₂, Γ)
    χ = T(π) * (one(T) - ν₁^2) / E
    Bnn = χ * sqrt((one(T) - H * ν₂^2) / (H * (one(T) - ν₁^2))) * σᵞ / 2
    Bmm = χ * σᵞ / 2
    Bll = χ / (sqrt(Γ) * (one(T) - ν₁))
    return TensND.Tens(Diagonal([Bll, Bmm, Bnn]), crack_basis(c))
end
