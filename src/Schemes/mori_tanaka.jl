# =============================================================================
#  mori_tanaka.jl — Mori-Tanaka scheme.
#
#  Each inclusion experiences the *average matrix strain* as far field. The
#  effective stiffness reads (for solid inclusions only)
#
#       C_MT = C₀ + (Σ_i fᵢ Nᵢ) : (f_m · I + Σ_i fᵢ A_dil^{(i)})⁻¹
#
#  with A_dil^{(i)} = (I + P_0 : (C_i − C_0))⁻¹ the dilute strain
#  concentration tensor and Nᵢ = (Cᵢ − C₀) : A_dil^{(i)}.
#
#  For cracks (`CrackDensity`) the strain concentration tensor is singular;
#  we add their density-weighted stiffness reduction
#  -C₀ : H : C₀ to the numerator (consistent with Kachanov 1992 in the
#  dilute crack limit) and skip the denominator term — physically the
#  crack volume → 0 limit anyway sends f_crack to zero.
# =============================================================================

"""
    _evaluate(rve, ::MoriTanaka, ::Val{p}; kw...) -> AbstractTens

Mori-Tanaka scheme for property `:p`
([Mori & Tanaka 1973](@cite mori1973);
[Christensen 1990](@cite christensen1990)). Dispatches on the order of
the matrix property tensor — 4th order for elasticity (`:C`), 2nd order
for conductivity (`:K`).
"""
function _evaluate(rve::RVE, ::MoriTanaka, ::Val{p}; kw...) where {p}
    P₀ = matrix_property(rve, p)
    return _mt_dispatch(rve, P₀, Val(p); kw...)
end

_mt_dispatch(rve, C₀::TensND.AbstractTens{4, 3}, ::Val{p}; kw...) where {p} =
    _mt_4(rve, C₀, Val(p); kw...)
_mt_dispatch(rve, K₀::TensND.AbstractTens{2, 3}, ::Val{p}; kw...) where {p} =
    _mt_2(rve, K₀, Val(p); kw...)

# ── 4th-order (elasticity) ──────────────────────────────────────────────────
function _mt_4(rve, C₀::TensND.AbstractTens{4, 3}, ::Val{p}; kw...) where {p}
    f_m  = matrix_volume_fraction(rve)
    Iref = _identity_like(C₀)
    A_avg = f_m * Iref          # ⟨A_dil⟩, matrix carries A_dil = I
    Nsum  = zero(C₀)
    for name in inclusion_phase_names(rve)
        a    = rve.amounts[name]
        geom = rve.phases[name].geometry
        if a isa VolumeFraction
            f = amount_value(a)
            P_i = phase_property(rve, name, p)
            A_dil = MFH_Core.strain_strain_loc(geom, P_i, C₀; kw...)
            A_avg += f * A_dil
            Nsum  += f * ((P_i - C₀) ⊡ A_dil)
        else  # CrackDensity — A_dil singular, fall back to dilute crack term
            Nsum += _phase_stiffness_contribution(rve, name, p, C₀; kw...)
        end
    end
    return C₀ + Nsum ⊡ inv(A_avg)
end

# ── 2nd-order (conductivity) ────────────────────────────────────────────────
function _mt_2(rve, K₀::TensND.AbstractTens{2, 3}, ::Val{p}; kw...) where {p}
    f_m  = matrix_volume_fraction(rve)
    Iref = _identity_like(K₀)
    A_avg = f_m * Iref
    Nsum  = zero(K₀)
    for name in inclusion_phase_names(rve)
        a    = rve.amounts[name]
        geom = rve.phases[name].geometry
        if a isa VolumeFraction
            f = amount_value(a)
            P_i = phase_property(rve, name, p)
            A_dil = MFH_Core.gradient_gradient_loc(geom, P_i, K₀; kw...)
            A_avg += f * A_dil
            Nsum  += f * ((P_i - K₀) ⋅ A_dil)
        else
            Nsum += _phase_stiffness_contribution(rve, name, p, K₀; kw...)
        end
    end
    return K₀ + Nsum ⋅ inv(A_avg)
end

# ── Identity tensor matching the algebra of the property tensor ─────────────

"""
    _identity_like(P) -> AbstractTens

Identity tensor of the same order/dimension as `P`. Used as the `A_dil`
weight of the matrix in average-strain schemes (Mori-Tanaka,
self-consistent, …).
"""
_identity_like(C::TensND.AbstractTens{4, 3}) =
    TensND.tens_Id4(Val(3), Val(eltype(C)))
_identity_like(K::TensND.AbstractTens{2, 3}) =
    TensND.tens_Id2(Val(3), Val(eltype(K)))
