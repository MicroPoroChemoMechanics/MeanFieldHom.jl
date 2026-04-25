# =============================================================================
#  voigt.jl — Voigt (uniform-strain) upper bound.
#
#       C_Voigt = Σ_i f_i  C_i
#
#  Cracks (`CrackDensity`) are ignored: their volume contribution → 0
#  in the penny limit. The matrix volume fraction is the implicit
#  complement `1 - Σ f_inc`.
# =============================================================================

"""
    _evaluate(rve, ::Voigt, ::Val{p}; kw...) -> AbstractTens

Voigt upper bound on the effective property `:p`:
``\\langle\\mathbb C\\rangle = \\sum_i f_i \\mathbb C_i``.

Phases carrying a [`CrackDensity`](@ref) instead of a
[`VolumeFraction`](@ref) are ignored (their volume contribution is
zero); use a Hill-tensor-aware scheme (e.g. [`Dilute`](@ref) or
[`MoriTanaka`](@ref)) to capture crack effects.

Reference: [Hill (1965)](@cite hill1965).
"""
function _evaluate(rve::RVE, ::Voigt, ::Val{p}; kw...) where {p}
    f_m  = matrix_volume_fraction(rve)
    Ceff = f_m * matrix_property(rve, p)
    for name in inclusion_phase_names(rve)
        a = rve.amounts[name]
        if a isa VolumeFraction
            Ceff += amount_value(a) * phase_property(rve, name, p)
        end
    end
    return Ceff
end
