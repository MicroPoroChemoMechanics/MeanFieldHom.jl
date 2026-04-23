# =============================================================================
#  contribution.jl — size-independent contribution tensors of the dilute
#  Eshelby problem (Kachanov-Sevostianov convention).
#
#  Given an inclusion of stiffness `C₁` in a matrix `C₀`, the dilute
#  effective stiffness correction is
#
#      ΔC_eff = f × N,      N = (C₁ - C₀) : A_εε,
#
#  where `A_εε` is the strain-strain localization tensor.  Dually, the
#  dilute effective compliance correction is
#
#      ΔS_eff = f × H,      H = (S₁ - S₀) : A_σσ,
#
#  with `S = C⁻¹` and `A_σσ` the stress-stress localization tensor.
#  The helpers `delta_stiffness` and `delta_compliance` apply the
#  volume fraction `f` (analogous to Budiansky density for cracks).
#
#  Conductivity analogues are implemented with 2-tensor algebra.
# =============================================================================

"""
    stiffness_contribution(incl, C₁, C₀; kw...) -> Tens{4,3}

Size-independent **stiffness contribution tensor**
`N = (C₁ - C₀) : A_εε` for an `AbstractInclusion` of stiffness `C₁`
in a matrix `C₀`.  For a dilute family of inclusions of volume fraction
`f`, the effective stiffness correction is
`ΔC_eff = f × N` — see [`delta_stiffness`](@ref).

See [Kachanov & Sevostianov (2018)](@cite kachanov2018).
"""
function stiffness_contribution(
        incl::AbstractInclusion,
        C₁::TensND.AbstractTens{4, 3},
        C₀::TensND.AbstractTens{4, 3};
        kw...
    )
    A = strain_strain_loc(incl, C₁, C₀; kw...)
    return (C₁ - C₀) ⊡ A
end

"""
    compliance_contribution(incl, C₁, C₀; kw...) -> Tens{4,3}

Size-independent **compliance contribution tensor**
`H = (S₁ - S₀) : A_σσ` for an `AbstractInclusion` of stiffness `C₁`
in a matrix `C₀` (`S = C⁻¹`).  For a dilute family, the effective
compliance correction is `ΔS_eff = f × H` — see [`delta_compliance`](@ref).

See [Kachanov & Sevostianov (2018)](@cite kachanov2018).
"""
function compliance_contribution(
        incl::AbstractInclusion,
        C₁::TensND.AbstractTens{4, 3},
        C₀::TensND.AbstractTens{4, 3};
        kw...
    )
    A_σσ = stress_stress_loc(incl, C₁, C₀; kw...)
    return (inv(C₁) - inv(C₀)) ⊡ A_σσ
end

"""
    delta_stiffness(N, f) -> Tens{4,3}

Dilute **effective stiffness correction** `ΔC = f × N` from the
size-independent contribution tensor `N` and the volume fraction `f`
of inclusions sharing that contribution.
"""
delta_stiffness(N::TensND.AbstractTens{4, 3}, f) = f * N

"""
    delta_compliance(H, f) -> Tens{4,3}

Dilute **effective compliance correction** `ΔS = f × H` from the
size-independent contribution tensor `H` and the volume fraction `f`.
(See also the crack-specific methods `delta_compliance(crack, H, ε)`
which use the Budiansky density convention and apply a geometric
prefactor.)
"""
delta_compliance(H::TensND.AbstractTens{4, 3}, f) = f * H

# =============================================================================
#  Conductivity contribution (2-tensor fields)
# =============================================================================

"""
    conductivity_contribution(incl, K₁, K₀; kw...) -> Tens{2,3}

Size-independent **conductivity contribution tensor**
`N_K = (K₁ - K₀) · A_∇∇` for an `AbstractInclusion` of conductivity
`K₁` in a matrix `K₀`.  Dilute effective correction:
`ΔK_eff = f × N_K`.
"""
function conductivity_contribution(
        incl::AbstractInclusion,
        K₁::TensND.AbstractTens{2, 3},
        K₀::TensND.AbstractTens{2, 3};
        kw...
    )
    A = gradient_gradient_loc(incl, K₁, K₀; kw...)
    return (K₁ - K₀) ⋅ A
end

"""
    resistivity_contribution(incl, K₁, K₀; kw...) -> Tens{2,3}

Size-independent **resistivity contribution tensor**
`H_R = (R₁ - R₀) · A_qq` for an `AbstractInclusion` (with `R = K⁻¹`).
Dilute effective correction: `ΔR_eff = f × H_R`.
"""
function resistivity_contribution(
        incl::AbstractInclusion,
        K₁::TensND.AbstractTens{2, 3},
        K₀::TensND.AbstractTens{2, 3};
        kw...
    )
    A_qq = flux_flux_loc(incl, K₁, K₀; kw...)
    return (inv(K₁) - inv(K₀)) ⋅ A_qq
end

"""
    delta_conductivity(N_K, f) -> Tens{2,3}

Dilute effective conductivity correction `ΔK = f × N_K`.
"""
delta_conductivity(N::TensND.AbstractTens{2, 3}, f) = f * N

"""
    delta_resistivity(H_R, f) -> Tens{2,3}

Dilute effective resistivity correction `ΔR = f × H_R`.  Generic
2-argument method; for cracks, see the 3-argument
`delta_resistivity(crack, R, ε)` with the Budiansky density prefactor.
"""
delta_resistivity(H::TensND.AbstractTens{2, 3}, f) = f * H
