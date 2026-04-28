# =============================================================================
#  contribution_helpers.jl — phase-level dilute contribution helpers shared
#  by Dilute / DiluteDual / Mori-Tanaka / Maxwell / PCW / SC / Diff schemes.
#
#  Each phase carries either a `VolumeFraction` (solid inclusion) or a
#  `CrackDensity` (flat crack); the contribution tensors and their
#  geometry-specific prefactors differ accordingly. The helpers below pick
#  the right combination at dispatch time.
# =============================================================================

# ── Stiffness / conductivity contribution ────────────────────────────────────
#
#  Solid inclusion (VolumeFraction)
#    elasticity   : Δ_solid = f · stiffness_contribution(geom, C_i, C₀)
#    conductivity : Δ_solid = f · conductivity_contribution(geom, K_i, K₀)
#
#  Flat crack (CrackDensity)
#    elasticity   : Δ_crack = (4π/3 or π) · ε · stiffness_contribution(crack, C₀)
#    conductivity : Δ_crack = (4π/3 or π) · ε · conductivity_contribution(crack, K₀)
#
#  Dispatch on the tensor order (`AbstractTens{4,3}` vs `AbstractTens{2,3}`)
#  picks the right `stiffness_contribution` / `conductivity_contribution`.

"""
    _phase_stiffness_contribution(rve, name, prop::Symbol, P₀; kw...)

Aggregate contribution of phase `name` to the dilute *stiffness*
correction at reference `P₀` for property `prop` (`:C` for elasticity,
`:K` for conductivity).
"""
function _phase_stiffness_contribution(
        rve::RVE, name::Symbol, prop::Symbol,
        P₀::TensND.AbstractTens{4, 3}; kw...
    )
    a = rve.amounts[name]
    geom = rve.phases[name].geometry
    sym = phase_symmetrize(rve, name)
    P₀_proj = _project_matrix(P₀, sym)
    if a isa VolumeFraction
        P_i = phase_property(rve, name, prop)
        N   = MFH_Core.stiffness_contribution(geom, P_i, P₀_proj; kw...)
        return amount_value(a) * _apply_symmetrize(N, sym)
    else  # CrackDensity
        K_int = _crack_interface_K4(rve, name)
        N = MFH_Core.stiffness_contribution(geom, P₀_proj;
                                             K_interface = K_int, kw...)
        return _apply_symmetrize(MFH_Core.delta_stiffness(geom, N, amount_value(a)), sym)
    end
end

function _phase_stiffness_contribution(
        rve::RVE, name::Symbol, prop::Symbol,
        P₀::TensND.AbstractTens{2, 3}; kw...
    )
    a = rve.amounts[name]
    geom = rve.phases[name].geometry
    sym = phase_symmetrize(rve, name)
    if a isa VolumeFraction
        P_i = phase_property(rve, name, prop)
        N   = MFH_Core.conductivity_contribution(geom, P_i, P₀; kw...)
        return amount_value(a) * _apply_symmetrize(N, sym)
    else
        α_int = _crack_interface_α(rve, name)
        N = MFH_Core.conductivity_contribution(geom, P₀;
                                                α_interface = α_int, kw...)
        return _apply_symmetrize(MFH_Core.delta_conductivity(geom, N, amount_value(a)), sym)
    end
end

# ── Helpers : pull optional interface-stiffness properties from the RVE ────
#
# A crack phase can carry an optional spring-like interface stiffness via
# either of two property keys :
#   * `:K_interface`  for elasticity   (a `Tens{2,3}` 3×3 symmetric)
#   * `:α_interface`  for conductivity (a `Real` scalar conductance)
# Returns `nothing` when the property is absent — the existing
# traction-free / free-flux pipeline is then used unchanged.

function _crack_interface_K4(rve::RVE, name::Symbol)
    props = rve.phases[name].properties
    return haskey(props, :K_interface) ? props[:K_interface] : nothing
end

function _crack_interface_α(rve::RVE, name::Symbol)
    props = rve.phases[name].properties
    return haskey(props, :α_interface) ? props[:α_interface] : nothing
end

# ── Compliance / resistivity contribution ────────────────────────────────────

"""
    _phase_compliance_contribution(rve, name, prop::Symbol, P₀; kw...)

Aggregate contribution of phase `name` to the dilute *compliance*
correction (resistivity for `:K`).
"""
function _phase_compliance_contribution(
        rve::RVE, name::Symbol, prop::Symbol,
        P₀::TensND.AbstractTens{4, 3}; kw...
    )
    a = rve.amounts[name]
    geom = rve.phases[name].geometry
    sym = phase_symmetrize(rve, name)
    if a isa VolumeFraction
        P_i = phase_property(rve, name, prop)
        H   = compliance_contribution(geom, P_i, P₀; kw...)
        return amount_value(a) * _apply_symmetrize(H, sym)
    else
        K_int = _crack_interface_K4(rve, name)
        H = compliance_contribution(geom, P₀; K_interface = K_int, kw...)
        return _apply_symmetrize(delta_compliance(geom, H, amount_value(a)), sym)
    end
end

function _phase_compliance_contribution(
        rve::RVE, name::Symbol, prop::Symbol,
        P₀::TensND.AbstractTens{2, 3}; kw...
    )
    a = rve.amounts[name]
    geom = rve.phases[name].geometry
    sym = phase_symmetrize(rve, name)
    if a isa VolumeFraction
        P_i = phase_property(rve, name, prop)
        R   = MFH_Core.resistivity_contribution(geom, P_i, P₀; kw...)
        return amount_value(a) * _apply_symmetrize(R, sym)
    else
        α_int = _crack_interface_α(rve, name)
        R = compliance_contribution(geom, P₀; α_interface = α_int, kw...)
        return _apply_symmetrize(delta_resistivity(geom, R, amount_value(a)), sym)
    end
end

# ── Strain-strain localization (Mori-Tanaka, SC, …) ──────────────────────────

"""
    _phase_dilute_concentration(rve, name, prop::Symbol, P₀; kw...) -> AbstractTens

Strain-strain (or gradient-gradient) dilute concentration tensor
``\\mathbb A_{\\varepsilon\\varepsilon}^{(i)}`` for phase `name` in the
reference medium `P₀`. Used by Mori-Tanaka and self-consistent kernels
when a per-phase tensor ``\\mathbb A_i`` is required (rather than just
the contribution sum).

For cracks, the strain concentration tensor is singular (the crack is
infinitely compliant in its normal direction); the helper returns the
zero tensor as a placeholder — schemes that need crack handling should
short-circuit on `geom isa AbstractCrack`.
"""
function _phase_dilute_concentration(
        rve::RVE, name::Symbol, prop::Symbol,
        P₀::TensND.AbstractTens{4, 3}; kw...
    )
    geom = rve.phases[name].geometry
    geom isa MFH_Core.AbstractCrack && return zero(P₀)   # caller must handle cracks separately
    P_i = phase_property(rve, name, prop)
    sym = phase_symmetrize(rve, name)
    P₀_proj = _project_matrix(P₀, sym)
    A = MFH_Core.strain_strain_loc(geom, P_i, P₀_proj; kw...)
    return _apply_symmetrize(A, sym)
end

function _phase_dilute_concentration(
        rve::RVE, name::Symbol, prop::Symbol,
        P₀::TensND.AbstractTens{2, 3}; kw...
    )
    geom = rve.phases[name].geometry
    geom isa MFH_Core.AbstractCrack && return zero(P₀)
    P_i = phase_property(rve, name, prop)
    sym = phase_symmetrize(rve, name)
    P₀_proj = _project_matrix(P₀, sym)
    A = MFH_Core.gradient_gradient_loc(geom, P_i, P₀_proj; kw...)
    return _apply_symmetrize(A, sym)
end
