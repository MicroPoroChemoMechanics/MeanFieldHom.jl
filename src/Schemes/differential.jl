# =============================================================================
#  differential.jl — Differential homogenisation scheme.
#
#  Integrates the Norris ODE
#       dC / df_i = (C_i - C) : A_dil^{(i)}(C)
#  step by step along a user-selectable trajectory, replicating the
#  Echoes formulation: at step k the per-phase increment Φ_k is the
#  solution of the linear system
#       (I - diag(path_j(k-1) · f_j_target)) · Φ_k = (path_i(k) - path_i(k-1)) · f_i_target
#  and the matrix is updated by an explicit Euler step using the
#  dilute scheme on each phase increment Φ_k.
#
#  Cracks: their density does not enter the linear system (no volume
#  competition); each step adds the density-weighted dilute crack
#  contribution to the running estimate.
# =============================================================================

"""
    _evaluate(rve, scheme::DifferentialScheme, ::Val{p}; kw...) -> AbstractTens

Differential homogenisation scheme for property `:p`
([Norris 1985](@cite norris1985)). Integrates the differential equation

```math
\\frac{\\mathrm d \\mathbb C}{\\mathrm d f_i} = (\\mathbb C_i - \\mathbb C):\\mathbb A_\\mathrm{dil}^{(i)}(\\mathbb C)
```

along the trajectory specified by `scheme.trajectory` over
`scheme.options.nsteps` steps. Cracks (`CrackDensity`) are integrated
separately: at each step a fraction `1/nsteps` of the target density is
applied through the dilute crack contribution evaluated at the current
matrix.
"""
function _evaluate(rve::RVE, scheme::DifferentialScheme, ::Val{p}; kw...) where {p}
    nsteps = get(scheme.options, :nsteps, 100)
    paths  = _resolve_paths(scheme.trajectory, rve, nsteps)
    P_init = matrix_property(rve, p)
    return _diff_integrate(rve, paths, nsteps, p, P_init; kw...)
end

# Order-aware contribution wrappers — pick stiffness_contribution / delta_stiffness
# for 4th-order (elasticity) and conductivity_contribution / delta_conductivity
# for 2nd-order (transport).

_solid_contrib(geom, P_i::TensND.AbstractTens{4, 3}, P₀::TensND.AbstractTens{4, 3}; kw...) =
    MFH_Core.stiffness_contribution(geom, P_i, P₀; kw...)
_solid_contrib(geom, P_i::TensND.AbstractTens{2, 3}, P₀::TensND.AbstractTens{2, 3}; kw...) =
    MFH_Core.conductivity_contribution(geom, P_i, P₀; kw...)

_crack_contrib(geom, P₀::TensND.AbstractTens{4, 3}; kw...) =
    MFH_Core.stiffness_contribution(geom, P₀; kw...)
_crack_contrib(geom, P₀::TensND.AbstractTens{2, 3}; kw...) =
    MFH_Core.conductivity_contribution(geom, P₀; kw...)

_apply_crack_density(geom, N, ε, ::TensND.AbstractTens{4, 3}) =
    MFH_Core.delta_stiffness(geom, N, ε)
_apply_crack_density(geom, N, ε, ::TensND.AbstractTens{2, 3}) =
    MFH_Core.delta_conductivity(geom, N, ε)

# ── Integration loop ────────────────────────────────────────────────────────

function _diff_integrate(rve::RVE{T},
                         paths::AbstractDict{Symbol},
                         nsteps::Int, prop::Symbol, P_init; kw...) where {T}
    P_curr = P_init
    # Gather solid (VolumeFraction) and crack (CrackDensity) phase names.
    solid_names = Symbol[]
    crack_names = Symbol[]
    targets_solid = T[]
    for name in inclusion_phase_names(rve)
        a = rve.amounts[name]
        if a isa VolumeFraction
            push!(solid_names, name)
            push!(targets_solid, amount_value(a))
        else
            push!(crack_names, name)
        end
    end

    n_solid = length(solid_names)
    for k in 1:nsteps
        # Solid increments via the linear system Φ_k = (I - diag(...))^{-1} ΔF_k
        if n_solid > 0
            f_prev = [paths[name][k]     * targets_solid[i] for (i, name) in enumerate(solid_names)]
            f_curr = [paths[name][k + 1] * targets_solid[i] for (i, name) in enumerate(solid_names)]
            ΔF = f_curr .- f_prev
            M  = LinearAlgebra.I(n_solid) - LinearAlgebra.Diagonal(f_prev)
            Φ  = M \ ΔF
            for (i, name) in enumerate(solid_names)
                Φi = Φ[i]
                Φi == 0 && continue
                geom = rve.phases[name].geometry
                P_i  = phase_property(rve, name, prop)
                P_curr += Φi * _solid_contrib(geom, P_i, P_curr; kw...)
            end
        end
        # Cracks: add 1/nsteps of their density-weighted contribution
        for name in crack_names
            geom = rve.phases[name].geometry
            ε    = amount_value(rve.amounts[name])
            N    = _crack_contrib(geom, P_curr; kw...)
            P_curr += _apply_crack_density(geom, N, ε / nsteps, P_curr)
        end
    end
    return P_curr
end
