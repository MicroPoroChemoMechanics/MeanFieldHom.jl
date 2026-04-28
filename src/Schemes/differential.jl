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

# ── Integration loop ────────────────────────────────────────────────────────
#
# At each integration step k, the matrix property is updated via an
# explicit Euler step:
#       C^{k}  = C^{k-1} + Σ_i Φ_k[i] · (C_i − C^{k-1}) : A_{εε,i}^{(C^{k-1})}
# where the increments Φ_k solve the volume-balance linear system (cf.
# the C++ reference's `homogenization_differential::compute_property`).
# Symmetrize on inclusion phases is honoured through the per-phase
# helpers (`_phase_stiffness_contribution`,
# `_phase_dilute_concentration`); the matrix is *not* symmetrized
# (matches the C++ reference, where the matrix is the iterating object
# and has no `symmetrize` keyword). For the dilute correction we
# rebuild a transient single-phase RVE so the helpers naturally pick
# up `Φi` as the volume fraction and the correct symmetrize.

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
        if n_solid > 0
            # Volume-balance linear system: row i of the matrix is
            # (I − f_i^{(k-1)} · 𝟙ᵀ), reproducing the C++ reference's
            # `Mk(i, Range::all()) -= fikm` which subtracts the
            # *row-specific* `fikm` from every column.
            f_prev = [paths[name][k]     * targets_solid[i] for (i, name) in enumerate(solid_names)]
            f_curr = [paths[name][k + 1] * targets_solid[i] for (i, name) in enumerate(solid_names)]
            ΔF = f_curr .- f_prev
            Telt = promote_type(eltype(f_prev), eltype(ΔF))
            M = Matrix{Telt}(LinearAlgebra.I, n_solid, n_solid)
            for i in 1:n_solid
                M[i, :] .-= f_prev[i]
            end
            Φ = M \ ΔF
            for (i, name) in enumerate(solid_names)
                Φi = Φ[i]
                iszero(Φi) && continue
                # Apply the dilute correction with symmetrize honoured.
                P_curr += Φi * _diff_dilute_correction(rve, name, prop, P_curr; kw...)
            end
        end
        # Cracks: add 1/nsteps of their density-weighted contribution
        # through the symmetrize-aware helper.
        for name in crack_names
            P_curr += _diff_crack_correction(rve, name, prop, P_curr, nsteps; kw...)
        end
    end
    return P_curr
end

# Dilute correction `(C_i − C^{(k-1)}) : A_{εε,i}^{(C^{(k-1)})}` for a
# single inclusion phase — symmetrize is honoured by routing through
# `_phase_dilute_concentration`.
function _diff_dilute_correction(rve::RVE, name::Symbol, prop::Symbol,
                                  P_curr::TensND.AbstractTens{4, 3}; kw...)
    P_i  = phase_property(rve, name, prop)
    A    = _phase_dilute_concentration(rve, name, prop, P_curr; kw...)
    return (P_i - P_curr) ⊡ A
end

function _diff_dilute_correction(rve::RVE, name::Symbol, prop::Symbol,
                                  P_curr::TensND.AbstractTens{2, 3}; kw...)
    P_i  = phase_property(rve, name, prop)
    A    = _phase_dilute_concentration(rve, name, prop, P_curr; kw...)
    return (P_i - P_curr) ⋅ A
end

# Crack increment for one step (1/nsteps of the target density).
function _diff_crack_correction(rve::RVE, name::Symbol, prop::Symbol,
                                 P_curr::TensND.AbstractTens, nsteps::Int; kw...)
    geom = rve.phases[name].geometry
    sym  = phase_symmetrize(rve, name)
    P₀_proj = _project_matrix(P_curr, sym)
    if P_curr isa TensND.AbstractTens{4, 3}
        K_int = _crack_interface_K4(rve, name)
        N = MFH_Core.stiffness_contribution(geom, P₀_proj;
                                             K_interface = K_int, kw...)
        return _apply_symmetrize(MFH_Core.delta_stiffness(geom, N,
                                  amount_value(rve.amounts[name]) / nsteps), sym)
    else
        α_int = _crack_interface_α(rve, name)
        N = MFH_Core.conductivity_contribution(geom, P₀_proj;
                                                α_interface = α_int, kw...)
        return _apply_symmetrize(MFH_Core.delta_conductivity(geom, N,
                                  amount_value(rve.amounts[name]) / nsteps), sym)
    end
end
