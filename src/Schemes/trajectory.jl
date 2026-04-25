# =============================================================================
#  trajectory.jl — DifferentialTrajectory hierarchy and per-step
#  fraction-path resolution.
#
#  Each trajectory describes how the *target* per-phase volume fractions
#  are reached cumulatively over `nsteps` integration steps:
#
#    path_i(k) ∈ [0, 1],   path_i(0) = 0,   path_i(nsteps) = 1.
#
#  The differential ODE-step then assigns the increment Φ_k to each
#  phase by solving a small linear system (see `differential.jl`).
# =============================================================================

"""
    _resolve_paths(traj::DifferentialTrajectory, rve::RVE, nsteps::Int)
        -> Dict{Symbol, Vector{Float64}}

Build the cumulative trajectory `path_i(k) ∈ [0,1]` for each non-matrix
phase of `rve` over `nsteps` steps, normalised so `path_i(0) = 0` and
`path_i(nsteps) = 1`. The trajectory is post-validated to ensure
monotonicity and the boundary values.
"""
function _resolve_paths end

# ── Proportional ────────────────────────────────────────────────────────────

function _resolve_paths(::Proportional, rve::RVE, nsteps::Int)
    paths = Dict{Symbol, Vector{Float64}}()
    for name in inclusion_phase_names(rve)
        paths[name] = collect(range(0.0, 1.0; length = nsteps + 1))
    end
    return paths
end

# ── Sequential ──────────────────────────────────────────────────────────────

function _resolve_paths(s::Sequential, rve::RVE, nsteps::Int)
    inc_names = inclusion_phase_names(rve)
    Set(s.order) == Set(inc_names) ||
        throw(ArgumentError("Sequential.order = $(s.order) must list every " *
                            "inclusion phase exactly once; got phases $(inc_names)"))
    nphases = length(s.order)
    nphases == 0 && return Dict{Symbol, Vector{Float64}}()
    paths = Dict{Symbol, Vector{Float64}}(name => zeros(nsteps + 1) for name in s.order)

    # Each phase gets a contiguous slice of the integration steps.
    base = nsteps ÷ nphases
    rem  = nsteps - base * nphases    # distribute leftover to the first `rem` phases
    cumstep = 0
    for (i, name) in enumerate(s.order)
        n_i = base + (i ≤ rem ? 1 : 0)
        # Phase i ramps 0 → 1 between steps `cumstep` and `cumstep + n_i`
        for k in 0:nsteps
            if k ≤ cumstep
                paths[name][k + 1] = 0.0
            elseif k ≥ cumstep + n_i
                paths[name][k + 1] = 1.0
            else
                paths[name][k + 1] = (k - cumstep) / n_i
            end
        end
        cumstep += n_i
    end
    return paths
end

# ── CustomPath ──────────────────────────────────────────────────────────────

function _resolve_paths(c::CustomPath, rve::RVE, nsteps::Int)
    inc_names = inclusion_phase_names(rve)
    paths = Dict{Symbol, Vector{Float64}}()
    for name in inc_names
        haskey(c.path, name) ||
            throw(ArgumentError("CustomPath: missing trajectory for phase :$(name)"))
        traj = c.path[name]
        length(traj) == nsteps + 1 ||
            throw(ArgumentError("CustomPath[:$(name)] has length $(length(traj)); " *
                                "expected nsteps + 1 = $(nsteps + 1)"))
        traj[1] == 0 ||
            throw(ArgumentError("CustomPath[:$(name)][1] must be 0; got $(traj[1])"))
        traj[end] == 1 ||
            throw(ArgumentError("CustomPath[:$(name)][end] must be 1; got $(traj[end])"))
        # Monotone non-decreasing
        for k in 1:nsteps
            traj[k + 1] - traj[k] ≥ -1.0e-12 ||
                throw(ArgumentError("CustomPath[:$(name)] is not monotone non-decreasing"))
        end
        paths[name] = collect(Float64.(traj))
    end
    return paths
end
