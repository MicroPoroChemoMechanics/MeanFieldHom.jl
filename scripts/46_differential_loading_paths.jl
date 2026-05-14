# =============================================================================
#  46_differential_loading_paths.jl
#
#  Demonstrates the **path-dependence of the differential homogenisation
#  scheme** (DEM) by computing the effective stiffness of a 3-phase
#  composite (matrix + 2 solid inclusions) along several incorporation
#  trajectories `τ -> f_α(τ)` :
#
#    1. `Proportional()`           — both phases grow linearly together.
#    2. `Sequential([:I1, :I2])`   — phase 1 first, then phase 2.
#    3. `Sequential([:I2, :I1])`   — phase 2 first, then phase 1.
#    4. `Path(:I1 => τ -> τ², :I2 => τ -> 2τ - τ²)`
#                                  — phase 2 frontloaded (concave), phase
#                                    1 backloaded (convex), both reach 1
#                                    at τ = 1.
#
#  All trajectories satisfy `f_α(0) = 0`, `f_α(1) = 1`, so they reach
#  the same **target** volume fractions at τ = 1.  But the resulting
#  effective stiffness `C^hom(τ = 1)` **differs** : DEM is genuinely
#  path-dependent because each infinitesimal volume increment is added
#  as a dilute inclusion in the *current* effective medium, which
#  itself depends on the prior incorporation history.
#
#  Output : `scripts/figures/46_differential_loading_paths.png` plotting
#  the bulk and shear moduli of `C^hom(τ)` along each path.
#
#  Usage  : julia --project scripts/46_differential_loading_paths.jl
# =============================================================================

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."); io = devnull)

using MeanFieldHom
using TensND
using OrdinaryDiffEq      # for `solve` access if user wants intermediate states
using Plots
using Printf

# ─── RVE : matrix + 2 solid inclusions ─────────────────────────────────────

const C_M = TensISO{3}(3 * 5.0, 2 * 2.0)        # matrix : k = 5,    μ = 2
const C_I1 = TensISO{3}(3 * 30.0, 2 * 12.0)      # phase 1 : stiff   (5× stiffer)
const C_I2 = TensISO{3}(3 * 0.5, 2 * 0.2)       # phase 2 : compliant (0.1× softer)

const F1, F2 = 0.2, 0.2    # target volume fractions

function build_rve()
    rve = RVE(:M)
    add_matrix!(rve, Ellipsoid(1.0), Dict(:C => C_M))
    add_phase!(rve, :I1, Ellipsoid(1.0), Dict(:C => C_I1); fraction = F1)
    add_phase!(rve, :I2, Ellipsoid(1.0), Dict(:C => C_I2); fraction = F2)
    return rve
end

# ─── Compute C^hom along τ for each trajectory ─────────────────────────────
#
# To get the full curve `C^hom(τ)` we ask the scheme for a high
# `nsteps` (saveat density along τ) and read the saved trajectory.

const NSTEPS = 200
const TAU = collect(range(0.0, 1.0; length = NSTEPS + 1))

function eval_path(traj)
    rve = build_rve()
    # The current public homogenize() returns C^hom(τ = 1) only.  To
    # recover the saved trajectory `C^hom(τ)` we re-build the ODE
    # integration here using the same RHS — keeps the demo
    # self-contained and explicit.
    paths = MeanFieldHom.Schemes._resolve_paths(traj, rve, NSTEPS)
    P_init = matrix_property(rve, :C)
    sym_tag = MeanFieldHom.Schemes._symmetry_tag(P_init)
    x0 = collect(TensND.get_data(P_init))
    ode_p = (
        rve = rve, prop = :C, paths = paths,
        solid_names = [:I1, :I2], crack_names = Symbol[],
        targets = Dict(:I1 => F1, :I2 => F2),
        sym_tag = sym_tag, proto = P_init, kw = NamedTuple(),
    )
    rhs! = (du, u, p, τ) -> MeanFieldHom.Schemes._diff_ode_rhs!(du, u, p, τ)
    prob = ODEProblem(rhs!, x0, (0.0, 1.0), ode_p)
    sol = solve(
        prob, Tsit5(); abstol = 1.0e-9, reltol = 1.0e-7,
        saveat = TAU, dense = false
    )
    # Extract bulk and shear moduli of every saved state.
    α = [u[1] for u in sol.u]   # = 3 k_eff
    β = [u[2] for u in sol.u]   # = 2 μ_eff
    return α ./ 3, β ./ 2
end

println("Computing four loading-path scenarios on the same target (f₁=$F1, f₂=$F2)…")

paths_to_run = (
    ("Proportional", Proportional()),
    ("Sequential :I1 → :I2", Sequential([:I1, :I2])),
    ("Sequential :I2 → :I1", Sequential([:I2, :I1])),
    (
        "Path (I1 ∝ τ²,  I2 ∝ 2τ−τ²)", Path(
            Dict(
                :I1 => τ -> τ^2,
                :I2 => τ -> 2τ - τ^2,
            )
        ),
    ),
)

results = Dict{String, Tuple{Vector{Float64}, Vector{Float64}}}()
for (name, traj) in paths_to_run
    println("  $name…")
    k, μ = eval_path(traj)
    results[name] = (k, μ)
    @printf "    k_eff(τ=1) = %.5f   μ_eff(τ=1) = %.5f\n" k[end] μ[end]
end

# ─── Plot ──────────────────────────────────────────────────────────────────

p_k = plot(
    xlabel = "τ  (fictitious incorporation time)",
    ylabel = "k_eff(τ)",
    title = "Bulk modulus along the trajectory",
    legend = :right
)
p_μ = plot(
    xlabel = "τ",
    ylabel = "μ_eff(τ)",
    title = "Shear modulus along the trajectory",
    legend = :right
)

colors = (:black, :red, :blue, :green)
for (i, (name, _)) in enumerate(paths_to_run)
    k, μ = results[name]
    plot!(p_k, TAU, k; label = name, color = colors[i], linewidth = 2)
    plot!(p_μ, TAU, μ; label = name, color = colors[i], linewidth = 2)
end

fig = plot(p_k, p_μ; layout = (1, 2), size = (1400, 600))
mkpath(joinpath(@__DIR__, "figures"))
out = joinpath(@__DIR__, "figures", "46_differential_loading_paths.png")
savefig(fig, out)
println("\nSaved : $out")

# ─── Numeric report ────────────────────────────────────────────────────────

println()
println("═══════════════════════════════════════════════════════════════════")
println(" Effective moduli at τ = 1 (same target volume fractions, different paths)")
println("═══════════════════════════════════════════════════════════════════")
@printf "  %-32s  %-12s  %-12s\n" "trajectory" "k_eff" "μ_eff"
for (name, _) in paths_to_run
    k, μ = results[name]
    @printf "  %-32s  %-12.5f  %-12.5f\n" name k[end] μ[end]
end
println()
println("The DEM scheme is genuinely **path-dependent** : different")
println("incorporation sequences `τ → f_α(τ)` reaching the same target")
println("volume fractions `f_α^∞` at τ=1 give different `C^hom(τ=1)`.")
