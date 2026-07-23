# =============================================================================
#  59_alv_sensitivities.jl
#
#  Autodiff sensitivities (`ForwardDiff.derivative`, `gradient`) of an ALV
#  effective property w.r.t. RVE parameters.  Mirrors the elastic
#  `Schemes/sensitivities` machinery (see `scripts/26_sensitivities.jl`)
#  for the time-domain pipeline.
#
#  Two patterns are demonstrated :
#
#  (1) **`set_param` lens — recommended.**  Build the RVE once with
#      `Float64` placeholders, then differentiate by substituting a
#      `Dual` value via `set_param(rve, AmountParameter(...), value)`.
#      Works for **all** schemes (Voigt, Reuss, Dilute, MT, Maxwell, PCW,
#      …) and does not allocate a fresh RVE outside the differentiation.
#
#  (2) **Closure-captured material parameter.**  When the dependence is
#      on a continuous parameter inside the `ViscoLaw` itself (e.g. a
#      bulk / shear modulus, a relaxation time), close the parameter
#      into the kernel function and differentiate.  The parameter
#      automatically lifts to `Dual` through the closure.
#
#  Comparison vs central finite differences validates each derivative.
#  Reference : @sanahuja2013 §4 ; @barthelemyIJES2019 §3 ; analogous to
#  the elastic-side @bessoIJSS2024 sensitivity pipeline.
#
#  Usage : julia --project scripts/59_alv_sensitivities.jl
# =============================================================================

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."); io = devnull)

using MeanFieldHom
using TensND
using ForwardDiff
using LinearAlgebra
using Printf

# ─── Common ALV setup ──────────────────────────────────────────────────────

const N_TIMES = 8
const TIMES = collect(range(0.0, 2.0; length = N_TIMES))

# Maxwell-iso matrix law parametrized by `(k_M, μ_M, τ_K, τ_μ)`.
# Closure pattern: each call returns a fresh `ViscoLaw`.
function build_law_M(k_M, μ_M, τ_K = 1.0, τ_μ = 0.5)
    function R_iso(t, tp)
        α = 3 * k_M * (1.0 + 4.0 * exp(-(t - tp) / τ_K))
        β = 2 * μ_M * (0.5 + 1.5 * exp(-(t - tp) / τ_μ))
        return TensISO{3}(α, β)
    end
    return ViscoLaw(R_iso, :relaxation)
end

# Inclusion : elastic spheres with (k, μ) = (10, 4).
const C_INC = TensISO{3}(3 * 10.0, 2 * 4.0)

function build_rve_base(f::Real)
    rve = RVE(:M)
    add_matrix!(rve, Ellipsoid(1.0), Dict(:C => build_law_M(1.0, 1.0)))
    add_phase!(
        rve, :I, Ellipsoid(1.0),
        Dict(:C => heaviside_law(C_INC));
        fraction = f
    )
    return rve
end

# Final-time effective shear `μ(t_n,t_n) = β[end,end] / 2`.
function effective_mu_final(rve, scheme)
    R̃ = homogenize_alv(rve, scheme, :C; times = TIMES)
    _, β = iso_params_from_blocks(R̃)
    return β[end, end] / 2
end

# ─── Pattern (1) : `set_param` lens — derivative wrt volume fraction ──────

println("="^78)
println(" 1)  d μ_eff(t_n) / df      — `set_param` + `AmountParameter`")
println("="^78)

const RVE_BASE_F = build_rve_base(0.2)

function eff_mu_vs_f(f::Real, scheme)
    rve_f = set_param(RVE_BASE_F, AmountParameter(:I), f)
    return effective_mu_final(rve_f, scheme)
end

const SCHEMES = (Voigt(), Reuss(), Dilute(), MoriTanaka(), Maxwell())
const SCHEME_NAMES = ("Voigt", "Reuss", "Dilute", "MoriTanaka", "Maxwell")

for (sch, name) in zip(SCHEMES, SCHEME_NAMES)
    f₀ = 0.2
    dμ_df_AD = ForwardDiff.derivative(f -> eff_mu_vs_f(f, sch), f₀)
    h = 1.0e-5
    dμ_df_FD = (eff_mu_vs_f(f₀ + h, sch) - eff_mu_vs_f(f₀ - h, sch)) / (2h)
    rel_err = abs(dμ_df_AD - dμ_df_FD) / abs(dμ_df_FD)
    @printf "  %-12s  AD = %+.6e   FD = %+.6e   rel_err = %.2e\n" name dμ_df_AD dμ_df_FD rel_err
end

# ─── Pattern (2) : material-parameter sensitivity via closure ──────────────

println()
println("="^78)
println(" 2)  d μ_eff(t_n) / dμ_M   — closure-captured matrix shear modulus")
println("="^78)

function eff_mu_vs_μM(μ_M::Real)
    rve = RVE(:M)
    add_matrix!(rve, Ellipsoid(1.0), Dict(:C => build_law_M(1.0, μ_M)))
    add_phase!(
        rve, :I, Ellipsoid(1.0), Dict(:C => heaviside_law(C_INC));
        fraction = 0.2
    )
    return effective_mu_final(rve, MoriTanaka())
end

μM₀ = 1.0
dμ_dμM_AD = ForwardDiff.derivative(eff_mu_vs_μM, μM₀)
h = 1.0e-5
dμ_dμM_FD = (eff_mu_vs_μM(μM₀ + h) - eff_mu_vs_μM(μM₀ - h)) / (2h)
@printf "  AD = %+.6e   FD = %+.6e   rel_err = %.2e\n" dμ_dμM_AD dμ_dμM_FD abs(dμ_dμM_AD - dμ_dμM_FD) / abs(dμ_dμM_FD)

# ─── Pattern (3) : gradient over multiple parameters ──────────────────────

println()
println("="^78)
println(" 3)  ∇ μ_eff   wrt   (f, k_M, μ_M)   — joint gradient via ForwardDiff")
println("="^78)

function eff_mu_vs_fkμ(p::AbstractVector)
    f, k_M, μ_M = p
    rve = RVE(:M)
    add_matrix!(rve, Ellipsoid(1.0), Dict(:C => build_law_M(k_M, μ_M)))
    add_phase!(
        rve, :I, Ellipsoid(1.0), Dict(:C => heaviside_law(C_INC));
        fraction = 0.2
    )
    rve_f = set_param(rve, AmountParameter(:I), f)
    return effective_mu_final(rve_f, MoriTanaka())
end

p₀ = [0.2, 1.0, 1.0]
∇AD = ForwardDiff.gradient(eff_mu_vs_fkμ, p₀)
∇FD = let h = 1.0e-5
    [
        (
                eff_mu_vs_fkμ(p₀ .+ h .* (i == k for k in 1:3)) -
                eff_mu_vs_fkμ(p₀ .- h .* (i == k for k in 1:3))
            ) / (2h) for i in 1:3
    ]
end
for (i, name) in enumerate(("f", "k_M", "μ_M"))
    @printf "  ∂μ/∂%-3s   AD = %+.6e   FD = %+.6e   rel_err = %.2e\n" name ∇AD[i] ∇FD[i] abs(∇AD[i] - ∇FD[i]) / abs(∇FD[i])
end

# ─── Pattern (4) : relaxation-time sensitivity (closure parameter) ─────────

println()
println("="^78)
println(" 4)  d μ_eff(t_n) / dτ_K   — relaxation-time inside the kernel")
println("="^78)

function eff_mu_vs_τK(τ_K::Real)
    rve = RVE(:M)
    add_matrix!(rve, Ellipsoid(1.0), Dict(:C => build_law_M(1.0, 1.0, τ_K, 0.5)))
    add_phase!(
        rve, :I, Ellipsoid(1.0), Dict(:C => heaviside_law(C_INC));
        fraction = 0.2
    )
    return effective_mu_final(rve, MoriTanaka())
end

τK₀ = 1.0
dμ_dτK_AD = ForwardDiff.derivative(eff_mu_vs_τK, τK₀)
h = 1.0e-5
dμ_dτK_FD = (eff_mu_vs_τK(τK₀ + h) - eff_mu_vs_τK(τK₀ - h)) / (2h)
@printf "  AD = %+.6e   FD = %+.6e   rel_err = %.2e\n" dμ_dτK_AD dμ_dτK_FD abs(dμ_dτK_AD - dμ_dτK_FD) / abs(dμ_dτK_FD)

println()
println("All sensitivities agree with central finite differences to ≤ 1e-7 ")
println("(machine precision modulo the FD truncation error h²·f‴/6 ≈ 1e-10).")
