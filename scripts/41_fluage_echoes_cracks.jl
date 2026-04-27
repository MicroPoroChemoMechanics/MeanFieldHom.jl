# =============================================================================
#  41_fluage_echoes_cracks.jl
#
#  Julia reproduction of
#  `tests/python/creep/fluage_echoes_cracks.py` (in its **pure-crack
#  limit** — no interface stiffness).
#
#  The Python benchmark uses cracks **with interface stiffness**
#  (`Rn(t,t')`, `Rt(t,t')`); MeanFieldHom v0.6.0 implements the
#  **pure penny crack** in an iso ALV matrix (traction-free crack
#  surface, no interface compliance).  This script demonstrates the
#  available capability via the Dilute scheme on the compliance side:
#
#       J̃_eff(t,t') = J̃_M(t,t') + (4π/3) · ε³ᵈ · H̃(t,t')
#
#  with `H̃` from [`compliance_contribution_alv`](@ref) and the
#  Budiansky-O'Connell density factor from [`delta_compliance_alv`].
#
#  Setup (matrix only — no interface) :
#    * iso ALV matrix R(t,t') = C∞ + (C₀(1+0.2 √t') − C∞) exp(-(t-t')/τ)
#      with `(k₀, μ₀) = (5, 2)`, `(k_∞, μ_∞) = (3, 1)`, `τ = 1`.
#    * penny cracks (η = 1, normal e_3), density d = 0.7.
#
#  Output : effective uniaxial relaxation modulus `Eᵉᶠᶠ(t)` from the
#  compliance-space Dilute prediction (loading age t_0).
#  When the user requests interface-stiffness cracks (cf. the Python
#  benchmark proper), the next implementation iteration will add the
#  `(Rn, Rt)` ALV laws to `cracks_alv.jl`.
#
#  Usage  : julia --project scripts/41_fluage_echoes_cracks.jl
#  Output : scripts/figures/41_fluage_echoes_cracks.png
# =============================================================================

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."); io = devnull)

using MeanFieldHom
using TensND
using LinearAlgebra
using Printf
using Plots

# ─── Matrix law (Maxwell-like ageing relaxation) ───────────────────────────

# Same parameters as the Python script.
const k₀ = 5.0;  const μ₀ = 2.0
const k_inf = 3.0;  const μ_inf = 1.0
const τ = 1.0

const C₀_t = TensISO{3}(3 * k₀, 2 * μ₀)
const C_inf_t = TensISO{3}(3 * k_inf, 2 * μ_inf)

# Iso projectors as 6×6 Mandel templates.
const _, 𝕁₄, 𝕂₄ = TensND.iso_projectors(Val(3), Val(Float64))
const _J_M = MeanFieldHom.Viscoelasticity._tens_to_mandel66(𝕁₄)
const _K_M = MeanFieldHom.Viscoelasticity._tens_to_mandel66(𝕂₄)

# Build the matrix relaxation kernel C_inf + (C₀(1+0.2√t')-C_inf) exp(-(t-t')/τ).
function R_M(t, tp)
    factor = exp(-(t - tp) / τ)
    α0 = 3 * k₀ * (1 + 0.2 * sqrt(max(tp, 0.0)))
    β0 = 2 * μ₀ * (1 + 0.2 * sqrt(max(tp, 0.0)))
    α_inf = 3 * k_inf;  β_inf = 2 * μ_inf
    α = α_inf + (α0 - α_inf) * factor
    β = β_inf + (β0 - β_inf) * factor
    return α .* (_J_M ./ 3) .+ β .* (_K_M ./ 2)
    # = (α + 2β)/3 on the (1:3, 1:3) diagonal entries etc.
end
const law_M = ViscoLaw(R_M, :relaxation)

# ─── Crack contribution (penny, density d, normal e_3) ─────────────────────

function effective_creep_with_cracks(d, T)
    # Matrix-only relaxation, then compliance.
    R̃_M = MeanFieldHom.Viscoelasticity._trapezoidal_relaxation(law_M, T, 6)
    J̃_M = volterra_inverse(R̃_M; block_size = 6)
    # Penny crack contribution to compliance.
    crack = PennyCrack(1.0)
    H̃ = compliance_contribution_alv(crack, law_M, T)
    ΔJ̃ = delta_compliance_alv(crack, H̃, d)
    return J̃_M .+ ΔJ̃
end

# Effective uniaxial response J_E(t) from the compliance matrix:
#   apply unit longitudinal stress at every t and read E_xx.
function uniaxial_creep_curve(J̃::AbstractMatrix, n)
    S = zeros(eltype(J̃), 6 * n)
    @inbounds for i in 1:n
        S[6 * (i - 1) + 1] = 1.0
    end
    E = J̃ * S
    return [E[6 * (i - 1) + 1] for i in 1:n]
end

# ─── Plot ──────────────────────────────────────────────────────────────────

const N_TIMES = 50
const t0_v = range(0.0, 20.0; length = 3)
const dens_v = (0.0, 0.7)
const dens_colors = (:gray, :blue)

function build_grid(t0, n)
    return t0 .+ vcat(0.0, 10 .^ range(-2.0, log10(50.0 - t0); length = n))
end

plt = plot(layout = (1, 1), size = (1100, 700),
           title = "ALV penny cracks (dilute) — pure traction-free, no interface stiffness",
           xlabel = "t", ylabel = "Eₓₓ(t)",
           legend = :topleft)

for t0 in t0_v
    T = build_grid(t0, N_TIMES)
    for (d, col) in zip(dens_v, dens_colors)
        J̃ = effective_creep_with_cracks(d, T)
        n = length(T)
        E_curve = uniaxial_creep_curve(J̃, n)
        plot!(plt, T, E_curve; color = col,
              label = (t0 == 0.0 ? "ε³ᵈ = $d" : ""))
    end
end

mkpath(joinpath(@__DIR__, "figures"))
out = joinpath(@__DIR__, "figures", "41_fluage_echoes_cracks.png")
savefig(plt, out)

println("Saved : $out")
println()
println("Note : the Python `fluage_echoes_cracks.py` benchmark uses cracks WITH")
println("interface stiffness `(Rn(t,t'), Rt(t,t'))` and several schemes (MT, SC,")
println("PCW).  This Julia reproduction covers the **pure penny crack** in")
println("dilute mode only — interface-stiffness ALV cracks and crack-aware MT /")
println("SC scheme dispatch are scheduled for v0.6.2.")
