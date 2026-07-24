# # n-layer confocal spheroid: equivalent conductivity of an imperfect-interface particle
#
# Barthélémy & Bignonnet (IJES 2020, §4) define the **equivalent
# particle** of a spheroidal inclusion with an imperfect interface as
# the homogeneous, perfectly-bonded spheroid carrying the same
# volume-averaged concentration tensors:
# ```math
# \uu k^{eq} = \langle\uu B\rangle_{\cal E} \cdot
#              \langle\uu A\rangle_{\cal E}^{-1}
# \qquad\text{(eq:defkeqAB)}
# ```
# — a genuinely SIZE-dependent quantity (unlike a perfect-interface
# homogeneous inclusion's concentration, which depends on shape only):
# the interface enters through a surface-to-volume ratio, so `kᵉ𝑞`
# interpolates between the perfect-interface limit (large particles,
# where the interface becomes negligible) and the fully-decoupled limit
# (small particles / strong interface resistance).
#
# `echoes_cpp/tests/python/spheroid_nlayers/spheroid_nlayers_keq.py`
# goes on to translate `kᵉ𝑞` into an "equivalent RADIUS" under several
# geometric conventions (confocal / similar-shape / equal-surface
# rescaling of a reference perfect-interface particle) via a separate
# closed-form APPROXIMATE model (Bonfoh-type). This script reproduces
# the foundational, EXACT quantity `kᵉ𝑞` itself directly from
# [`LayeredSpheroid`](@ref)'s confocal transfer-matrix solution — using
# no formula beyond what `test/LayeredSpheroids/test_conductivity.jl`
# already validates against the classical Eshelby limit and the
# fully-insulated-core limit.

import Pkg                                                          #jl
Pkg.activate(joinpath(@__DIR__, ".."); io = devnull)                 #jl

using MeanFieldHom
using TensND
using Printf
using Plots
gr()

const km = 1.0
const k1 = 10.0
const a_phys = 1.0   # physical axis semi-axis (fixed) — kᵉ𝑞 varies with ω at fixed absolute size

# ## Equivalent conductivity `kᵉ𝑞 = ⟨B⟩·⟨A⟩⁻¹` (axial / transverse components)

function k_eq(ω, interface; Nseries = 8)
    K1 = TensISO{3}(k1); K0 = TensISO{3}(km)
    axis = ω > 1 ? (1.0, 0.0, 0.0) : (0.0, 0.0, 1.0)
    s = LayeredSpheroid((a_phys,), (a_phys / ω,), (K1,); interfaces = (interface,), Nseries, axis)
    A = gradient_gradient_loc(s, K1, K0)
    B = flux_gradient_loc(s, K1, K0)
    i_axial = ω > 1 ? 1 : 3
    i_trans = ω > 1 ? 2 : 1
    return B[i_axial, i_axial] / A[i_axial, i_axial], B[i_trans, i_trans] / A[i_trans, i_trans]
end

# ## Sweep over aspect ratio, for two Kapitza resistances
#
# `ρ → 0` recovers the perfect-interface Eshelby value at every `ω`
# (shown as a dotted reference curve); larger `ρ` pulls `kᵉ𝑞` toward
# `kₘ` (the interface increasingly decouples the particle from the
# matrix).

const omegas_prolate = exp10.(range(0.01, 2; length = 30))
const omegas_oblate = 1 ./ omegas_prolate
const rhos = (0.05, 0.5)

function _sweep(omegas, ρ)
    a = Float64[]; t = Float64[]
    for ω in omegas
        ka, kt = k_eq(ω, MeanFieldHom.KapitzaInterface(ρ))
        push!(a, ka); push!(t, kt)
    end
    return a, t
end

function _sweep_perfect(omegas)
    a = Float64[]; t = Float64[]
    for ω in omegas
        ka, kt = k_eq(ω, MeanFieldHom.PerfectInterface())
        push!(a, ka); push!(t, kt)
    end
    return a, t
end

a_perfect_pro, t_perfect_pro = _sweep_perfect(omegas_prolate)
a_perfect_ob, t_perfect_ob = _sweep_perfect(omegas_oblate)

println("Equivalent conductivity kᵉ𝑞 vs. aspect ratio — a=$a_phys, k₁/kₘ=$(k1 / km)")
println("─"^70)
@printf "Perfect-interface reference (shape-independent): kᵉ𝑞 ≡ k₁ = %.4f, in every direction\n" k1
for ρ in rhos
    ka, kt = _sweep(omegas_prolate, ρ)
    @printf "ρ=%.2f   ω=%6.2f : k_eq_axial/km=%.4f  k_eq_trans/km=%.4f  (interface pulls both toward kₘ=%.1f)\n" ρ omegas_prolate[15] ka[15] kt[15] km
end
println()

# ## Graphical output

p1 = plot(
    xscale = :log10, xlabel = "ω (prolate)", ylabel = "k_eq / km",
    title = "Prolate", legend = :topleft,
)
p2 = plot(
    xscale = :log10, xlabel = "ω (oblate)", ylabel = "k_eq / km",
    title = "Oblate", legend = :topleft,
)
plot!(p1, omegas_prolate, a_perfect_pro ./ km; label = "axial (perfect)", color = :black, linestyle = :dot, lw = 2)
plot!(p1, omegas_prolate, t_perfect_pro ./ km; label = "trans (perfect)", color = :gray, linestyle = :dot, lw = 2)
plot!(p2, omegas_oblate, a_perfect_ob ./ km; label = "axial (perfect)", color = :black, linestyle = :dot, lw = 2)
plot!(p2, omegas_oblate, t_perfect_ob ./ km; label = "trans (perfect)", color = :gray, linestyle = :dot, lw = 2)
colors = (:steelblue, :darkorange)
for (i, ρ) in enumerate(rhos)
    ka_pro, kt_pro = _sweep(omegas_prolate, ρ)
    ka_ob, kt_ob = _sweep(omegas_oblate, ρ)
    plot!(p1, omegas_prolate, ka_pro ./ km; label = "axial, ρ=$ρ", color = colors[i], lw = 2)
    plot!(p1, omegas_prolate, kt_pro ./ km; label = "trans, ρ=$ρ", color = colors[i], linestyle = :dash, lw = 2)
    plot!(p2, omegas_oblate, ka_ob ./ km; label = "axial, ρ=$ρ", color = colors[i], lw = 2)
    plot!(p2, omegas_oblate, kt_ob ./ km; label = "trans, ρ=$ρ", color = colors[i], linestyle = :dash, lw = 2)
end
p_full = plot(
    p1, p2; layout = (1, 2), size = (1200, 500),
    plot_title = "Equivalent particle conductivity — Kapitza interface, k₁/kₘ=$(k1 / km)",
)
p_full

const figdir = joinpath(@__DIR__, "figures")                        #jl
isdir(figdir) || mkdir(figdir)                                       #jl
figpath = joinpath(figdir, "34_spheroid_equivalent_conductivity.png") #jl
savefig(p_full, figpath)                                              #jl
display(p_full)                                                       #jl
@printf "\nSaved : %s\n" figpath                                      #jl
