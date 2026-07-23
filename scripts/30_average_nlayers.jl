# =============================================================================
#  30_average_nlayers.jl
#
#  Volume-averaged strain and stress localization tensors of an isotropic
#  n-layer composite sphere.  Mirrors the verification idea of
#  `tests/python/echoes_tests/average_nlayers.py` (random per-layer
#  moduli, random reference matrix) and adds a per-layer bar chart of
#  the bulk and shear localization factors `(α_k, β_k)`.
#
#  The MeanFieldHom.jl API used here:
#   * `LayeredSphere(radii, moduli)` — geometry + per-layer stiffness;
#   * `strain_strain_loc(sphere, C₀; layer=k)` — per-layer iso A_k;
#   * `stiffness_contribution(sphere, C₀)`    — size-independent N;
#   * `layer_volume_fraction(sphere, k)`      — f_k.
#
#  The whole-inclusion average <A_εε>_Ω = Σ_k f_k A_k (perfect
#  interfaces) is reconstructed and cross-checked against the dilute-
#  scheme identity  C_eff = C₀ + f·N  ⇔  N = <(C_k − C₀) : A_εε>.
# =============================================================================

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."); io = devnull)

using MeanFieldHom
using TensND
using Random
using Printf
using LinearAlgebra
using Plots

# ─── Helpers ─────────────────────────────────────────────────────────────────

# (E, ν) → (3K, 2μ) for direct TensISO{3} construction.
function _stiff_Enu(E::Real, ν::Real)
    K = E / (3 * (1 - 2ν))
    μ = E / (2 * (1 + ν))
    return TensISO{3}(3K, 2μ)
end

_iso_Kμ(C::TensISO{4, 3}) = (TensND.get_data(C)[1] / 3, TensND.get_data(C)[2] / 2)

# Convert layer fractions (positive, sum = 1 expected) to ascending radii
# given a target outer radius `R`.  Layer k occupies r_{k-1}..r_k with
# (r_k³ - r_{k-1}³) / R³ = f_k.
function _radii_from_fractions(f::AbstractVector{<:Real}, R::Real)
    @assert all(>(0), f) "layer fractions must be > 0"
    f_norm = f ./ sum(f)
    cum = zero(R)
    radii = similar(f_norm, typeof(R))
    for k in eachindex(f_norm)
        cum += f_norm[k] * R^3
        radii[k] = cbrt(cum)
    end
    return Tuple(radii)
end

# ─── Setup — random n-layer sphere + random reference ────────────────────────

Random.seed!(20260426)

const n = 10
const R = 5.0

# Random layer fractions (Dirichlet-ish: independent uniforms re-normalized).
const fractions = rand(n) .+ 0.05
const fractions_n = fractions ./ sum(fractions)

# Random per-layer moduli (E, ν).
const E_rand = 1.0 .+ 9.0 .* rand(n)
const ν_rand = 0.05 .+ 0.4 .* rand(n)
const C_layers = ntuple(k -> _stiff_Enu(E_rand[k], ν_rand[k]), n)

# Random reference (matrix) modulus.
const C_ref = _stiff_Enu(1.0 + 9.0 * rand(), 0.05 + 0.4 * rand())

const radii = _radii_from_fractions(fractions_n, R)
const sphere = LayeredSphere(radii, C_layers)

# ─── Print summary ───────────────────────────────────────────────────────────

println("Random n-layer sphere — n = $n, R = $R")
println("─"^78)
println("Per-layer (E, ν) :")
for k in 1:n
    @printf "  k=%2d  E = %6.3f  ν = %5.3f\n" k E_rand[k] ν_rand[k]
end
println()

K_ref, μ_ref = _iso_Kμ(C_ref)
@printf "Reference matrix : K_ref = %.5f, μ_ref = %.5f\n" K_ref μ_ref
println()

println("Geometry :")
@printf "  outer radius          : %.5f\n" radii[end]
@printf "  Σ layer_fraction      : %.10f  (should be 1)\n" sum(
    layer_volume_fraction(sphere, k) for k in 1:n
)
println("  layer_radius          : ", join(map(r -> @sprintf("%.4f", r), radii), ", "))
println(
    "  layer_volume_fraction : ", join(
        map(k -> @sprintf("%.5f", layer_volume_fraction(sphere, k)), 1:n), ", "
    )
)
println()

# ─── Per-layer localization tensors A_k (iso 4-tensor) ───────────────────────

A_per_layer = ntuple(k -> strain_strain_loc(sphere, C_ref; layer = k), n)
αβ = ntuple(k -> TensND.get_data(A_per_layer[k]), n)
α_k = [αβ[k][1] for k in 1:n]   # bulk localization
β_k = [αβ[k][2] for k in 1:n]   # shear localization

println("Per-layer localization factors :")
@printf "  %-4s  %-12s  %-12s\n" "k" "α_k (bulk)" "β_k (shear)"
for k in 1:n
    @printf "  %2d    %+10.6f    %+10.6f\n" k α_k[k] β_k[k]
end
println()

# Whole-inclusion average <A_εε>_Ω = Σ_k f_k A_k  (perfect interfaces).
A_whole_α = sum(layer_volume_fraction(sphere, k) * α_k[k] for k in 1:n)
A_whole_β = sum(layer_volume_fraction(sphere, k) * β_k[k] for k in 1:n)
@printf "Whole-inclusion average :  α = %+12.8f   β = %+12.8f\n" A_whole_α A_whole_β

# Cross-check: stiffness_contribution N satisfies  N = <(C_k − C₀) : A_k>.
N_tensor = stiffness_contribution(sphere, C_ref)
N_α, N_β = TensND.get_data(N_tensor)
N_α_check = sum(
    layer_volume_fraction(sphere, k) *
        (TensND.get_data(C_layers[k])[1] - TensND.get_data(C_ref)[1]) *
        α_k[k] for k in 1:n
)
N_β_check = sum(
    layer_volume_fraction(sphere, k) *
        (TensND.get_data(C_layers[k])[2] - TensND.get_data(C_ref)[2]) *
        β_k[k] for k in 1:n
)
@printf "stiffness_contribution :  N_α = %+12.6f  (rec. %+12.6f, diff=%.2e)\n" N_α N_α_check abs(N_α - N_α_check)
@printf "                         N_β = %+12.6f  (rec. %+12.6f, diff=%.2e)\n" N_β N_β_check abs(N_β - N_β_check)
println()

# ─── Reset reference to last-layer modulus (mirror Python l. 59-60) ──────────
println("Resetting reference to C_layers[end] :")
C_ref2 = C_layers[end]
A_layer_last = strain_strain_loc(sphere, C_ref2; layer = n)
α_last, β_last = TensND.get_data(A_layer_last)
@printf "  layer_eE(n) with C_ref = C_layers[n] :  α = %+10.6f  β = %+10.6f\n" α_last β_last
println("  (when reference = layer modulus, the layer is invisible to itself.)")
println()

# ─── Graphical output ────────────────────────────────────────────────────────

const figdir = joinpath(@__DIR__, "figures")
isdir(figdir) || mkdir(figdir)

p1 = bar(
    1:n, α_k;
    xlabel = "layer k", ylabel = "α_k (bulk)",
    title = "Bulk localization per layer",
    legend = false, color = :steelblue
)
hline!(
    p1, [A_whole_α]; lw = 2, color = :red, linestyle = :dash,
    label = "<α> = $(round(A_whole_α; digits = 4))"
)
p1 = plot!(p1; legend = :topright)

p2 = bar(
    1:n, β_k;
    xlabel = "layer k", ylabel = "β_k (shear)",
    title = "Shear localization per layer",
    legend = false, color = :darkorange
)
hline!(
    p2, [A_whole_β]; lw = 2, color = :red, linestyle = :dash,
    label = "<β> = $(round(A_whole_β; digits = 4))"
)
p2 = plot!(p2; legend = :topright)

p3 = bar(
    1:n, [layer_volume_fraction(sphere, k) for k in 1:n];
    xlabel = "layer k", ylabel = "f_k",
    title = "Layer volume fractions",
    legend = false, color = :seagreen
)

p_full = plot(
    p1, p2, p3; layout = (1, 3), size = (1500, 450),
    plot_title = "n-layer sphere ($n layers, R=$R) — average localizations"
)

figpath = joinpath(figdir, "30_average_nlayers.png")
savefig(p_full, figpath)
@printf "\nSaved : %s\n" figpath
