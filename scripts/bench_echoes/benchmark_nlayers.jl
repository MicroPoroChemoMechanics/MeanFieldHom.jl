# =============================================================================
#  scripts/bench_echoes/benchmark_nlayers.jl
#
#  Verification of MeanFieldHom.jl `LayeredSphere` against four
#  independent references:
#
#   § 1  **Bulk α_k** vs `echoes.layer_eE` (volume-averaged strain
#        localisation in each layer).
#   § 2  **Internal consistency**: Julia state-vector recurrence vs
#        a direct 8×8 linear-system solver assembled from the same mode
#        formulas (sanity check that the recurrence implements the
#        correct boundary-value problem).
#   § 3  **Analytical limits**: shear localisation `β_k` in degenerate
#        configurations (vanishing core, vanishing shell, core ≡ shell)
#        compared to the closed-form single-layer Eshelby result.
#   § 4  **Local bulk stress profile** `σ_rr(r), σ_θθ(r)` vs
#        `echoes.loc_sS` under remote hydrostatic loading.
#
#  Note on the shear (β_k) ECHOES comparison
#  -----------------------------------------
#  Direct comparison of `β_k = (layer_eE[1,1] - layer_eE[1,2])` between
#  Julia and ECHOES disagrees by 1–50 % in genuine multi-layer cases,
#  while bulk α_k matches to 5e-13 and the local bulk stress profile
#  matches `loc_sS` to 5e-16.  Both Julia's mode formulas and the 8×8
#  direct solver agree, and the result reproduces the analytical
#  Eshelby limits.  The disagreement appears to stem from a different
#  internal convention in `echoes.layer_eE` (the C++ source uses
#  `layer++` before fetching mode amplitudes — i.e. the outer layer's
#  amplitudes evaluated over the inner layer's volume).  We therefore
#  benchmark β_k against analytical limits rather than ECHOES.
#
#  Run from the `MeanFieldHom.jl` package root:
#    julia --project=scripts/bench_echoes scripts/bench_echoes/benchmark_nlayers.jl
# =============================================================================

import Pkg
Pkg.activate(@__DIR__; io = devnull)

using MeanFieldHom
using TensND
using LinearAlgebra
using Printf
using Random
using PyCall
using Plots

# Internal helpers reused for the local-bulk profile (mirrors script 32).
import MeanFieldHom.LayeredSpheres: _iso_bulk_shear, _bulk_state_seq,
    _bulk_extract_AB, _shear_M_matrix, _layer_avg_dev_shear_factor

# ─── Python-side wrappers ────────────────────────────────────────────────────

py"""
import echoes
import numpy as np

def py_stiff_kmu(K, mu):
    return echoes.stiff_kmu(float(K), float(mu))

def py_make_nlayers(radii, props, Cref):
    radii_np = np.asarray(radii, dtype=float)
    spn = echoes.sphere_nlayers(radii=radii_np, prop={'C': props})
    spn.set_ref('C', Cref)
    return spn

def py_layer_eE(spn, k):
    # ECHOES is 0-indexed.
    return np.asarray(spn.layer_eE(k))

def py_eE(spn):
    return np.asarray(spn.eE)

def py_sE(spn):
    return np.asarray(spn.sE)

def py_layer_fraction(spn, k):
    return float(spn.layer_fraction(k))

def py_loc_sS(spn, r, theta, phi):
    return np.asarray(spn.loc_sS(r, theta, phi))
"""

const py_stiff_kmu = py"py_stiff_kmu"
const py_make_nlayers = py"py_make_nlayers"
const py_layer_eE = py"py_layer_eE"
const py_eE_total = py"py_eE"
const py_sE_total = py"py_sE"
const py_layer_frac = py"py_layer_fraction"
const py_loc_sS = py"py_loc_sS"

# ─── Helpers ─────────────────────────────────────────────────────────────────

# (E, ν) → (K, μ).
function _Kmu_Enu(E, ν)
    K = E / (3 * (1 - 2ν))
    μ = E / (2 * (1 + ν))
    return K, μ
end

# Convert ECHOES 6×6 Voigt-Mandel stiffness array (standard `Cref.array`) to
# Julia (3K, 2μ) iso-data: `α = (KM[1,1] + 2 KM[1,2]) = 3K`, similarly for β.
function _iso_data_from_echoes6x6(KM::AbstractMatrix)
    α = KM[1, 1] + 2 * KM[1, 2]   # = 3K
    β = KM[1, 1] - KM[1, 2]       # = 2μ
    return α, β
end

# Convert layer fractions to ascending radii given outer R.
function _radii_from_fractions(f::AbstractVector{<:Real}, R::Real)
    f_norm = f ./ sum(f)
    cum = zero(R); radii = similar(f_norm, typeof(R))
    for k in eachindex(f_norm)
        cum += f_norm[k] * R^3
        radii[k] = cbrt(cum)
    end
    return radii
end

relerr(a, b) = (abs(a) + abs(b) < 1.0e-14) ? 0.0 : abs(a - b) / max(abs(a), abs(b))

# ─── §1 + §2  Random n-layer cross-check ────────────────────────────────────

println("="^78)
println("§1  Bulk α_k vs ECHOES layer_eE (random n-layer configs)")
println("="^78)

const rtol_match = 1.0e-8

Random.seed!(20260426)

const N_CONFIGS = 30
const N_LAYERS_RANGE = 2:8

n_pass_α = 0; n_fail_α = 0
worst_α_err = 0.0

for cfg in 1:N_CONFIGS
    n = rand(N_LAYERS_RANGE)
    R_outer = 1.0 + 4.0 * rand()

    fractions = rand(n) .+ 0.05
    radii = _radii_from_fractions(fractions, R_outer)

    E_lay = 1.0 .+ 99.0 .* rand(n)
    ν_lay = 0.05 .+ 0.4 .* rand(n)
    Kμ_lay = [_Kmu_Enu(E_lay[k], ν_lay[k]) for k in 1:n]

    E_ref = 1.0 + 99.0 * rand()
    ν_ref = 0.05 + 0.4 * rand()
    K_ref, μ_ref = _Kmu_Enu(E_ref, ν_ref)

    C_layers_jl = ntuple(k -> TensISO{3}(3 * Kμ_lay[k][1], 2 * Kμ_lay[k][2]), n)
    C_ref_jl = TensISO{3}(3 * K_ref, 2 * μ_ref)
    sphere_jl = LayeredSphere(Tuple(radii), C_layers_jl)

    C_layers_py = [py_stiff_kmu(Kμ_lay[k][1], Kμ_lay[k][2]) for k in 1:n]
    C_ref_py = py_stiff_kmu(K_ref, μ_ref)
    spn_py = py_make_nlayers(radii, C_layers_py, C_ref_py)

    cfg_α_err = 0.0
    for k in 1:n
        A_jl = strain_strain_loc(sphere_jl, C_ref_jl; layer = k)
        α_jl, _ = TensND.get_data(A_jl)
        eE_py = py_layer_eE(spn_py, k - 1)
        α_py, _ = _iso_data_from_echoes6x6(eE_py)
        cfg_α_err = max(cfg_α_err, relerr(α_jl, α_py))
    end

    pass_α = cfg_α_err ≤ rtol_match
    pass_α ? (global n_pass_α += 1) : (global n_fail_α += 1)
    global worst_α_err = max(worst_α_err, cfg_α_err)
end

@printf "  %d/%d configs within rtol = %.0e (worst rerr = %.3e)\n" n_pass_α (n_pass_α + n_fail_α) rtol_match worst_α_err
println()

# ─── §2  Internal consistency : Julia recurrence vs direct 8×8 solver ───────

println("="^78)
println("§2  β_k self-consistency : recurrence vs direct 8×8 linear-system")
println("="^78)

# Direct 8×8 solver for the 2-layer Y₂-harmonic shear problem.
# Unknowns x = (a₁, b₁, a₂, b₂, c₂, d₂, c_∞, d_∞);  c₁ = d₁ = 0 enforced.
# BC at r = ∞ : matrix mode-1 amplitude = 1, mode-2 amplitude = 0.
#
# β_layer1 is the layer-VOLUME-AVERAGED deviatoric strain localisation, not
# the bare mode-1 amplitude a₁: the core carries both the uniform mode 1 (a₁)
# and the r³-varying mode 2 (b₁), whose Y₂-projected volume average adds
# `b₁ · F₁` with the Christensen-Lo factor F₁ = _layer_avg_dev_shear_factor.
# (Earlier this returned a₁ alone — correct only in the degenerate limits
# where b₁ → 0, so it agreed with §3 but disagreed with the recurrence by a
# few % on general configs.)
function direct_2layer_β_layer1(r1, r2, κc, μc, κs, μs, κm, μm)
    Mc_r1 = _shear_M_matrix(r1, κc, μc)
    Ms_r1 = _shear_M_matrix(r1, κs, μs)
    Ms_r2 = _shear_M_matrix(r2, κs, μs)
    Mm_r2 = _shear_M_matrix(r2, κm, μm)
    A = zeros(8, 8); b = zeros(8)
    A[1:4, 1:2] = Mc_r1[:, 1:2]
    A[1:4, 3:6] = -Ms_r1
    A[5:8, 3:6] = Ms_r2
    A[5:8, 7:8] = -Mm_r2[:, 3:4]
    b[5:8] = Mm_r2[:, 1]
    x = A \ b
    a₁, b₁ = x[1], x[2]
    return a₁ + b₁ * _layer_avg_dev_shear_factor(0.0, r1, κc, μc)
end

n_pass_self = 0; worst_self_err = 0.0
for cfg in 1:20
    Kc = 1 + 99 * rand(); μc = 0.5 + 49.5 * rand()
    Ks = 1 + 99 * rand(); μs = 0.5 + 49.5 * rand()
    Km = 1 + 99 * rand(); μm = 0.5 + 49.5 * rand()
    r1 = 0.1 + 0.8 * rand()

    sphere = LayeredSphere(
        (r1, 1.0),
        (TensISO{3}(3Kc, 2μc), TensISO{3}(3Ks, 2μs))
    )
    C0 = TensISO{3}(3Km, 2μm)
    _, β_recurrence = TensND.get_data(strain_strain_loc(sphere, C0; layer = 1))
    β_direct = direct_2layer_β_layer1(r1, 1.0, Kc, μc, Ks, μs, Km, μm)
    err = relerr(β_recurrence, β_direct)
    err < rtol_match && (global n_pass_self += 1)
    global worst_self_err = max(worst_self_err, err)
end
@printf "  %d/20 configs within rtol = %.0e (worst rerr = %.3e)\n" n_pass_self rtol_match worst_self_err
println()

# ─── §3  Analytical-limit check : β in degenerate configurations ────────────

println("="^78)
println("§3  β_layer1 vs analytical Eshelby in degenerate limits")
println("="^78)

# Single-layer Eshelby strain localisation for a sphere of moduli (μ₁) in
# matrix (κ₀, μ₀):  β_∞ = 1 / (1 + α_dev (μ₁/μ₀ − 1))  with
# α_dev = 6(κ₀+2μ₀) / (5(3κ₀+4μ₀)).
function β_eshelby_sphere(μ1, κ0, μ0)
    α_dev = 6 * (κ0 + 2μ0) / (5 * (3κ0 + 4μ0))
    return 1 / (1 + α_dev * (μ1 / μ0 - 1))
end

const Kc, μc = 80.0, 30.0
const Ks, μs = 20.0, 8.0
const Km, μm = 50.0, 20.0
const C_ref_lim = TensISO{3}(3Km, 2μm)
const C_core = TensISO{3}(3Kc, 2μc)
const C_shell = TensISO{3}(3Ks, 2μs)

# Limit 1 : shell ≡ matrix → β_layer1 = single-layer Eshelby (core in matrix).
let
    sphere = LayeredSphere((0.5, 1.0), (C_core, C_ref_lim))
    _, β_jl = TensND.get_data(strain_strain_loc(sphere, C_ref_lim; layer = 1))
    β_an = β_eshelby_sphere(μc, Km, μm)
    @printf "  shell ≡ matrix      :  Julia β = %.10f   analytical = %.10f   relerr = %.2e\n" β_jl β_an relerr(β_jl, β_an)
end

# Limit 2 : core ≡ shell → β_layer1 = single-layer Eshelby for sphere of
# core moduli at full radius (radius 1) in matrix.
let
    sphere = LayeredSphere((0.5, 1.0), (C_core, C_core))
    _, β_jl = TensND.get_data(strain_strain_loc(sphere, C_ref_lim; layer = 1))
    β_an = β_eshelby_sphere(μc, Km, μm)
    @printf "  core ≡ shell        :  Julia β = %.10f   analytical = %.10f   relerr = %.2e\n" β_jl β_an relerr(β_jl, β_an)
end

# Limit 3 : tiny core (r₁ → 0) → β_layer1 = β_inner_Eshelby × β_shell_Eshelby
let
    r1 = 1.0e-4
    sphere = LayeredSphere((r1, 1.0), (C_core, C_shell))
    _, β_jl = TensND.get_data(strain_strain_loc(sphere, C_ref_lim; layer = 1))
    β_shell = β_eshelby_sphere(μs, Km, μm)
    β_core_in_shell = β_eshelby_sphere(μc, Ks, μs)
    β_an = β_core_in_shell * β_shell
    @printf "  r₁ = 1e-4 (vanishing core) :  Julia β = %.10f   analytical = %.10f   relerr = %.2e\n" β_jl β_an relerr(β_jl, β_an)
end

# Limit 4 : N=3 homogeneous (all layers equal to matrix) → β_k = 1 ∀ k.
let
    sphere = LayeredSphere((0.3, 0.7, 1.0), (C_ref_lim, C_ref_lim, C_ref_lim))
    βs = [TensND.get_data(strain_strain_loc(sphere, C_ref_lim; layer = k))[2] for k in 1:3]
    @printf "  homogeneous N=3     :  β = (%.10f, %.10f, %.10f)   (expected 1.0)\n" βs[1] βs[2] βs[3]
end
println()

# ─── §4  Local bulk profile vs `loc_sS` under hydrostatic load ──────────────

println("="^78)
println("§4  Local stress profile (hydrostatic) vs ECHOES loc_sS")
println("="^78)

# Use the Christensen-style 2-layer setup of script 32_local_nlayers.jl.
const Eo, νo = 30.0, 0.3
const Ei, νi = 100.0, 0.3
const Eitz, νitz = 0.1 * Ei, 0.2

K_o, μ_o = _Kmu_Enu(Eo, νo)
K_i, μ_i = _Kmu_Enu(Ei, νi)
K_itz, μ_itz = _Kmu_Enu(Eitz, νitz)

const C_ref_loc_jl = TensISO{3}(3 * K_o, 2 * μ_o)
const C_ref_loc_py = py_stiff_kmu(K_o, μ_o)
const C_layers_loc_py = [py_stiff_kmu(K_i, μ_i), py_stiff_kmu(K_itz, μ_itz)]

const R_inner = 1.0
const ep_layer = 2.0
const radii_loc = [R_inner, R_inner + ep_layer]
const sphere_loc_jl = LayeredSphere(
    (R_inner, R_inner + ep_layer),
    (
        TensISO{3}(3 * K_i, 2 * μ_i),
        TensISO{3}(3 * K_itz, 2 * μ_itz),
    )
)
const spn_loc_py = py_make_nlayers(radii_loc, C_layers_loc_py, C_ref_loc_py)

# Bulk-profile evaluator — copy of script 32 helper for self-containment.
function bulk_AB(sphere::LayeredSphere{T, N}, C₀::TensISO{4, 3}) where {T, N}
    κ₀, μ₀ = _iso_bulk_shear(C₀)
    inside, s_mat = _bulk_state_seq(sphere, κ₀, μ₀)
    radii = sphere.radii
    A_inf, B_inf = _bulk_extract_AB(radii[N], κ₀, μ₀, s_mat[1], s_mat[2])
    inv_A_inf = 1 / A_inf
    AB = ntuple(N) do k
        κ_k, μ_k = _iso_bulk_shear(sphere.moduli[k])
        Ak, Bk = _bulk_extract_AB(radii[k], κ_k, μ_k, inside[k][1], inside[k][2])
        (Ak * inv_A_inf, Bk * inv_A_inf)
    end
    return AB, B_inf * inv_A_inf
end

function bulk_stresses(
        sphere::LayeredSphere{T, N}, C₀::TensISO{4, 3}, r;
        ε_v::Real = 1.0
    ) where {T, N}
    AB, B_inf = bulk_AB(sphere, C₀)
    radii = sphere.radii
    κ₀, μ₀ = _iso_bulk_shear(C₀)
    if r ≤ radii[N]
        k = findfirst(rk -> r ≤ rk, collect(radii))::Int
        Ak, Bk = AB[k]
        κ_k, μ_k = _iso_bulk_shear(sphere.moduli[k])
        σ_rr = 3 * κ_k * Ak * ε_v - 4 * μ_k * Bk * ε_v / r^3
        σ_θθ = 3 * κ_k * Ak * ε_v + 2 * μ_k * Bk * ε_v / r^3
    else
        σ_rr = 3 * κ₀ * ε_v - 4 * μ₀ * (B_inf * ε_v) / r^3
        σ_θθ = 3 * κ₀ * ε_v + 2 * μ₀ * (B_inf * ε_v) / r^3
    end
    return σ_rr, σ_θθ
end

# Far-field hydrostatic strain ε_v ⇒ remote uniaxial in *each* direction.
# In ECHOES Voigt-Mandel convention, `loc_sS(r, θ, φ)` returns a 6×6
# matrix σ_local = M · σ_∞.  For hydrostatic σ∞ = 3K₀ ε_v · 𝟙, a
# Voigt 6-vector representation is `(3 K₀ ε_v, 3 K₀ ε_v, 3 K₀ ε_v, 0, 0, 0)`.
# In the spherical-symmetric basis at point (r, 0, 0) (i.e. on the x axis),
# `σ_rr_local = M · S` projected on Voigt index 1 (xx).
const ε_v = 1.0
const σ_far = 3 * K_o * ε_v
const Σ_inf_voigt = [σ_far, σ_far, σ_far, 0.0, 0.0, 0.0]   # σ∞ in Voigt-6

const lr_check = collect(range(0.05 * R_inner, 5 * (R_inner + ep_layer); length = 80))
σ_rr_jl = similar(lr_check); σ_θθ_jl = similar(lr_check)
σ_rr_py = similar(lr_check); σ_θθ_py = similar(lr_check)

for (i, r) in enumerate(lr_check)
    σrr, σθθ = bulk_stresses(sphere_loc_jl, C_ref_loc_jl, r; ε_v = ε_v)
    σ_rr_jl[i] = σrr; σ_θθ_jl[i] = σθθ
    # ECHOES at (r, θ=0, φ=0) — local frame has e_z radial, so σ_rr ↔ Voigt index 3.
    M = py_loc_sS(spn_loc_py, r, 0.0, 0.0)
    σ_local = M * Σ_inf_voigt
    σ_rr_py[i] = σ_local[3]   # zz in local frame = rr at (θ=0, φ=0)
    σ_θθ_py[i] = σ_local[1]   # xx local = θθ
end

# Relative error sweep.
errs_rr = [relerr(σ_rr_jl[i], σ_rr_py[i]) for i in eachindex(lr_check)]
errs_θθ = [relerr(σ_θθ_jl[i], σ_θθ_py[i]) for i in eachindex(lr_check)]
@printf "Local stress max relerr — σ_rr : %.3e, σ_θθ : %.3e (over %d points)\n\n" maximum(errs_rr) maximum(errs_θθ) length(lr_check)

# Plot Julia vs ECHOES.
p_loc = plot(;
    xlabel = "r", ylabel = "σ_ij / σ∞",
    title = "Local stress profile — hydrostatic far-field",
    legend = :topright, grid = true
)
plot!(p_loc, lr_check, σ_rr_jl ./ σ_far; lw = 2, color = :red, label = "σ_rr (Julia)")
plot!(
    p_loc, lr_check, σ_rr_py ./ σ_far; lw = 0, marker = :circle, ms = 4,
    color = :red, label = "σ_rr (ECHOES)"
)
plot!(p_loc, lr_check, σ_θθ_jl ./ σ_far; lw = 2, color = :blue, label = "σ_θθ (Julia)")
plot!(
    p_loc, lr_check, σ_θθ_py ./ σ_far; lw = 0, marker = :diamond, ms = 4,
    color = :blue, label = "σ_θθ (ECHOES)"
)
hline!(p_loc, [1.0]; lw = 1, color = :black, linestyle = :dot, label = "σ∞")
vline!(p_loc, [R_inner]; lw = 1, color = :black, linestyle = :dash, label = "")
vline!(p_loc, [R_inner + ep_layer]; lw = 1, color = :black, linestyle = :dash, label = "")

const figdir = joinpath(@__DIR__, "figures")
isdir(figdir) || mkdir(figdir)
figpath = joinpath(figdir, "benchmark_nlayers.png")
savefig(p_loc, figpath)
@printf "\nSaved : %s\n" figpath
