# # n-layer confocal spheroid: imperfect-interface effective conductivity
#
# Effective conductivity of a two-phase composite — an isotropic matrix
# reinforced by aligned [`LayeredSpheroid`](@ref) particles carrying a
# single imperfect thermal interface (Kapitza contact resistance) — as a
# function of a dimensionless interface (Biot-type) parameter. This is
# the benchmark configuration of
# [Kushch, Sevostianov & Belyaev (2015)](@cite kushch2015) ("Effective
# conductivity of spheroidal particle composite with imperfect
# interfaces"), also treated as the equivalent-inclusion application of
# [Barthélémy & Bignonnet (2020)](@cite barthelemyBignonnetIJES2020).
#
# The example demonstrates the central result of the equivalent-particle
# theory: **the full confocal series solution of the imperfect-interface
# spheroid is homogenization-equivalent to a single homogeneous,
# perfectly-bonded ellipsoid** carrying the anisotropic conductivity
# ```math
# \mathbf{k}^{eq} = \langle\mathbf B\rangle \cdot \langle\mathbf A\rangle^{-1}
# ```
# (Barthélémy & Bignonnet, eq. for ``k^{eq}``). We compute the effective
# conductivity **both ways** — through the layered spheroid and through
# its equivalent ellipsoid — for two mean-field schemes (Mori–Tanaka and
# Dilute) and check that the two curves coincide.
#
# The API exercised here is the conduction-only counterpart of the
# elastic [`LayeredSphere`](@ref)/[`RVE`](@ref) machinery:
# - [`LayeredSpheroid`](@ref)`(axis_radii, disk_radii, moduli; interfaces)`
#   — confocal geometry + per-layer isotropic conductivity;
# - [`MeanFieldHom.KapitzaInterface`](@ref)`(resistance)` — the "LC"
#   (low-conducting) interface, `[T] = ρ · qₙ`;
# - [`gradient_gradient_loc`](@ref) / [`flux_gradient_loc`](@ref) — the
#   volume-averaged concentration `⟨A⟩` and contribution `⟨B⟩` tensors
#   that define `kᵉ𝑞`;
# - `RVE` + `add_matrix!`/`add_phase!` + [`homogenize`](@ref).

import Pkg                                                          #jl
Pkg.activate(joinpath(@__DIR__, ".."); io = devnull)                 #jl

using MeanFieldHom
using TensND
using LinearAlgebra
using Printf
using Plots
gr()

# ## Setup
#
# Aligned prolate spheroids (aspect ratio `ω`, symmetry axis `ê₁`),
# volume fraction `c = 0.5`, conductivity contrast `k₁/kₘ = 10³`,
# isotropic matrix `kₘ = 1`. A single-layer particle (`N = 1`) has a
# single interface, at the particle/matrix boundary.

const km = 1.0
const k1 = 1.0e3
const c_frac = 0.5
const Nseries = 5
const axis = (1.0, 0.0, 0.0)          # prolate symmetry axis
const K1 = TensISO{3}(k1)
const Km = TensISO{3}(km)

# `k̃₁₁` is the effective conductivity along the particles' symmetry axis.

_axial(K) = K[1, 1]

# Build the layered-spheroid particle for aspect ratio `ω` and Kapitza
# resistance `ρ`.

function _particle(ω, ρ; N = Nseries)
    interface = MeanFieldHom.KapitzaInterface(ρ)
    return LayeredSpheroid((ω,), (1.0,), (K1,); interfaces = (interface,), Nseries = N, axis)
end

# ### Path 1 — full layered-spheroid series solution

function khom_series(ω, ρ, scheme; N = Nseries)
    s = _particle(ω, ρ; N)
    rve = RVE(:M)
    add_matrix!(rve, Ellipsoid(1.0, 1.0, 1.0), Dict(:K => Km))
    add_phase!(rve, :I, s, Dict(:K => K1); fraction = c_frac)
    return _axial(homogenize(rve, scheme, :K))
end

# ### Path 2 — equivalent homogeneous ellipsoid `kᵉ𝑞 = ⟨B⟩⟨A⟩⁻¹`
#
# `⟨A⟩` and `⟨B⟩` are transversely-isotropic (diagonal in the
# axial/transverse frame), so `kᵉ𝑞` is assembled component-wise and fed
# to an ordinary [`Ellipsoid`](@ref) of the same shape.

function k_eq_tensor(ω, ρ; N = Nseries)
    s = _particle(ω, ρ; N)
    A = gradient_gradient_loc(s, K1, Km)
    B = flux_gradient_loc(s, K1, Km)
    keq_axial = B[1, 1] / A[1, 1]
    keq_trans = B[2, 2] / A[2, 2]
    return TensND.TensTI{2}(keq_trans, keq_axial, axis)
end

function khom_equiv(ω, ρ, scheme; N = Nseries)
    Keq = k_eq_tensor(ω, ρ; N)
    rve = RVE(:M)
    add_matrix!(rve, Ellipsoid(1.0, 1.0, 1.0), Dict(:K => Km))
    add_phase!(rve, :I, Ellipsoid(ω, 1.0, 1.0), Dict(:K => Keq); fraction = c_frac)
    return _axial(homogenize(rve, scheme, :K))
end

# ## Interface-parameter sweep
#
# `h̃ = kₘ·b / ρ` (Biot-type number, with the equatorial radius `b = 1`):
# `h̃ → ∞` (small `ρ`) is the perfect-interface limit, `h̃ → 0`
# (large `ρ`) the fully insulated one.

const h_tilde = exp10.(range(-2, 3; length = 60))
const omegas = (3.0, 10.0)
const schemes = (MoriTanaka(), Dilute())
const scheme_names = ("Mori–Tanaka", "Dilute")

# Compute, for every `(ω, scheme)`, both the series and the
# equivalent-inclusion curves, and track the largest discrepancy.

series = Dict{Tuple{Float64, Int}, Vector{Float64}}()
equiv = Dict{Tuple{Float64, Int}, Vector{Float64}}()
max_gap = 0.0
for ω in omegas, (si, scheme) in enumerate(schemes)
    vs = similar(h_tilde)
    ve = similar(h_tilde)
    for (j, h) in enumerate(h_tilde)
        ρ = km / h                    # b = 1
        vs[j] = khom_series(ω, ρ, scheme)
        ve[j] = khom_equiv(ω, ρ, scheme)
    end
    series[(ω, si)] = vs
    equiv[(ω, si)] = ve
    global max_gap = max(max_gap, maximum(abs.(vs .- ve) ./ abs.(vs)))
end

println("Kushch-style imperfect-interface sweep — c = $c_frac, k₁/kₘ = $k1")
println("─"^70)
@printf "Largest relative gap  |k_series − k_equiv| / |k_series|  over the whole\n"
@printf "sweep (both schemes, both aspect ratios): %.2e\n" max_gap
println("→ the layered spheroid and its equivalent ellipsoid homogenize identically.")
println()

# ## Graphical output
#
# Solid lines: full layered-spheroid series. Open circles: equivalent
# homogeneous ellipsoid `kᵉ𝑞`. They overlie each other for every scheme
# and aspect ratio.

p = plot(
    xscale = :log10, xlabel = "h̃ = kₘ·b / ρ", ylabel = "k₁₁ᵉᶠᶠ / kₘ",
    title = "Effective axial conductivity vs. Kapitza interface parameter\n" *
        "(c=$c_frac, k₁/kₘ=$k1) — series vs. equivalent inclusion",
    legend = :topleft, size = (860, 580),
)
colors = (:steelblue, :darkorange)
styles = (:solid, :dash)
for (i, ω) in enumerate(omegas), (si, scheme) in enumerate(schemes)
    plot!(
        p, h_tilde, series[(ω, si)] ./ km;
        label = "series, ω=$(round(Int, ω)), $(scheme_names[si])",
        color = colors[i], linestyle = styles[si], lw = 2,
    )
    scatter!(
        p, h_tilde[1:4:end], equiv[(ω, si)][1:4:end] ./ km;
        label = "k_eq inclusion, ω=$(round(Int, ω)), $(scheme_names[si])",
        color = colors[i], markershape = :circle, markersize = 4,
        markerstrokewidth = 1, markeralpha = 0.0, markerstrokecolor = colors[i],
    )
end
p

const figdir = joinpath(@__DIR__, "figures")                        #jl
isdir(figdir) || mkdir(figdir)                                       #jl
figpath = joinpath(figdir, "32_spheroid_nlayers_conductivity.png")    #jl
savefig(p, figpath)                                                   #jl
display(p)                                                            #jl
@printf "\nSaved : %s\n" figpath                                      #jl
