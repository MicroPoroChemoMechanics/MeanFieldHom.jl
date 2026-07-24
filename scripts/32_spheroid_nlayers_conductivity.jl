# # n-layer confocal spheroid: imperfect-interface effective conductivity
#
# Effective conductivity of a Mori-Tanaka composite reinforced by
# [`LayeredSpheroid`](@ref) inclusions carrying a single imperfect
# thermal interface (Kapitza resistance or surface conductance), as a
# function of a dimensionless interface (Biot-type) parameter — the
# physical setting studied by
# [Kushch, Sevostianov & Belyaev (2015)](@cite kushch2015)
# ("Effective conductivity of spheroidal particle composite with
# imperfect interfaces") and mirroring the validation idea of
# `echoes_cpp/tests/python/spheroid_nlayers/spheroid_nlayers_test_Kushch.py`.
#
# The API exercised here is the conduction-only counterpart of the
# elastic [`LayeredSphere`](@ref)/[`RVE`](@ref) machinery:
# - [`LayeredSpheroid`](@ref)`(axis_radii, disk_radii, moduli; interfaces)`
#   — confocal geometry + per-layer isotropic conductivity;
# - [`MeanFieldHom.KapitzaInterface`](@ref)`(resistance)` — the "LC"
#   (low-conducting) interface, `[T] = ρ · qₙ`;
# - [`MeanFieldHom.SurfaceConductiveInterface`](@ref)`(conductance)` —
#   the "HC" (highly-conducting) interface, `[qₙ] = -divₛ(β ∇ₛT)`;
# - `RVE` + `add_matrix!`/`add_phase!` + [`homogenize`](@ref)`(rve,
#   MoriTanaka(), :K)`.
#
# Two physically exact limits anchor the interface-parameter sweep and
# double as a regression check (already verified to machine/near-machine
# precision in `test/LayeredSpheroids/test_conductivity.jl`):
# - `ρ → 0` (LC) or `β → 0` (HC) recover the **perfect interface**;
# - `ρ → ∞` (LC) recovers the **fully insulated** (impermeable) core,
#   identical to a vanishing-conductivity homogeneous inclusion.

import Pkg                                                          #jl
Pkg.activate(joinpath(@__DIR__, ".."); io = devnull)                 #jl

using MeanFieldHom
using TensND
using Printf
using Plots
gr()

# ## Setup
#
# A single-layer prolate spheroid (aspect ratio `ω`), concentration
# `c = 0.5`, conductivity contrast `k₁/kₘ = 10³`, matrix `kₘ = 1`
# (isotropic), embedded via Mori-Tanaka. The interface sits at the
# spheroid/matrix boundary (`N = 1` layer ⟹ a single interface).

const km = 1.0
const k1 = 1.0e3
const c_frac = 0.5
const Nseries = 8

_matrix_ellipsoid() = Ellipsoid(1.0, 1.0, 1.0)

function _mt_axial_conductivity(ω::Real, interface; Nseries = Nseries)
    K1 = TensISO{3}(k1)
    Km = TensISO{3}(km)
    axis = ω > 1 ? (1.0, 0.0, 0.0) : (0.0, 0.0, 1.0)   # match `Spheroid(ω)`'s own axis convention
    b = 1.0
    a = b * ω
    s = LayeredSpheroid((a,), (b,), (K1,); interfaces = (interface,), Nseries, axis)
    rve = RVE(:M)
    add_matrix!(rve, _matrix_ellipsoid(), Dict(:K => Km))
    add_phase!(rve, :I, s, Dict(:K => K1); fraction = c_frac)
    Keff = homogenize(rve, MoriTanaka(), :K)
    idx = ω > 1 ? 1 : 3   # axial direction, wherever `axis` put it
    return Keff[idx, idx]
end

# ## Interface parameter sweep
#
# `h̃ = kₘ·b / ρ` (Biot-type number) for the Kapitza (LC) interface —
# `h̃ → ∞` is the perfect-interface limit, `h̃ → 0` the fully insulated
# one. Two aspect ratios and two truncation orders `𝒩` illustrate both
# the physical trend and the series' fast convergence (see script `33`
# for a dedicated convergence study).

const h_tilde = exp10.(range(-2, 3; length = 60))
const omegas = (3.0, 1 / 3)
const Ns_compare = (3, Nseries)

results = Dict{Tuple{Float64, Int}, Vector{Float64}}()
for ω in omegas, N in Ns_compare
    vals = map(h_tilde) do h
        ρ = km / h   # b = 1
        _mt_axial_conductivity(ω, MeanFieldHom.KapitzaInterface(ρ); Nseries = N)
    end
    results[(ω, N)] = vals
end

# ## Exact anchors — perfect interface / fully insulated core
#
# Independent of `h̃`: the two limits the curves must approach at the
# right and left ends of the sweep.

function _mt_perfect(ω)
    K1 = TensISO{3}(k1); Km = TensISO{3}(km)
    axis = ω > 1 ? (1.0, 0.0, 0.0) : (0.0, 0.0, 1.0)
    s = LayeredSpheroid((ω,), (1.0,), (K1,); Nseries, axis)
    rve = RVE(:M)
    add_matrix!(rve, _matrix_ellipsoid(), Dict(:K => Km))
    add_phase!(rve, :I, s, Dict(:K => K1); fraction = c_frac)
    Keff = homogenize(rve, MoriTanaka(), :K)
    return Keff[(ω > 1 ? 1 : 3), (ω > 1 ? 1 : 3)]
end

function _mt_insulated(ω)
    K0 = TensISO{3}(1.0e-9); Km = TensISO{3}(km)
    axis = ω > 1 ? (1.0, 0.0, 0.0) : (0.0, 0.0, 1.0)
    s = LayeredSpheroid((ω,), (1.0,), (K0,); Nseries, axis)
    rve = RVE(:M)
    add_matrix!(rve, _matrix_ellipsoid(), Dict(:K => Km))
    add_phase!(rve, :I, s, Dict(:K => K0); fraction = c_frac)
    Keff = homogenize(rve, MoriTanaka(), :K)
    return Keff[(ω > 1 ? 1 : 3), (ω > 1 ? 1 : 3)]
end

println("Kushch-style imperfect-interface sweep — c = $c_frac, k₁/kₘ = $k1")
println("─"^70)
for ω in omegas
    @printf "ω = %6.3f   k_perfect = %.6f   k_insulated = %.6f\n" ω _mt_perfect(ω) _mt_insulated(ω)
end
println()

# ## Graphical output

p = plot(
    xscale = :log10, xlabel = "h̃ = kₘ·b / ρ", ylabel = "k₁₁ᵉᶠᶠ / kₘ",
    title = "MT effective axial conductivity vs. Kapitza interface parameter\n(c=$c_frac, k₁/kₘ=$k1)",
    legend = :topleft, size = (800, 550),
)
colors = (:steelblue, :darkorange)
for (i, ω) in enumerate(omegas)
    for N in Ns_compare
        plot!(
            p, h_tilde, results[(ω, N)] ./ km;
            label = "ω=$(round(ω; digits = 2)), N=$N", color = colors[i],
            linestyle = N == Ns_compare[1] ? :dash : :solid, lw = 2,
        )
    end
    hline!(p, [_mt_perfect(ω) / km]; color = colors[i], linestyle = :dot, lw = 1, label = nothing)
    hline!(p, [_mt_insulated(ω) / km]; color = colors[i], linestyle = :dot, lw = 1, label = nothing)
end
p

const figdir = joinpath(@__DIR__, "figures")                        #jl
isdir(figdir) || mkdir(figdir)                                       #jl
figpath = joinpath(figdir, "32_spheroid_nlayers_conductivity.png")    #jl
savefig(p, figpath)                                                   #jl
display(p)                                                            #jl
@printf "\nSaved : %s\n" figpath                                      #jl
