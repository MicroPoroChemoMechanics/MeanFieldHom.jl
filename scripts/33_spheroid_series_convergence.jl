# # n-layer confocal spheroid: series truncation & the quadrature vs. BigFloat choice
#
# Barthélémy & Bignonnet (IJES 2020) truncate the spheroidal-harmonic
# series at `𝒩` terms (odd degrees `1, 3, …, 2𝒩-1`) and show that the
# coupling-integral summation `Iᵢⱼ = Σₖ γ₂ₖ Wₖ(q)` needs a working
# precision of `≳ 0.8·(2𝒩-1)` decimal digits — the coefficients `γ₂ₖ`
# grow like `10^{0.8n}` while `Iᵢⱼ` itself is `O(1/n)`, so plain
# `Float64` (16 digits) silently loses accuracy once `𝒩 ≳ 10`. Their
# reference implementation therefore evaluates the coupling matrices
# with `mpmath` at a precision set from `𝒩` (see
# `echoes_cpp/interface/python/py_inclusions/spheroid_nlayers.py`, and
# `spheroid_nlayers_converge_series.py` for the convergence figure
# reproduced below).
#
# `MeanFieldHom.jl` sidesteps the issue rather than reproducing it: the
# default [`coupling_matrices`](@ref MeanFieldHom.LayeredSpheroids.coupling_matrices)
# backend integrates the SAME closed-form integrals directly by Gauss
# quadrature (`QuadGK`), evaluating only Legendre VALUES (never the
# monomial expansion `γ₂ₖ`) at real `x ∈ [-1,1]` — a perfectly smooth,
# well-conditioned integrand for any `𝒩`, so plain `Float64` stays
# accurate well past the point where the original series needs
# multi-precision. The faithful BigFloat port of the original
# monomial-coefficient series (`method = :series`) is also shipped, as
# an independent validation oracle (`test/LayeredSpheroids/test_coupling.jl`
# cross-checks the two to machine precision) and to reproduce this
# exact figure.

import Pkg                                                          #jl
Pkg.activate(joinpath(@__DIR__, ".."); io = devnull)                 #jl

using MeanFieldHom
using MeanFieldHom.LayeredSpheroids: spheroid_ba_ratios, coupling_matrices
using TensND
using Printf
using Plots
gr()

# ## Convergence of `b_{N+1,1}` vs. truncation order 𝒩
#
# Mirrors `spheroid_nlayers_converge_series.py`: a single-layer
# prolate spheroid with an INSULATING core (`k₁ = 0`) and a unit
# surface-conductive (HC) interface, matrix `kₘ = 1`. For each aspect
# ratio `ω`, the relative change of `(b/a)ₐ` and `(b/a)ₜ` between
# `𝒩-1` and `𝒩` is tracked as `𝒩` grows — this is exactly the
# quantity the paper's Fig. `relerror` plots (its `b_{N+1,1}` IS the
# `(b/a)` ratio here since `a_{N+1,1} = 1` throughout).

const Km = TensISO{3}(1.0)
const K0core = TensISO{3}(0.0)
const omegas = (1.1, 2.0, 1.0e1, 1.0e2, 1.0e3, 1.0e4)
const Nmax = 15

function _ba_sequence(ω)
    s(N) = LayeredSpheroid(
        (ω,), (1.0,), (K0core,);
        interfaces = (MeanFieldHom.SurfaceConductiveInterface(1.0),),
        Nseries = N, axis = (1.0, 0.0, 0.0),
    )
    return [spheroid_ba_ratios(s(N), Km) for N in 2:Nmax]
end

rel_err_a = Dict{Float64, Vector{Float64}}()
rel_err_t = Dict{Float64, Vector{Float64}}()
for ω in omegas
    seq = _ba_sequence(ω)
    ba = first.(seq)
    bt = last.(seq)
    rel_err_a[ω] = [abs((ba[i] - ba[i - 1]) / ba[i - 1]) for i in 2:length(ba)]
    rel_err_t[ω] = [abs((bt[i] - bt[i - 1]) / bt[i - 1]) for i in 2:length(bt)]
end

println("Convergence of (b/a)_axial, (b/a)_trans vs 𝒩 (insulating core, unit HC interface)")
println("─"^78)
for ω in omegas
    @printf "ω=%8.1e   rel.diff @ 𝒩=5: axial=%.2e trans=%.2e   @ 𝒩=10: axial=%.2e trans=%.2e\n" ω rel_err_a[ω][4] rel_err_t[ω][4] rel_err_a[ω][9] rel_err_t[ω][9]
end
println()

# ## Quadrature vs. BigFloat series: agreement at the 𝒩 where the
# original algorithm would need extra precision
#
# `𝒩 = 12` needs `0.8×(2×12-1) ≈ 18` decimal digits — already past
# double precision (eq:maxlogcij). The `:quadrature` path (default,
# used above) never forms the ill-conditioned sum; the `:series` path
# (BigFloat, set to the paper's precision rule) is the historical
# approach. Both must still agree, since they compute the same
# integrals two different ways.

q_test = 3.0
for N in (8, 12, 16)
    Iq, Jq, Kq, Lq = coupling_matrices(q_test, N; method = :quadrature)
    Is, Js, Ks, Ls = coupling_matrices(q_test, N; method = :series)
    err = maximum(abs.(Iq .- Is)) / maximum(abs.(Is))
    @printf "𝒩=%2d  (needs %2d digits)  max rel. diff (I matrix, quadrature vs series) = %.2e\n" N round(Int, 0.8 * (2N - 1)) err
end
println()

# ## Graphical output

p1 = plot(
    xlabel = "N (truncation order)", ylabel = "rel. diff. on (b/a)ₐ between N-1 and N",
    yscale = :log10, ylims = (1.0e-16, 1.0e-1), title = "Axial", legend = :topright,
)
p2 = plot(
    xlabel = "N (truncation order)", ylabel = "rel. diff. on (b/a)ₜ between N-1 and N",
    yscale = :log10, ylims = (1.0e-16, 1.0e-1), title = "Transverse", legend = :topright,
)
markers = (:circle, :diamond, :utriangle, :square, :star5, :cross)
for (i, ω) in enumerate(omegas)
    lbl = ω ≥ 10 ? @sprintf("ω=%.0e", ω) : @sprintf("ω=%.1f", ω)
    ns = 3:Nmax
    scatter!(p1, ns, max.(rel_err_a[ω], 1.0e-17); label = lbl, marker = markers[i], ms = 4)
    scatter!(p2, ns, max.(rel_err_t[ω], 1.0e-17); label = lbl, marker = markers[i], ms = 4)
end
p_full = plot(p1, p2; layout = (1, 2), size = (1200, 500))
p_full

const figdir = joinpath(@__DIR__, "figures")                        #jl
isdir(figdir) || mkdir(figdir)                                       #jl
figpath = joinpath(figdir, "33_spheroid_series_convergence.png")      #jl
savefig(p_full, figpath)                                              #jl
display(p_full)                                                       #jl
@printf "\nSaved : %s\n" figpath                                      #jl
