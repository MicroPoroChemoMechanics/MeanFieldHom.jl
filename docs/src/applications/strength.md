# Quasi-brittle strength of cement paste and mortar

This chapter implements the multi-scale strength-upscaling model of
[pichler2011](@cite) — and mirrors the corresponding chapter of the Echoes book
[echoes](@cite). It predicts the macroscopic uniaxial **compressive strength**
``f_c`` of a mortar from the strength of a single hydrate phase, by upscaling
through three nested homogenization scales.

| Scale | RVE | Phases | Scheme |
|:-----:|:----|:-------|:------:|
| 1 | Hydrate foam (HF) | oriented hydrate needles + water + air | Self-Consistent |
| 2 | Cement paste (CP) | HF matrix + clinker grains | Mori-Tanaka |
| 3 | Mortar (MO) | CP matrix + sand grains | Mori-Tanaka |

!!! note "One autodiff pass instead of a hand-rolled chain rule"
    The strength criterion needs ``\partial \mathbb C^{\rm hom}_{MO}/\partial
    \mu_{\rm hyd}`` — the sensitivity of the mortar stiffness to the hydrate shear
    modulus, propagated through all three scales. In Echoes this is assembled by
    an explicit chain rule over the five transversely-isotropic parameters of each
    intermediate tensor (`homogenize_derivative` per index, per scale). In
    `MeanFieldHom` the same quantity is a **single [ForwardDiff](https://github.com/JuliaDiff/ForwardDiff.jl)
    pass** through the entire nested chain (multi-bin SC + two MT stages): making
    the θ = 0 family's shear modulus a `Dual` and reading the partial of the final
    `C_mo` gives value and derivative at once. The two routes agree to the
    validation tolerances of `scripts/bench_echoes/benchmark_pichler.jl` (moduli
    within 1 %, ``f_c`` within 2 %).

The code below is the model of
[`scripts/common/pichler_model.jl`](https://github.com/MicroPoroChemoMechanics/MeanFieldHom.jl/blob/main/scripts/common/pichler_model.jl),
reproduced inline.

## Constants and volume fractions

Densities (relative to water) and phase moduli are those of
[pichler2011](@cite). Powers' hydration model gives the volume fractions of
clinker, water, hydrates and air as functions of ``w/c`` and the hydration
degree ``\alpha`` (``\alpha_{\max} = \min(1, w/c / 0.42)``); the sand fraction
follows from the sand-to-cement ratio ``s/c``.

```@example strength
using MeanFieldHom
using TensND
using ForwardDiff
using LinearAlgebra
using Plots
gr()  # headless backend; GKSwstype is set to "100" in make.jl

const d_clin, d_hyd, d_san = 3.15, 2.073, 2.648
const K_clin, μ_clin = 116.7, 53.8
const K_hyd_ref, μ_hyd_ref = 18.7, 11.8
const K_san, μ_san = 37.8, 44.3

# Water/air are regularized to a tiny positive stiffness: it selects the
# percolating SC branch and keeps the iteration smooth (deliberate deviation
# from the exact zero — the source of the ~2% f_c gap with Echoes).
const TINY = 1.0e-3
const NTHETA = 20
const ω_aspect = 1.0e4

f_clin(wc, α) = (1 - α) / (1 + d_clin * wc)
f_w(wc, α) = d_clin * (wc - 0.42α) / (1 + d_clin * wc)
f_hyd(wc, α) = 1.42 * d_clin / d_hyd * α / (1 + d_clin * wc)
fh_san(wc, sc) = sc / d_san / (1 / d_clin + wc + sc / d_san)
αmax(wc) = min(1.0, wc / 0.42)
nothing # hide
```

## Scale 1 — hydrate foam (multi-bin self-consistent)

The hydrate crystals are near-perfect **needles** (prolate spheroids,
``\omega = 10^4``) with a random orientation distribution. The distribution is
discretized into `NTHETA` polar bins (`polar_orientation_bins`); each bin is one
`Spheroid` family tilted by its polar angle and carrying a
`TISymmetrize` — the *exact* azimuthal average about the global axis
``e_z``, preserving the non-major-symmetric content of the concentration tensors.
The foam (needles + water + air) is homogenized self-consistently. The shear
modulus of the **first** family (``\theta = 0``) is kept as a separate argument
`μ_b0` so it can later be seeded with a `ForwardDiff.Dual`.

```@example strength
function build_hf(wc, α_p, μ_b0; N = NTHETA, ω = ω_aspect)
    fclin, fw, fhyd = f_clin(wc, α_p), f_w(wc, α_p), f_hyd(wc, α_p)
    fair = max(0.0, 1 - fclin - fw - fhyd)
    fthyd, ftw, ftair = fhyd / (1 - fclin), fw / (1 - fclin), fair / (1 - fclin)

    T = typeof(μ_b0)
    ez = (0.0, 0.0, 1.0)
    rve = RVE(:M; T = T)
    # Zero-volume matrix = SC seed only (Σ inclusion fractions = 1).
    add_matrix!(
        rve, Ellipsoid(1.0),
        Dict(:C => TensISO{3}(convert(T, 3K_hyd_ref), convert(T, 2μ_hyd_ref)));
        symmetrize = :iso
    )
    for (i, bin) in enumerate(polar_orientation_bins(N))
        μ_h = i == 1 ? μ_b0 : convert(T, μ_hyd_ref)
        add_phase!(
            rve, Symbol(:HYD, i),
            Spheroid(ω; euler_angles = (bin.θ, 0.0, 0.0)),
            Dict(:C => TensISO{3}(convert(T, 3K_hyd_ref), 2μ_h));
            fraction = fthyd * bin.weight, symmetrize = TISymmetrize(ez)
        )
    end
    C_tiny = TensISO{3}(convert(T, 3TINY), convert(T, 2TINY))
    add_phase!(rve, :W, Ellipsoid(1.0), Dict(:C => C_tiny); fraction = ftw, symmetrize = :iso)
    add_phase!(rve, :AIR, Ellipsoid(1.0), Dict(:C => C_tiny); fraction = ftair, symmetrize = :iso)
    return homogenize(
        rve, SelfConsistent(; abstol = 1.0e-8, maxiters = 1000, damping = 0.5),
        :C; select_best = true
    )
end
nothing # hide
```

## Scales 2 and 3 — cement paste and mortar (Mori-Tanaka)

The converged foam is a `TensTI{4}`; `best_fit_ti` extracts its physical
5-parameter TI stiffness (Echoes does the same with `tensor(C_hf.array, TI)`)
before the two Mori-Tanaka stages: clinker grains in the foam matrix (CP), then
sand grains in the paste matrix (MO).

```@example strength
function build_cp(wc, α_p, C_hf)
    T = eltype(C_hf)
    rve = RVE(:HF; T = T)
    add_matrix!(rve, Ellipsoid(1.0), Dict(:C => C_hf))
    add_phase!(
        rve, :CLIN, Ellipsoid(1.0),
        Dict(:C => TensISO{3}(convert(T, 3K_clin), convert(T, 2μ_clin)));
        fraction = f_clin(wc, α_p)
    )
    return homogenize(rve, MoriTanaka(), :C)
end

function build_mo(wc, sc, C_cp)
    T = eltype(C_cp)
    rve = RVE(:CP; T = T)
    add_matrix!(rve, Ellipsoid(1.0), Dict(:C => C_cp))
    add_phase!(
        rve, :SAN, Ellipsoid(1.0),
        Dict(:C => TensISO{3}(convert(T, 3K_san), convert(T, 2μ_san)));
        fraction = fh_san(wc, sc)
    )
    return homogenize(rve, MoriTanaka(), :C)
end

# Full chain: scalar μ_b0 in, C_mo array out.
function multiscale_C_mo(wc, α_p, sc, μ_b0; N = NTHETA, ω = ω_aspect)
    C_hf = best_fit_ti(build_hf(wc, α_p, μ_b0; N = N, ω = ω), (0.0, 0.0, 1.0))
    C_cp = build_cp(wc, α_p, C_hf)
    C_mo = build_mo(wc, sc, C_cp)
    return get_array(C_mo)
end
nothing # hide
```

## Strength criterion and the autodiff sensitivity

The compliance pull-back ``\mathbf M = \mathbf S_{MO}:\partial\mathbb
C_{MO}/\partial\mu_{\rm hyd}:\mathbf S_{MO}`` gives the axial term ``M_{3333}``,
and [pichler2011](@cite)'s criterion reads

```math
\frac{f_c}{\sigma^{\rm ult}_{\rm hyd}} = \frac{1}{\sqrt{\,M_{3333}\,2\mu_{\rm hyd}^2 / f_\theta\,}},
```

where ``f_\theta`` is the mortar volume fraction of the perturbed (``\theta=0``)
hydrate family. The sensitivity ``\partial\mathbb C_{MO}/\partial\mu_{\rm hyd}``
is obtained by seeding ``\mu_{b0}`` with a `Dual` and reading the partial of the
final array — one pass through the whole three-scale chain.

```@example strength
extract_kμ(arr) = k_mu(TensND.proj_tens(Val(:ISO), arr)[1])

function pichler_strength(arr_C_mo, arr_dC, μh, f_θ)
    K_mo, μ_mo = extract_kμ(arr_C_mo)
    S = zeros(eltype(arr_dC), 3, 3, 3, 3)
    for i in 1:3, j in 1:3, k in 1:3, l in 1:3
        S[i, j, k, l] = (i == j) * (k == l) / (9K_mo) +
            (((i == k) * (j == l) + (i == l) * (j == k)) - 2 * (i == j) * (k == l) / 3) / (4μ_mo)
    end
    M = zero(eltype(arr_dC))
    @inbounds for a in 1:3, b in 1:3, c in 1:3, d in 1:3
        M += S[3, 3, a, b] * arr_dC[a, b, c, d] * S[c, d, 3, 3]
    end
    return 1 / sqrt(abs(M) * 2μh^2 / f_θ)
end

function compute_point(wc, α_p; sc = 0.0, N = NTHETA, ω = ω_aspect)
    TagT = typeof(ForwardDiff.Tag(multiscale_C_mo, Float64))
    μ_dual = ForwardDiff.Dual{TagT}(float(μ_hyd_ref), 1.0)
    arr = multiscale_C_mo(wc, α_p, sc, μ_dual; N = N, ω = ω)
    arr_C = ForwardDiff.value.(arr)
    arr_dC = ForwardDiff.partials.(arr, 1)

    K_mo, μ_mo = extract_kμ(arr_C)
    E_mo, _ = E_nu(iso_stiffness(K_mo, μ_mo))
    f_θ = f_hyd(wc, α_p) * (1 - fh_san(wc, sc)) * polar_orientation_bins(N)[1].weight
    # criterion uses the iso parameter 2μ: d/d(2μ) = (1/2) d/dμ
    fc = pichler_strength(arr_C, arr_dC ./ 2, μ_hyd_ref, f_θ)
    return (; K_mo, μ_mo, E_mo, fc)
end

r = compute_point(0.5, αmax(0.5) * (1 - 1e-12))
(k = round(r.K_mo, digits = 3), μ = round(r.μ_mo, digits = 3),
    E = round(r.E_mo, digits = 2), fc = round(r.fc, digits = 4))
```

## Results — strength and stiffness vs hydration degree

For pure cement paste (``s/c = 0``), the effective bulk and shear moduli and the
normalized compressive strength ``f_c/\sigma^{\rm ult}_{\rm hyd}`` all rise
monotonically with hydration and fall with ``w/c`` — the trend reported by
[pichler2011](@cite).

```@example strength
wc_list = (0.25, 0.35, 0.5, 0.65)

pk = plot(; xlabel = "α", ylabel = "k_mortar (GPa)", legend = :topleft, framestyle = :box)
pμ = plot(; xlabel = "α", ylabel = "μ_mortar (GPa)", legend = false, framestyle = :box)
pfc = plot(; xlabel = "α", ylabel = "f_c / σ_ult", legend = false, framestyle = :box)

for wc in wc_list
    αs = range(0.05, αmax(wc) * (1 - 1e-9); length = 14)
    ks, μs, fcs, ak = Float64[], Float64[], Float64[], Float64[]
    for α in αs
        try
            r = compute_point(wc, α)
            push!(ak, α); push!(ks, r.K_mo); push!(μs, r.μ_mo); push!(fcs, r.fc)
        catch
        end
    end
    plot!(pk, ak, ks; label = "w/c=$wc", lw = 2)
    plot!(pμ, ak, μs; lw = 2)
    plot!(pfc, ak, fcs; lw = 2)
end

plot(pk, pμ, pfc; layout = (1, 3), size = (1050, 340),
    plot_title = "Cement paste (s/c = 0)")
```

## Effect of the sand content

Adding sand (increasing ``s/c``) dilutes the load-bearing hydrate foam and lowers
the normalized strength at full hydration, the effect being sharper at low
``w/c``:

```@example strength
psc = plot(; xlabel = "sand-to-cement ratio s/c", ylabel = "f_c / σ_ult",
    legend = :topright, framestyle = :box, size = (720, 440))
for wc in (0.35, 0.5, 0.65)
    scs = range(0.0, 5.0; length = 12)
    vals, kept = Float64[], Float64[]
    for sc in scs
        try
            fc = compute_point(wc, αmax(wc) * (1 - 1e-9); sc = sc).fc
            fc > 0 && (push!(kept, sc); push!(vals, fc))
        catch
        end
    end
    plot!(psc, kept, vals; label = "w/c=$wc", lw = 2, marker = :circle, ms = 3)
end
psc
```

The absolute strength in MPa is recovered by multiplying by
``\sigma^{\rm ult}_{\rm hyd}``, which [pichler2011](@cite) calibrate to ≈ 70–90 MPa
for typical C-S-H.
