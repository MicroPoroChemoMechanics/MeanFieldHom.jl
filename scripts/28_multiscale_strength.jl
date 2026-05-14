# =============================================================================
#  28_multiscale_strength.jl
#
#  Three-scale upscaling of cement-paste / mortar elasticity and quasi-brittle
#  strength using the multi-orientation hydrate-foam model of
#
#    Pichler, B. and Hellmich, C. (2011), "Upscaling quasi-brittle strength of
#    cement paste and mortar : a multi-scale engineering mechanics model",
#    Cement and Concrete Research 41, 467-476.
#    https://doi.org/10.1016/j.cemconres.2011.01.010
#
#  Scales :
#    1. Hydrate Foam (HF)   : Self-Consistent — flat oblate hydrates
#                              spread over NTHETA orientations + water + air.
#                              Each hydrate family is given a TI symmetrize
#                              around the global axis ez (uniform azimuthal
#                              distribution at fixed polar angle θ_i), so the
#                              homogenised C_hf is isotropic when all families
#                              share the same shear modulus.
#    2. Cement Paste (CP)   : Mori-Tanaka — HF + clinker.
#    3. Mortar       (MO)   : Mori-Tanaka — CP + sand.
#
#  Strength criterion (Pichler & Hellmich 2011) :
#       M  = S_mo : dC_mo/dμ_at_θ : S_mo            (compliance pull-back)
#       fc = 1 / √( M[2,2,2,2] · 2 μ² / f_θ )
#  where μ is the hydrate shear modulus, f_θ is the volume fraction of the
#  perturbed hydrate family in the mortar, and dC_mo / dμ_at_θ is the partial
#  derivative of the homogenised stiffness with respect to the shear modulus
#  of one specific orientation family (others held fixed).
#
#  The sensitivity API of MeanFieldHom v0.4.0 lets us compute that partial in
#  one shot through ForwardDiff : a single Dual variable on the shear modulus
#  of family θ propagates through all three scales (SC + MT + MT), and the
#  TI symmetrize on each family keeps the SC iteration well-typed. No manual
#  chain rule of partial Jacobians is required.
# =============================================================================

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."); io = devnull)

using MeanFieldHom
using ForwardDiff
using TensND
using Printf
using LinearAlgebra
using Plots

# ── Physical constants ─────────────────────────────────────────────────────
const ρ_w = 1.0
const ρ_clin = 3.15;  const d_clin = ρ_clin / ρ_w
const ρ_hyd = 2.073; const d_hyd = ρ_hyd / ρ_w
const ρ_san = 2.648; const d_san = ρ_san / ρ_w

const K_clin, μ_clin = 116.7, 53.8
const K_hyd_ref, μ_hyd_ref = 18.7, 11.8
const K_san, μ_san = 37.8, 44.3
const TINY = 1.0e-3   # numerical regularisation for water and air (≈ zero stiffness)

# Hydrate spheroid aspect ratio and angular discretisation (Pichler 2011)
const NTHETA = 20
const ω_aspect = 1.0e4

# ── Powers volume fractions ────────────────────────────────────────────────
f_clin(wc, α) = (1 - α) / (1 + d_clin * wc)
f_w(wc, α) = d_clin * (wc - 0.42 * α) / (1 + d_clin * wc)
f_hyd(wc, α) = 1.42 * d_clin / d_hyd * α / (1 + d_clin * wc)
fh_san(wc, sc) = sc / d_san / (1 / d_clin + wc + sc / d_san)
αmax(wc) = min(1.0, wc / 0.42)

# ── Angular discretisation : N polar bins on (0, π/2) ──────────────────────
function disc_theta(N::Int)
    thm = vcat(0.0, [π / 2 * (i - 0.5) / (N - 1) for i in 1:(N - 1)])
    theta = [π / 2 * (i - 1) / (N - 1) for i in 1:N]
    thp = vcat([π / 2 * (i - 0.5) / (N - 1) for i in 1:(N - 1)], π / 2)
    return thm, theta, thp
end

# ── Hydrate foam : SC of multi-orientation hydrates + water + air ──────────
#
# Each hydrate orientation family is registered as a separate phase with
# `symmetrize = TISymmetrize((0,0,1))` : the localisation tensor for a
# fixed polar angle θ_i is averaged over the azimuthal angle around ez,
# producing a TI(ez) contribution. The resulting C_hf at a uniform shear
# modulus μ_ref is rotationally isotropic (in the limit of a fine bin grid).
#
# `μ_at_θ` is the shear modulus of the family at the chosen perturbation
# index `θ_idx`, which can be Dual; the others stay at `μ_ref`.
function build_hf(wc, α_p, μ_hyd; ω::Real = ω_aspect)
    fclin = f_clin(wc, α_p)
    fw = f_w(wc, α_p)
    fhyd = f_hyd(wc, α_p)
    fair = max(0.0, 1 - fclin - fw - fhyd)
    fthyd_t = fhyd / (1 - fclin)
    ftw_t = fw / (1 - fclin)
    ftair_t = fair / (1 - fclin)

    T = typeof(μ_hyd)
    rve = RVE(:M; T = T)
    add_matrix!(
        rve, Ellipsoid(1.0),
        Dict(
            :C => TensISO{3}(
                convert(T, 3 * K_hyd_ref),
                convert(T, 2 * μ_hyd_ref)
            )
        )
    )
    # Single hydrate phase as a flat oblate spheroid with iso-symmetrize :
    # the localization tensor is averaged over the uniform spatial
    # distribution of orientations, producing an isotropic homogenised
    # tensor.
    # Prolate hydrate spheroid (axes (1, 1, ω) with ω >> 1 = needle-shape).
    # Iso symmetrize gives a uniform spatial distribution of orientations.
    geom_hyd = Spheroid(ω)
    add_phase!(
        rve, :HYD, geom_hyd,
        Dict(:C => TensISO{3}(convert(T, 3 * K_hyd_ref), 2 * μ_hyd));
        fraction = fthyd_t, symmetrize = :iso
    )
    add_phase!(
        rve, :W, Ellipsoid(1.0),
        Dict(:C => TensISO{3}(convert(T, 3 * TINY), convert(T, 2 * TINY)));
        fraction = ftw_t
    )
    add_phase!(
        rve, :AIR, Ellipsoid(1.0),
        Dict(:C => TensISO{3}(convert(T, 3 * TINY), convert(T, 2 * TINY)));
        fraction = ftair_t
    )
    return homogenize(
        rve, SelfConsistent(;
            abstol = 1.0e-8, maxiters = 1000,
            damping = 0.5
        ),
        :C; select_best = true
    )
end

# ── Cement paste : MT(HF, clinker) ─────────────────────────────────────────
function build_cp(wc, α_p, C_hf::TensND.AbstractTens)
    fclin = f_clin(wc, α_p)
    T = eltype(C_hf)
    rve = RVE(:HF; T = T)
    add_matrix!(rve, Ellipsoid(1.0), Dict(:C => C_hf))
    add_phase!(
        rve, :CLIN, Ellipsoid(1.0),
        Dict(:C => TensISO{3}(convert(T, 3 * K_clin), convert(T, 2 * μ_clin)));
        fraction = fclin
    )
    return homogenize(rve, MoriTanaka(), :C)
end

# ── Mortar : MT(CP, sand) ──────────────────────────────────────────────────
function build_mo(wc, sc, C_cp::TensND.AbstractTens)
    fsan = fh_san(wc, sc)
    T = eltype(C_cp)
    rve = RVE(:CP; T = T)
    add_matrix!(rve, Ellipsoid(1.0), Dict(:C => C_cp))
    add_phase!(
        rve, :SAN, Ellipsoid(1.0),
        Dict(:C => TensISO{3}(convert(T, 3 * K_san), convert(T, 2 * μ_san)));
        fraction = fsan
    )
    return homogenize(rve, MoriTanaka(), :C)
end

# ── Full multi-scale chain : input scalar μ_hyd, output stiffness array ──
function multiscale_C_mo(wc, α_p, sc, μ_hyd)
    C_hf = build_hf(wc, α_p, μ_hyd)
    C_cp = build_cp(wc, α_p, C_hf)
    C_mo = build_mo(wc, sc, C_cp)
    return get_array(C_mo)
end

# ── Iso bulk and shear moduli of a (nearly) iso 4-tensor ───────────────────
function extract_kμ(arr::AbstractArray)
    K = sum(arr[i, i, j, j] for i in 1:3, j in 1:3) / 9
    full_trace = sum(arr[i, j, i, j] for i in 1:3, j in 1:3)
    μ = (full_trace - 3K) / 10
    return K, μ
end
extract_E(K, μ) = 9K * μ / (3K + μ)

# ── Strength criterion (Pichler & Hellmich 2011) ───────────────────────────
#
# At a value point where C_mo is iso, S_mo is also iso. M = S_mo : dC : S_mo
# is computed by tensor double-dot :  M_ijkl = S_ijab · dC_abcd · S_cdkl.
# The Pichler component is M[2,2,2,2] (axial-axial in compliance pull-back).
function pichler_strength(
        arr_C_mo::AbstractArray,
        arr_dC::AbstractArray,
        μh::Real,
        f_θ::Real
    )
    K_mo, μ_mo = extract_kμ(arr_C_mo)
    # Iso compliance S_mo = (1/9K) δ_ij δ_kl + (1/4μ)(δ_ik δ_jl + δ_il δ_jk - 2/3 δ_ij δ_kl)
    Sα = 1 / (3 * 3 * K_mo)        # spherical projector coefficient
    Sβ = 1 / (2 * 2 * μ_mo)        # deviatoric projector coefficient
    arr_S = zeros(eltype(arr_dC), 3, 3, 3, 3)
    for i in 1:3, j in 1:3, k in 1:3, l in 1:3
        sph = Sα * (i == j) * (k == l)
        dev = Sβ * (
            ((i == k) * (j == l) + (i == l) * (j == k)) -
                2 * (i == j) * (k == l) / 3
        )
        # Note : S_ijkl = (1/9K) δ_ij δ_kl + (1/4μ) (...)
        arr_S[i, j, k, l] = (i == j) * (k == l) / (9 * K_mo) +
            (
            ((i == k) * (j == l) + (i == l) * (j == k)) -
                2 * (i == j) * (k == l) / 3
        ) / (4 * μ_mo)
    end
    M = zeros(eltype(arr_dC), 3, 3, 3, 3)
    @inbounds for i in 1:3, j in 1:3, k in 1:3, l in 1:3
        s = zero(eltype(arr_dC))
        for a in 1:3, b in 1:3, c in 1:3, d in 1:3
            s += arr_S[i, j, a, b] * arr_dC[a, b, c, d] * arr_S[c, d, k, l]
        end
        M[i, j, k, l] = s
    end
    # Pichler & Hellmich 2011 use the axial-axial component of M (along the
    # symmetry axis ez = e_3) in their strength criterion. In tensor index
    # notation this is M_3333.
    M_axial = M[3, 3, 3, 3]
    return 1 / sqrt(abs(M_axial) * 2 * μh^2 / f_θ)
end

# ── IFT-based partial derivative dC_hf / dμ_HYD₀_TI[bin0] ───────────────────
#
# The Pichler-Hellmich strength criterion uses the partial of the
# mortar stiffness with respect to the shear modulus of *one*
# orientation bin of the hydrate phase (bin 0, axis ‖ ez), pulled
# through the three-scale chain. We compute this partial via the
# implicit-function theorem on a single multi-bin TI SC step taken
# from the iso fixed point `C_hf` (computed by the single-iso SC
# above): at the SC fixed point `C_hf*`, the analytical IFT formula
# `dC_hf*/dp = (I − ∂F/∂C)⁻¹ ∂F/∂p` is exact. Both Jacobians are
# obtained by ForwardDiff on the multi-bin SC step `F`; a small linear
# solve gives the partial in TI(ez) Walpole basis. The chain through
# CP and MO is then propagated by another ForwardDiff pass on a
# scalar perturbation along that direction.

function _disc_theta(N)
    thm = vcat(0.0, [π / 2 * (i - 0.5) / (N - 1) for i in 1:(N - 1)])
    theta = [π / 2 * (i - 1) / (N - 1) for i in 1:N]
    thp = vcat([π / 2 * (i - 0.5) / (N - 1) for i in 1:(N - 1)], π / 2)
    return thm, theta, thp
end

_iso_to_walpole(α, β) =
    [(α + 2β) / 3, (2α + β) / 3, sqrt(2.0) * (α - β) / 3, β, β]

# ── Mandel 6×6 helpers ───────────────────────────────────────────────────────
# Manual Mandel encoding for 4-tensors with minor symmetry (but possibly NOT
# major symmetry).  Bypasses TensND types so that summing per-bin TI(rotated)
# contributions stays in 6×6 matrix arithmetic — Tensors.tomandel only handles
# SymmetricTensor (sym 6×6) or general Tensor (rank-deficient 9×9).

const _MANDEL_IDX = ((1, 1), (2, 2), (3, 3), (2, 3), (1, 3), (1, 2))
@inline _mandel_scl(::Type{T}, k::Int) where {T} =
    k <= 3 ? one(T) : sqrt(T(2))

function _arr_to_mandel66(arr::AbstractArray{T, 4}) where {T}
    M = zeros(T, 6, 6)
    @inbounds for a in 1:6, b in 1:6
        i, j = _MANDEL_IDX[a]
        k, l = _MANDEL_IDX[b]
        M[a, b] = _mandel_scl(T, a) * _mandel_scl(T, b) * arr[i, j, k, l]
    end
    return M
end

function _mandel66_to_arr(M::AbstractMatrix{T}) where {T}
    arr = zeros(T, 3, 3, 3, 3)
    @inbounds for a in 1:6, b in 1:6
        i, j = _MANDEL_IDX[a]
        k, l = _MANDEL_IDX[b]
        v = M[a, b] / (_mandel_scl(T, a) * _mandel_scl(T, b))
        arr[i, j, k, l] = v
        arr[j, i, k, l] = v
        arr[i, j, l, k] = v
        arr[j, i, l, k] = v
    end
    return arr
end

# Iso projection of a TensTI{4,T,5}(ez) — analytical formula on Walpole 5-vec.
function _walpole_ez_to_iso(c_data::AbstractVector{T}) where {T}
    ℓ₁, ℓ₂, ℓ₃, ℓ₅, ℓ₆ = c_data[1], c_data[2], c_data[3], c_data[4], c_data[5]
    sq2 = sqrt(T(2))
    α_iso = (ℓ₁ + 2 * ℓ₂ + 2 * sq2 * ℓ₃) / 3
    full_trace = ℓ₁ + 2 * ℓ₂ + 2 * ℓ₅ + 2 * ℓ₆
    β_iso = (full_trace - α_iso) / 5
    return α_iso, β_iso
end

# Project a 4-array onto TI(ez) Walpole 5-vec (major-symmetric component).
function _arr_to_walpole_ez(arr::AbstractArray{T, 4}) where {T}
    sq2 = sqrt(T(2))
    ℓ₁ = arr[3, 3, 3, 3]
    ℓ₂ = (arr[1, 1, 1, 1] + arr[1, 1, 2, 2] + arr[2, 2, 1, 1] + arr[2, 2, 2, 2]) / 2
    ℓ₃a = (arr[3, 3, 1, 1] + arr[3, 3, 2, 2]) / sq2
    ℓ₄a = (arr[1, 1, 3, 3] + arr[2, 2, 3, 3]) / sq2
    ℓ₃ = (ℓ₃a + ℓ₄a) / 2
    ℓ₅ = 2 * arr[1, 2, 1, 2]
    ℓ₆ = 2 * arr[1, 3, 1, 3]
    return [ℓ₁, ℓ₂, ℓ₃, ℓ₅, ℓ₆]
end

# ── Custom multi-bin TI SC step in 6×6 Mandel ─────────────────────────────────
# Bypasses MFH's `_sc_step` (which can't handle per-bin TI(rotated_axis)
# contributions because the resulting `+(TensTI{4,5}(axis_i), TensTI{4,5}(axis_j))`
# triggers an axis-mismatch assertion in TensND, and the unstructured fallback
# hits a rank-deficient 9×9 LU on inversion).  Each phase's localisation is
# computed via the analytical TI-coaxial Hill tensor (Dual-friendly), converted
# to a 6×6 Mandel matrix, and the SC step is closed in 6×6 matrix arithmetic.
# Mirrors the C++ reference's per-bin `set_reference("C", X2_TI(theta_i))` line
# in `homogenization_scheme.h::evaluate`.
function _F_walpole(c_data::AbstractVector, μ_b0, wc, α_p; N = NTHETA)
    T = promote_type(eltype(c_data), typeof(μ_b0))

    cd = ntuple(i -> convert(T, c_data[i]), 5)
    C_0_TI = TensND.TensTI{4, T, 5}(
        cd,
        (convert(T, 0.0), convert(T, 0.0), convert(T, 1.0))
    )
    α_iso0, β_iso0 = _walpole_ez_to_iso(collect(cd))
    C_0_iso = TensND.TensISO{3}(α_iso0, β_iso0)

    fclin = f_clin(wc, α_p); fw = f_w(wc, α_p); fhyd = f_hyd(wc, α_p)
    fair = max(0.0, 1 - fclin - fw - fhyd)
    fthyd_t = fhyd / (1 - fclin); ftw_t = fw / (1 - fclin); ftair_t = fair / (1 - fclin)
    thm, thetas, thp = _disc_theta(N)

    A_avg = zeros(T, 6, 6)
    CA_avg = zeros(T, 6, 6)

    # HYD bins — per-bin TI(spheroid_axis) projection of the reference, analytical
    # TI-coaxial Hill (since projected matrix and spheroid share the bin axis).
    for i in 1:N
        f_i = fthyd_t * (cos(thm[i]) - cos(thp[i]))
        θ = thetas[i]
        bin_axis = (convert(T, sin(θ)), zero(T), convert(T, cos(θ)))
        spheroid = Spheroid(ω_aspect; euler_angles = (θ, 0.0, 0.0))

        C_ref_i = MeanFieldHom.Schemes._apply_symmetrize(
            C_0_TI, MeanFieldHom.Schemes.TISymmetrize(bin_axis)
        )

        μ_h = (i == 1) ? μ_b0 : convert(T, μ_hyd_ref)
        C_HYD_i = TensND.TensISO{3}(convert(T, 3 * K_hyd_ref), 2 * μ_h)

        A_dil_i = strain_strain_loc(spheroid, C_HYD_i, C_ref_i)
        A_arr = TensND.get_array(A_dil_i)
        CA_i = C_HYD_i ⊡ A_dil_i
        CA_arr = TensND.get_array(CA_i)

        A_avg .+= f_i .* _arr_to_mandel66(A_arr)
        CA_avg .+= f_i .* _arr_to_mandel66(CA_arr)
    end

    # Spherical phases (W, AIR) — iso projection of the reference.
    sphere = Ellipsoid(1.0)
    C_W = TensND.TensISO{3}(convert(T, 3 * TINY), convert(T, 2 * TINY))
    A_W = strain_strain_loc(sphere, C_W, C_0_iso)
    A_W_arr = TensND.get_array(A_W)
    CA_W = C_W ⊡ A_W
    CA_W_arr = TensND.get_array(CA_W)
    A_W_KM = _arr_to_mandel66(A_W_arr)
    CA_W_KM = _arr_to_mandel66(CA_W_arr)
    A_avg .+= ftw_t .* A_W_KM
    CA_avg .+= ftw_t .* CA_W_KM
    A_avg .+= ftair_t .* A_W_KM
    CA_avg .+= ftair_t .* CA_W_KM

    F_KM = CA_avg / A_avg
    F_arr = _mandel66_to_arr(F_KM)
    return _arr_to_walpole_ez(F_arr)
end

function _dCh_dμb0_walpole(α_hf, β_hf, wc, α_p)
    c0 = _iso_to_walpole(α_hf, β_hf)
    # Linearise at one SC step from the iso `C_hf` (mirroring the C++
    # reference's `homogenize(... maxnb=1)` before `homogenize_derivative`).
    c1 = _F_walpole(c0, μ_hyd_ref, wc, α_p)
    J_c = ForwardDiff.jacobian(c -> _F_walpole(c, μ_hyd_ref, wc, α_p), c1)
    df_dμ = ForwardDiff.derivative(μ -> _F_walpole(c1, μ, wc, α_p), μ_hyd_ref)
    # Divide by 2 to convert from `d/dμ` (our parameterisation) to
    # `d/d(2μ)` (the iso parameter index used in the reference's strength
    # criterion `M[3,3,3,3] · 2 μ²`).
    return ((LinearAlgebra.I - J_c) \ df_dμ) ./ 2
end

# ── Compute one (wc, α) point : moduli + strength criterion ───────────────
function compute_point(wc, α_p; sc = 0.0)
    # Value
    arr_C_mo = multiscale_C_mo(wc, α_p, sc, μ_hyd_ref)
    K_mo, μ_mo = extract_kμ(arr_C_mo)
    E_mo = extract_E(K_mo, μ_mo)

    # Bin-0 fraction of the hydrate phase in the mortar
    fhyd_in_mortar = f_hyd(wc, α_p) * (1 - fh_san(wc, sc))
    w_bin0 = 1.0 - cos((π / 2) * 0.5 / (NTHETA - 1))
    f_θ = fhyd_in_mortar * w_bin0

    # IFT partial in TI(ez) Walpole basis at the iso C_hf fixed point.
    C_hf = build_hf(wc, α_p, μ_hyd_ref)
    α_hf, β_hf = TensND.get_data(C_hf)
    dCh_dμ = _dCh_dμb0_walpole(α_hf, β_hf, wc, α_p)

    # Push the partial through CP (MT) and MO (MT) via one scalar
    # ForwardDiff pass along the direction `dCh_dμ`.
    f_chain = function (t)
        Cw = _iso_to_walpole(α_hf, β_hf) .+ t .* dCh_dμ
        Tt = typeof(t)
        Cwt = ntuple(i -> convert(Tt, Cw[i]), 5)
        C_hf_perturbed = TensND.TensTI{4, Tt, 5}(
            Cwt,
            (zero(Tt), zero(Tt), one(Tt))
        )
        C_cp = build_cp(wc, α_p, C_hf_perturbed)
        C_mo = build_mo(wc, sc, C_cp)
        return get_array(C_mo)
    end
    arr_dC_mo = ForwardDiff.derivative(f_chain, 0.0)

    fc = pichler_strength(arr_C_mo, arr_dC_mo, μ_hyd_ref, f_θ)
    return (; K_mo, μ_mo, E_mo, fc)
end

# ── Sweep + plot (figure structure of Pichler & Hellmich 2011, Fig. 4) ────
println("="^78)
println("Multi-scale upscaling of cement-paste / mortar (Pichler-Hellmich 2011)")
println("(NTHETA = $NTHETA, ω = $ω_aspect)")
println("="^78)

const wcs = [0.157, 0.25, 0.35, 0.5, 0.65, 0.8]
const sc_default = 0.0
const N_α = 20
const α_min = 0.005   # minimum hydration degree to plot (curves go to zero as α→0)

p1 = plot(;
    xlabel = "α", ylabel = "k_mortar (GPa)",
    xlims = (0, 1), ylims = (0, 35), legend = :topleft
)
p2 = plot(;
    xlabel = "α", ylabel = "μ_mortar (GPa)",
    xlims = (0, 1), ylims = (0, 20), legend = false
)
p3 = plot(;
    xlabel = "α", ylabel = "f_c / σ_ult",
    xlims = (0, 1), ylims = (0, 2), legend = false
)
p4 = plot(;
    xlabel = "f_c / σ_ult", ylabel = "E_mortar (GPa)",
    xlims = (0, 2), ylims = (0, 50), legend = false
)

for wc in wcs
    # Clamp αmax just below 1 to avoid the matrix-fraction = -ε rounding
    # warning at αmax exactly (sum of inclusion fractions = 1 + machine eps).
    αs = collect(range(α_min, αmax(wc) * (1 - 1.0e-12); length = N_α))
    K_arr = Float64[]; μ_arr = Float64[]; E_arr = Float64[]; fc_arr = Float64[]
    α_kept = Float64[]
    for α_p in αs
        try
            r = compute_point(wc, α_p; sc = sc_default)
            push!(K_arr, r.K_mo); push!(μ_arr, r.μ_mo); push!(E_arr, r.E_mo)
            push!(fc_arr, r.fc); push!(α_kept, α_p)
        catch e
            @warn "compute_point failed at wc=$wc α=$α_p" exception = (e, catch_backtrace())
        end
    end
    if !isempty(K_arr)
        plot!(p1, α_kept, K_arr, lw = 2, label = "wc = $wc")
        plot!(p2, α_kept, μ_arr, lw = 2)
        plot!(p3, α_kept, fc_arr, lw = 2)
        plot!(p4, fc_arr, E_arr, lw = 2)
    end
end

p_full = plot(
    p1, p2, p3, p4; layout = (2, 2), size = (1000, 800),
    plot_title = "Multi-scale strength upscaling — MeanFieldHom v0.4"
)

figdir = joinpath(@__DIR__, "figures")
isdir(figdir) || mkdir(figdir)
figpath = joinpath(figdir, "28_multiscale_strength.png")
savefig(p_full, figpath)
@printf "\nSaved : %s\n" figpath

# Tabular summary at one wc
println("\n[Tabular] wc = 0.50")
@printf "  %5s   %10s   %10s   %10s   %10s\n" "α" "k_mo" "μ_mo" "fc" "E_mo"
println("  " * "─"^60)
let wc = 0.5
    αs = collect(filter(α -> α > 0.05, range(0.05, αmax(wc); length = 6)))
    for α_p in αs
        try
            r = compute_point(wc, α_p)
            @printf "  %.3f   %10.4f   %10.4f   %10.4f   %10.4f\n" α_p r.K_mo r.μ_mo r.fc r.E_mo
        catch e
            @printf "  %.3f   (failed : %s)\n" α_p typeof(e)
        end
    end
end

println("\nDone.")
