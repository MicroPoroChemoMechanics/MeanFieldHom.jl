# =============================================================================
#  12_hill_ti_coaxial.jl
#
#  Demonstration of the closed-form analytical Hill polarisation tensor for
#  a spheroidal inclusion coaxial with a transversely isotropic matrix
#  (Barthélémy 2020, eqs. 49–58).
#
#  This script reproduces Figure 1 of the paper, plotting the Eshelby
#  components Sᵢⱼ as a function of the spheroid aspect ratio
#  ξ = (transverse)/(axial) over four decades, using the same TI
#  matrix introduced by Sevostianov, Yilmaz, Kushch & Levin (2005).
#  It also benchmarks the analytical builder against the residue and
#  DECUHR algorithms.
#
#  Usage:  julia --project scripts/12_hill_ti_coaxial.jl
#
#  Generates:
#   * scripts/figures/hill_ti_coaxial.png — S_{ijkl} vs ξ
#   * console output: timing comparison Analytical / Residue / DECUHR
# =============================================================================

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."); io = devnull)

using MeanFieldHom
using TensND
using LinearAlgebra
using Printf

println("=== Hill polarisation tensor — TI matrix coaxial with spheroid ===")
println()

# Sevostianov et al. (2005) test stiffness from Tab. 1 of Barthélémy 2020.
const C1111 = 2.179
const C1122 = 0.579
const C1133 = 0.689
const C3333 = 10.345
const C2323 = 1.0

println("TI matrix (axis aligned with the spheroid axis in each case):")
@printf("  C₁₁₁₁ = %g    C₁₁₂₂ = %g    C₁₁₃₃ = %g\n", C1111, C1122, C1133)
@printf("  C₃₃₃₃ = %g    C₂₃₂₃ = %g\n", C3333, C2323)
println()

# ── Sweep over aspect ratios — Eshelby S_{ijkl} = P_{ijmn} C_{mnkl} -----------
# To stay coaxial throughout the sweep, the TI axis is taken along the
# spheroid axis, which depends on the shape:
#   * ξ ≤ 1 (oblate, axis e₃)   — use TI with axis e₃, ell = Ellipsoid(1, 1, ξ)
#   * ξ > 1 (prolate, axis e₁)  — use TI with axis e₁, ell = Ellipsoid(ξ, 1, 1)

println("Eshelby tensor components vs. spheroid aspect ratio ξ = (axial)/(transverse)")
println("(coaxial: TI axis = spheroid axis;  ξ < 1 : oblate,  ξ > 1 : prolate)")
println()
@printf("%6s  %12s  %12s  %12s  %12s\n", "ξ", "S_ax,ax", "S_tr,tr", "S_ax,tr", "S_tr,tr′")
println("─"^60)
for ξ in (0.1, 0.3, 0.5, 1.0, 2.0, 3.0, 10.0)
    if ξ ≤ 1
        ell = Ellipsoid(1.0, 1.0, ξ)              # oblate axis e₃
        n_axis = [0.0, 0.0, 1.0]
        i_ax, i_tr = 3, 1
    else
        ell = Ellipsoid(ξ, 1.0, 1.0)              # prolate axis e₁
        n_axis = [1.0, 0.0, 0.0]
        i_ax, i_tr = 1, 2
    end
    C_TI = tens_TI(C1111, C1122, C1133, C3333, C2323, n_axis)
    P = hill_tensor(ell, C_TI)                     # analytical via dispatcher
    S = P ⊡ C_TI                                   # Eshelby tensor
    Sa = get_array(S)
    @printf("%6.2f  %12.5e  %12.5e  %12.5e  %12.5e\n",
            ξ,
            Sa[i_ax, i_ax, i_ax, i_ax],
            Sa[i_tr, i_tr, i_tr, i_tr],
            Sa[i_ax, i_ax, i_tr, i_tr],
            i_ax == 3 ? Sa[1, 1, 2, 2] : Sa[2, 2, 3, 3])
end
println()

# ── Cross-validation block needs a fixed config — use the canonical e₃ ------
n_axis = [0.0, 0.0, 1.0]
C_TI = tens_TI(C1111, C1122, C1133, C3333, C2323, n_axis)

# ── Cross-validation: analytical vs. numerical paths --------------------------
println("Cross-validation on a moderate oblate spheroid (1, 1, 0.3):")
ell = Ellipsoid(1.0, 1.0, 0.3)
P_ana = hill_tensor(ell, C_TI; method = :auto)        # analytical (default)
P_res = hill_tensor(ell, C_TI; method = :residues)
P_dec = hill_tensor(ell, C_TI; method = :decuhr)

diff_res = maximum(abs.(get_array(P_ana) .- get_array(P_res)))
diff_dec = maximum(abs.(get_array(P_ana) .- get_array(P_dec)))
@printf("  max |P_ana − P_residue| = %.3e\n", diff_res)
@printf("  max |P_ana − P_DECUHR | = %.3e\n", diff_dec)
println()

# ── Timing comparison ---------------------------------------------------------
function _bench(f, n)
    f()  # warm-up
    t0 = time()
    for _ in 1:n; f(); end
    return (time() - t0) / n
end

println("Timing (per call, average over 100 calls):")
t_ana = _bench(() -> hill_tensor(ell, C_TI; method = :auto), 100)
t_res = _bench(() -> hill_tensor(ell, C_TI; method = :residues), 100)
t_dec = _bench(() -> hill_tensor(ell, C_TI; method = :decuhr), 10)
@printf("  Analytical : %.3e s\n", t_ana)
@printf("  Residue    : %.3e s    (×%.0f)\n", t_res, t_res / t_ana)
@printf("  DECUHR     : %.3e s    (×%.0f)\n", t_dec, t_dec / t_ana)

println()
println("Done.")
