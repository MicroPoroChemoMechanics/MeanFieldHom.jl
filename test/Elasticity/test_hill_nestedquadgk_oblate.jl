# PR-4 validation: power change-of-variable in nestedquadgk for flat ellipsoids.
#
# Before this PR the nestedquadgk backend used a plain z ∈ [0, 1] outer
# integral.  For flat ellipsoids (ω ≪ 1) the integrand has a steep gradient
# near z = 1 that requires many subdivisions — or fails to converge at default
# maxiters.  The new path applies z(u) = 1 − (1 − u)^α with
# α = max(1, log₁₀(1/ω)), concentrating Gauss-Kronrod nodes near z = 1.
#
# Note: the Ellipsoid constructor already sorts semi-axes in descending order
# for real-valued types (via _sort_axes_and_basis), so η = a₂/a₁ ≤ 1 and
# ω = a₃/a₁ ≤ 1 are guaranteed on entry — the axis-sort inside the backend
# is always a no-op for Float64 inputs.
using Test
using MeanFieldHom
using TensND

const ATOL4 = 5.0e-6   # looser than spheroidal due to near-singular kernel at z≈1
const basis4 = TensND.CanonicalBasis{3, Float64}()

# Triclinic reference stiffness (same hand-picked matrix as test_anisotropic.jl)
const _KM_obl = [
    210.0 80.0 75.0 5.0 4.0 3.0;
     80.0 195.0 90.0 -2.0 3.0 -1.0;
     75.0 90.0 220.0 1.0 -2.0 2.0;
      5.0 -2.0 1.0 60.0 2.5 1.5;
      4.0 3.0 -2.0 2.5 65.0 -1.0;
      3.0 -1.0 2.0 1.5 -1.0 55.0
]


@testset "nestedquadgk — oblate ellipsoid ω = 0.01 (α = 2)" begin
    # Stored as (1.0, 1.0, 0.01) after constructor sort; ω = 0.01 → α = 2.
    # The change-of-variable concentrates quadrature near z = 1.
    C = TensND.invKM(_KM_obl, basis4)
    ell = Ellipsoid(1.0, 1.0, 0.01)

    P_nqg = hill_tensor(ell, C; method = :nestedquadgk)
    P_dec = hill_tensor(ell, C; method = :decuhr, reltol = 1.0e-9)

    diff = maximum(abs(P_nqg[i, j, k, l] - P_dec[i, j, k, l])
                   for i in 1:3, j in 1:3, k in 1:3, l in 1:3)
    @test diff < ATOL4
    @test all(isfinite(P_nqg[i, j, k, l]) for i in 1:3, j in 1:3, k in 1:3, l in 1:3)
end


@testset "nestedquadgk — very oblate ellipsoid ω = 0.001 (α = 3)" begin
    # ω = 0.001 → α = 3; tightest practical test for the change-of-variable.
    C = TensND.invKM(_KM_obl, basis4)
    ell = Ellipsoid(1.0, 1.0, 0.001)

    P_nqg = hill_tensor(ell, C; method = :nestedquadgk)
    P_dec = hill_tensor(ell, C; method = :decuhr, reltol = 1.0e-9)

    diff = maximum(abs(P_nqg[i, j, k, l] - P_dec[i, j, k, l])
                   for i in 1:3, j in 1:3, k in 1:3, l in 1:3)
    @test diff < ATOL4
    @test all(isfinite(P_nqg[i, j, k, l]) for i in 1:3, j in 1:3, k in 1:3, l in 1:3)
end


@testset "nestedquadgk — regular triaxial (α = 1, regression guard)" begin
    # ω = 1/3 → α = max(1, 0.48) = 1 (identity change-of-variable).
    # Checks that the new path does not regress accuracy on well-conditioned cases.
    C = TensND.invKM(_KM_obl, basis4)
    ell = Ellipsoid(3.0, 2.0, 1.0)

    P_nqg = hill_tensor(ell, C; method = :nestedquadgk, reltol = 1.0e-12)
    P_res = hill_tensor(ell, C; method = :residues)

    scale = maximum(abs(P_nqg[i, j, k, l])
                    for i in 1:3, j in 1:3, k in 1:3, l in 1:3)
    diff = maximum(abs(P_nqg[i, j, k, l] - P_res[i, j, k, l])
                   for i in 1:3, j in 1:3, k in 1:3, l in 1:3)
    @test diff < 1.0e-7 * scale
end
