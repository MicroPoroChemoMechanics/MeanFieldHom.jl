using Test
using MeanFieldHom.LayeredSpheroids: legendre_odd

# =============================================================================
#  test_legendre.jl — associated Legendre recurrences (`legendre.jl`).
#
#  Checked against closed-form low-degree polynomials (P₁, P₃, P₅, P₁¹,
#  Q₁, Q₁¹) and, for the higher-degree / Q¹ entries with no simple
#  closed form at hand, against a central finite-difference derivative
#  — cheap and independent of the recurrence itself.
# =============================================================================

@testset "legendre_odd — P0 (plain Legendre) against closed forms" begin
    x = 1.7
    P0, dP0 = legendre_odd(:P0, x, 3)   # degrees 1, 3, 5
    @test P0[1] ≈ x                         # P₁(x) = x
    @test P0[2] ≈ (5x^3 - 3x) / 2            # P₃(x)
    @test P0[3] ≈ (63x^5 - 70x^3 + 15x) / 8  # P₅(x)
    @test dP0[1] ≈ 1.0
    @test dP0[2] ≈ (15x^2 - 3) / 2
end

@testset "legendre_odd — Q0 against Q₁(x) = x·arccoth(x) - 1" begin
    x = 1.7
    ax = atanh(1 / x)
    Q0, dQ0 = legendre_odd(:Q0, x, 1)
    @test Q0[1] ≈ x * ax - 1
    @test dQ0[1] ≈ ax - x / (x^2 - 1)
end

@testset "legendre_odd — P1 (q branch) against P₁¹(q) = √(q²-1)" begin
    x = 1.7
    P1, dP1 = legendre_odd(:P1, x, 1)
    @test P1[1] ≈ sqrt(x^2 - 1)
    @test dP1[1] ≈ x / sqrt(x^2 - 1)
end

@testset "legendre_odd — P1p (p branch) against P₁¹(p) = -√(1-p²)" begin
    p = 0.3
    P1p, dP1p = legendre_odd(:P1p, p, 1)
    @test P1p[1] ≈ -sqrt(1 - p^2)
    @test dP1p[1] ≈ -p / (-sqrt(1 - p^2))
end

@testset "legendre_odd — derivative tables match finite differences" begin
    # Coarse sanity check only (higher degrees grow fast, so a fixed step
    # trades truncation vs. cancellation error): the precise validation is
    # the closed-form checks above and the coupling-matrix cross-check
    # against the independent BigFloat series in `test_coupling.jl`.
    h = 1.0e-5
    for x in (1.3, 1.7, 3.3), Nseries in (1, 2, 4)
        for kind in (:P0, :Q0, :P1, :Q1)
            vals_p, ders = legendre_odd(kind, x, Nseries)
            vals_hi, _ = legendre_odd(kind, x + h, Nseries)
            vals_lo, _ = legendre_odd(kind, x - h, Nseries)
            fd = (vals_hi .- vals_lo) ./ (2h)
            @test all(isapprox.(fd, ders; rtol = 1.0e-4, atol = 1.0e-6))
        end
        p = 0.3
        P1p, dP1p = legendre_odd(:P1p, p, Nseries)
        P1p_hi, _ = legendre_odd(:P1p, p + h, Nseries)
        P1p_lo, _ = legendre_odd(:P1p, p - h, Nseries)
        fd = (P1p_hi .- P1p_lo) ./ (2h)
        @test all(isapprox.(fd, dP1p; rtol = 1.0e-4, atol = 1.0e-6))
    end
end

@testset "legendre_odd — generic over Complex and BigFloat" begin
    q = im * 2.3   # oblate confocal parameter q = iτ
    P0, _ = legendre_odd(:P0, q, 3)
    @test eltype(P0) <: Complex
    @test P0[1] ≈ q

    xb = big"1.7"
    P0b, _ = legendre_odd(:P0, xb, 3)
    @test eltype(P0b) == BigFloat
    @test P0b[1] ≈ xb
end
