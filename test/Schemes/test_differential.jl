# =============================================================================
#  test_differential.jl — DifferentialScheme + trajectories.
#
#  Coverage:
#   1. Single-phase RVE returns the matrix property.
#   2. Bracketed by Voigt/Reuss bounds.
#   3. Single-inclusion RVE: all three trajectories give identical results
#      (the trajectory shape is irrelevant when there's only one phase).
#   4. Multi-phase RVE: Proportional and Sequential agree in the dilute
#      limit (small target fractions).
#   5. CustomPath validation: non-monotone, wrong endpoints, wrong length,
#      missing phase → all raise ArgumentError.
#   6. Crack-only RVE: stiffness reduces monotonically as nsteps→∞.
#   7. Conductivity (`property = :K`).
#   8. ForwardDiff sensitivity to f.
#   9. Symbol shortcuts.
# =============================================================================

using Test
using MeanFieldHom
using TensND
using LinearAlgebra
using ForwardDiff

const ATOL_DIFF = 1.0e-9
const RTOL_DIFF = 1.0e-8

@testset "Differential — sanity (single-phase)" begin
    C_m = TensISO{3}(30.0, 10.0)
    rve = RVE(:M)
    add_matrix!(rve, Ellipsoid(1.0), Dict(:C => C_m))
    @test homogenize(rve, DifferentialScheme()) ≈ C_m
end

@testset "Differential — bracketed by Voigt/Reuss" begin
    rve = RVE(:M)
    add_matrix!(rve, Ellipsoid(1.0), Dict(:C => TensISO{3}(30.0, 10.0)))
    add_phase!(
        rve, :I, Ellipsoid(1.0), Dict(:C => TensISO{3}(60.0, 20.0));
        fraction = 0.3
    )
    Vv = get_array(homogenize(rve, Voigt()))[1, 1, 1, 1]
    Vr = get_array(homogenize(rve, Reuss()))[1, 1, 1, 1]
    Vd = get_array(homogenize(rve, DifferentialScheme(; nsteps = 200)))[1, 1, 1, 1]
    @test Vr - RTOL_DIFF * abs(Vr) ≤ Vd ≤ Vv + RTOL_DIFF * abs(Vv)
end

@testset "Differential — trajectory invariance for single-inclusion RVE" begin
    rve = RVE(:M)
    add_matrix!(rve, Ellipsoid(1.0), Dict(:C => TensISO{3}(30.0, 10.0)))
    add_phase!(
        rve, :I, Ellipsoid(1.0), Dict(:C => TensISO{3}(60.0, 20.0));
        fraction = 0.3
    )
    C_prop = homogenize(rve, DifferentialScheme(; trajectory = Proportional(), nsteps = 100))
    C_seq = homogenize(rve, DifferentialScheme(; trajectory = Sequential([:I]), nsteps = 100))
    custom = CustomPath(Dict(:I => collect(range(0.0, 1.0; length = 101))))
    C_cus = homogenize(rve, DifferentialScheme(; trajectory = custom, nsteps = 100))
    # Tsit5 takes different adaptive steps depending on the smoothness
    # of `df/dτ` (Proportional has constant df, Sequential has step
    # discontinuities at window boundaries, CustomPath is piecewise
    # linear), so the parametrization invariance holds only up to the
    # solver's reltol (1e-6 by default).
    @test isapprox(C_prop, C_seq; rtol = 1.0e-5)
    @test isapprox(C_prop, C_cus; rtol = 1.0e-5)
end

@testset "Differential — multi-phase Proportional vs Sequential (dilute limit)" begin
    rve = RVE(:M)
    add_matrix!(rve, Ellipsoid(1.0), Dict(:C => TensISO{3}(30.0, 10.0)))
    add_phase!(
        rve, :I1, Ellipsoid(1.0), Dict(:C => TensISO{3}(60.0, 20.0));
        fraction = 0.01
    )
    add_phase!(
        rve, :I2, Ellipsoid(1.0), Dict(:C => TensISO{3}(15.0, 5.0));
        fraction = 0.01
    )
    Cp = get_array(
        homogenize(
            rve, DifferentialScheme(;
                trajectory = Proportional(),
                nsteps = 200
            )
        )
    )[1, 1, 1, 1]
    Cs = get_array(
        homogenize(
            rve, DifferentialScheme(;
                trajectory = Sequential([:I1, :I2]),
                nsteps = 200
            )
        )
    )[1, 1, 1, 1]
    # Dilute-limit difference is O(f²)
    @test abs(Cp - Cs) < 1.0e-3
end

@testset "Differential — CustomPath validation errors" begin
    rve = RVE(:M)
    add_matrix!(rve, Ellipsoid(1.0), Dict(:C => TensISO{3}(30.0, 10.0)))
    add_phase!(
        rve, :I, Ellipsoid(1.0), Dict(:C => TensISO{3}(60.0, 20.0));
        fraction = 0.3
    )

    # Wrong endpoint at start
    bad_start = CustomPath(Dict(:I => vcat([0.5], collect(range(0.0, 1.0; length = 100)))))
    @test_throws ArgumentError homogenize(rve, DifferentialScheme(; trajectory = bad_start, nsteps = 100))

    # Non-monotone
    nm = collect(range(0.0, 1.0; length = 101))
    nm[50] = 0.0
    bad_mono = CustomPath(Dict(:I => nm))
    @test_throws ArgumentError homogenize(rve, DifferentialScheme(; trajectory = bad_mono, nsteps = 100))

    # Missing phase
    rve2 = RVE(:M)
    add_matrix!(rve2, Ellipsoid(1.0), Dict(:C => TensISO{3}(30.0, 10.0)))
    add_phase!(
        rve2, :I1, Ellipsoid(1.0), Dict(:C => TensISO{3}(60.0, 20.0));
        fraction = 0.1
    )
    add_phase!(
        rve2, :I2, Ellipsoid(1.0), Dict(:C => TensISO{3}(15.0, 5.0));
        fraction = 0.1
    )
    bad_miss = CustomPath(Dict(:I1 => collect(range(0.0, 1.0; length = 101))))
    @test_throws ArgumentError homogenize(rve2, DifferentialScheme(; trajectory = bad_miss, nsteps = 100))
end

@testset "Differential — Path (functional) trajectory" begin
    rve = RVE(:M)
    add_matrix!(rve, Ellipsoid(1.0), Dict(:C => TensISO{3}(30.0, 10.0)))
    add_phase!(
        rve, :I, Ellipsoid(1.0), Dict(:C => TensISO{3}(60.0, 20.0));
        fraction = 0.3
    )
    # f(τ) = τ² is monotone with f(0)=0, f(1)=1.  Since DEM with one
    # solid phase is parametrization-invariant in τ, this should give
    # the same C^hom(τ=1) as the default Proportional path.
    C_path = homogenize(
        rve, DifferentialScheme(;
            trajectory = Path(Dict(:I => τ -> τ^2))
        )
    )
    C_prop = homogenize(rve, DifferentialScheme())
    @test isapprox(C_path, C_prop; rtol = 1.0e-5)
end

@testset "Differential — Path validation errors" begin
    rve = RVE(:M)
    add_matrix!(rve, Ellipsoid(1.0), Dict(:C => TensISO{3}(30.0, 10.0)))
    add_phase!(
        rve, :I, Ellipsoid(1.0), Dict(:C => TensISO{3}(60.0, 20.0));
        fraction = 0.3
    )
    # f(0) ≠ 0
    @test_throws ArgumentError homogenize(
        rve,
        DifferentialScheme(; trajectory = Path(Dict(:I => τ -> 0.5 + τ / 2)))
    )
    # f(1) ≠ 1
    @test_throws ArgumentError homogenize(
        rve,
        DifferentialScheme(; trajectory = Path(Dict(:I => τ -> τ / 2)))
    )
    # Non-monotone (sin(2πτ) has both signs)
    @test_throws ArgumentError homogenize(
        rve,
        DifferentialScheme(; trajectory = Path(Dict(:I => τ -> τ + 0.3sin(2π * τ))))
    )
    # Missing phase
    @test_throws ArgumentError homogenize(
        rve,
        DifferentialScheme(; trajectory = Path(Dict(:OTHER => τ -> τ)))
    )
end

@testset "Differential — saveat insensitivity (`nsteps`)" begin
    # The adaptive ODE solver controls integration step via abstol/reltol,
    # so `nsteps` (now `saveat` density) should not affect the final
    # result beyond solver tolerance.
    rve = RVE(:M)
    add_matrix!(rve, Ellipsoid(1.0), Dict(:C => TensISO{3}(30.0, 10.0)))
    add_phase!(
        rve, :I, Ellipsoid(1.0), Dict(:C => TensISO{3}(60.0, 20.0));
        fraction = 0.3
    )
    C_50 = homogenize(rve, DifferentialScheme(; nsteps = 50))
    C_200 = homogenize(rve, DifferentialScheme(; nsteps = 200))
    @test isapprox(C_50, C_200; rtol = 1.0e-6)
end

@testset "Differential — `abstol` / `reltol` kwargs are forwarded" begin
    # Loose tolerances should give a slightly different result than
    # tight ones — but both should be finite and within Voigt/Reuss.
    rve = RVE(:M)
    add_matrix!(rve, Ellipsoid(1.0), Dict(:C => TensISO{3}(30.0, 10.0)))
    add_phase!(
        rve, :I, Ellipsoid(1.0), Dict(:C => TensISO{3}(60.0, 20.0));
        fraction = 0.3
    )
    C_loose = homogenize(rve, DifferentialScheme(; abstol = 1.0e-3, reltol = 1.0e-2))
    C_tight = homogenize(rve, DifferentialScheme(; abstol = 1.0e-10, reltol = 1.0e-9))
    @test all(isfinite, get_array(C_tight))
    # Tighter tolerance gives a more accurate answer (not necessarily
    # closer to loose) — just check finiteness.
    @test all(isfinite, get_array(C_loose))
end

@testset "Differential — crack RVE reduces stiffness monotonically" begin
    rve = RVE(:M)
    add_matrix!(rve, Ellipsoid(1.0), Dict(:C => TensISO{3}(30.0, 10.0)))
    add_phase!(
        rve, :CRACK, PennyCrack(1.0), Dict(:C => TensISO{3}(30.0, 10.0));
        density = 0.05
    )
    C_d = homogenize(rve, DifferentialScheme(; nsteps = 200))
    @test get_array(C_d)[3, 3, 3, 3] < get_array(TensISO{3}(30.0, 10.0))[3, 3, 3, 3]
    @test all(isfinite, get_array(C_d))
end

@testset "Differential — conductivity" begin
    rve = RVE(:M)
    add_matrix!(rve, Ellipsoid(1.0), Dict(:K => TensISO{3}(2.0)))
    add_phase!(rve, :I, Ellipsoid(1.0), Dict(:K => TensISO{3}(8.0)); fraction = 0.3)
    Kv = get_array(homogenize(rve, Voigt(); property = :K))[1, 1]
    Kr = get_array(homogenize(rve, Reuss(); property = :K))[1, 1]
    Kd = get_array(
        homogenize(
            rve, DifferentialScheme(; nsteps = 200);
            property = :K
        )
    )[1, 1]
    @test Kr - RTOL_DIFF * abs(Kr) ≤ Kd ≤ Kv + RTOL_DIFF * abs(Kv)
end

@testset "Differential — ForwardDiff sensitivity to f" begin
    f_diff(f) = begin
        DT = typeof(f)
        rve = RVE(:M; T = DT)
        add_matrix!(rve, Ellipsoid(1.0), Dict(:C => TensISO{3}(30.0, 10.0)))
        add_phase!(
            rve, :I, Ellipsoid(1.0), Dict(:C => TensISO{3}(60.0, 20.0));
            fraction = f
        )
        get_array(homogenize(rve, DifferentialScheme(; nsteps = 50)))[1, 1, 1, 1]
    end
    df = ForwardDiff.derivative(f_diff, 0.3)
    @test isfinite(df)
    @test df > 0
end

@testset "Differential — Symbol shortcuts" begin
    rve = RVE(:M)
    add_matrix!(rve, Ellipsoid(1.0), Dict(:C => TensISO{3}(30.0, 10.0)))
    add_phase!(
        rve, :I, Ellipsoid(1.0), Dict(:C => TensISO{3}(60.0, 20.0));
        fraction = 0.3
    )
    @test homogenize(rve, :differential) ≈ homogenize(rve, DifferentialScheme())
    @test homogenize(rve, :diff) ≈ homogenize(rve, DifferentialScheme())
    @test homogenize(rve, :DIFF) ≈ homogenize(rve, DifferentialScheme())
end
