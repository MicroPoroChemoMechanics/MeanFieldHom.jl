# =============================================================================
#  test_self_consistent.jl — SelfConsistent + AsymmetricSelfConsistent.
#
#  Coverage:
#   1. Single-phase RVE returns the matrix property exactly.
#   2. Iso 2-phase composite : SC bracketed by Voigt/Reuss.
#   3. SC fixed-point self-consistency : `step(C_eff) ≈ C_eff` to abstol.
#   4. ASC ≡ SC when the matrix is the soft phase (stiffness-form path).
#   5. ASC handles inclusion-soft RVE through the compliance-form path.
#   6. Conductivity (`property = :K`) — same recipes via gradient_gradient_loc.
#   7. ForwardDiff sensitivity through the volume fraction.
#   8. NewtonDefault (built-in, dependency-free) agrees with AndersonDefault.
#   9. Symbol shortcuts.
#
#  NonlinearSolve.jl-backed algorithms (NewtonRaphson, TrustRegion, …)
#  and AutoNonlinear are covered separately in
#  `test_self_consistent_nls.jl`.
# =============================================================================

using Test
using MeanFieldHom
using TensND
using LinearAlgebra
using ForwardDiff

const ATOL_SC = 1.0e-9
const RTOL_SC = 1.0e-8

@testset "SelfConsistent — sanity (single-phase)" begin
    C_m = TensISO{3}(30.0, 10.0)
    rve = RVE(:M)
    add_matrix!(rve, Ellipsoid(1.0), Dict(:C => C_m))
    @test homogenize(rve, SelfConsistent()) ≈ C_m
    @test homogenize(rve, AsymmetricSelfConsistent()) ≈ C_m
end

@testset "SelfConsistent — bracketed by Voigt/Reuss" begin
    rve = RVE(:M)
    add_matrix!(rve, Ellipsoid(1.0), Dict(:C => TensISO{3}(30.0, 10.0)))
    add_phase!(
        rve, :I, Ellipsoid(1.0), Dict(:C => TensISO{3}(60.0, 20.0));
        fraction = 0.3
    )

    Vv = get_array(homogenize(rve, Voigt()))[1, 1, 1, 1]
    Vr = get_array(homogenize(rve, Reuss()))[1, 1, 1, 1]
    Vsc = get_array(homogenize(rve, SelfConsistent()))[1, 1, 1, 1]
    @test Vr - RTOL_SC * abs(Vr) ≤ Vsc ≤ Vv + RTOL_SC * abs(Vv)
end

@testset "SelfConsistent — fixed-point self-consistency" begin
    rve = RVE(:M)
    add_matrix!(rve, Ellipsoid(1.0), Dict(:C => TensISO{3}(30.0, 10.0)))
    add_phase!(
        rve, :I, Ellipsoid(1.0), Dict(:C => TensISO{3}(60.0, 20.0));
        fraction = 0.3
    )
    C_eff = homogenize(rve, SelfConsistent(; abstol = 1.0e-12, maxiters = 200))
    # Verify : one more SC step on C_eff itself should return ≈ C_eff
    step_once = MeanFieldHom.Schemes._sc_step(rve, C_eff, :C)
    @test maximum(abs.(get_array(step_once) .- get_array(C_eff))) < 1.0e-9
end

@testset "AsymmetricSelfConsistent — matches SC when matrix is soft" begin
    # Inclusion stiffer than matrix → stiffness form preferred → ASC ≡ SC
    rve = RVE(:M)
    add_matrix!(rve, Ellipsoid(1.0), Dict(:C => TensISO{3}(30.0, 10.0)))
    add_phase!(
        rve, :I, Ellipsoid(1.0), Dict(:C => TensISO{3}(60.0, 20.0));
        fraction = 0.3
    )
    @test homogenize(rve, AsymmetricSelfConsistent()) ≈
        homogenize(rve, SelfConsistent())
end

@testset "AsymmetricSelfConsistent — uses compliance form when matrix is stiff" begin
    # Soft inclusion in stiff matrix → ASC switches to compliance-form
    rve = RVE(:M)
    add_matrix!(rve, Ellipsoid(1.0), Dict(:C => TensISO{3}(60.0, 20.0)))
    add_phase!(
        rve, :I, Ellipsoid(1.0), Dict(:C => TensISO{3}(3.0, 1.0));
        fraction = 0.3
    )
    C_asc = homogenize(rve, AsymmetricSelfConsistent())
    # Bracketed by Voigt/Reuss
    Vv = get_array(homogenize(rve, Voigt()))[1, 1, 1, 1]
    Vr = get_array(homogenize(rve, Reuss()))[1, 1, 1, 1]
    Va = get_array(C_asc)[1, 1, 1, 1]
    @test Vr - RTOL_SC * abs(Vr) ≤ Va ≤ Vv + RTOL_SC * abs(Vv)
end

@testset "SelfConsistent — conductivity" begin
    rve = RVE(:M)
    add_matrix!(rve, Ellipsoid(1.0), Dict(:K => TensISO{3}(2.0)))
    add_phase!(
        rve, :I, Ellipsoid(1.0), Dict(:K => TensISO{3}(8.0));
        fraction = 0.3
    )
    Kv = get_array(homogenize(rve, Voigt(); property = :K))[1, 1]
    Kr = get_array(homogenize(rve, Reuss(); property = :K))[1, 1]
    Ksc = get_array(homogenize(rve, SelfConsistent(); property = :K))[1, 1]
    @test Kr - RTOL_SC * abs(Kr) ≤ Ksc ≤ Kv + RTOL_SC * abs(Kv)
end

@testset "SelfConsistent — ForwardDiff sensitivity to f" begin
    f_sc(f) = begin
        DT = typeof(f)
        rve = RVE(:M; T = DT)
        add_matrix!(rve, Ellipsoid(1.0), Dict(:C => TensISO{3}(30.0, 10.0)))
        add_phase!(
            rve, :I, Ellipsoid(1.0), Dict(:C => TensISO{3}(60.0, 20.0));
            fraction = f
        )
        get_array(homogenize(rve, SelfConsistent()))[1, 1, 1, 1]
    end
    df = ForwardDiff.derivative(f_sc, 0.3)
    @test isfinite(df)
    @test df > 0
end

@testset "SelfConsistent — NewtonDefault works out of the box" begin
    # Since v0.7.0 ForwardDiff is a strong dependency and the built-in
    # `NewtonDefault` SC solver ships with the package — quadratic
    # convergence on iso / TI / ortho canonical components.
    rve = RVE(:M)
    add_matrix!(rve, Ellipsoid(1.0), Dict(:C => TensISO{3}(30.0, 10.0)))
    add_phase!(
        rve, :I, Ellipsoid(1.0), Dict(:C => TensISO{3}(60.0, 20.0));
        fraction = 0.3
    )
    C_anderson = homogenize(rve, SelfConsistent(; algorithm = AndersonDefault()))
    C_newton = homogenize(rve, SelfConsistent(; algorithm = NewtonDefault()))
    @test isapprox(C_anderson, C_newton; atol = 1.0e-6, rtol = 1.0e-6)
end

@testset "SelfConsistent — NewtonDefault ForwardDiff sensitivity (non-matrix phase)" begin
    # Regression test: `x0 = matrix_property(rve, p)` can stay `Float64`
    # while `step(x0)` promotes to `Dual` internally, whenever the
    # differentiated parameter lives on a phase OTHER than the one `x0`
    # is built from (an inclusion modulus, or any volume fraction). This
    # used to crash `_solve_sc(::NewtonDefault, …)` (`Tref` was derived
    # from `eltype(p0)` alone); it must now match a central finite
    # difference for every such parameter.
    rve = RVE(:M)
    add_matrix!(rve, Ellipsoid(1.0), Dict(:C => TensISO{3}(30.0, 10.0)))
    add_phase!(
        rve, :I, Ellipsoid(1.0), Dict(:C => TensISO{3}(60.0, 20.0));
        fraction = 0.3
    )
    idxC = C -> get_array(C)[1, 1, 1, 1]

    d_incl_modulus = derivative(
        rve, SelfConsistent(; algorithm = NewtonDefault()),
        property(:I, :C, :bulk); indexer = idxC
    )
    d_fraction = derivative(
        rve, SelfConsistent(; algorithm = NewtonDefault()),
        amount(:I); indexer = idxC
    )

    function f_modulus(K_I)
        r = RVE(:M)
        add_matrix!(r, Ellipsoid(1.0), Dict(:C => TensISO{3}(30.0, 10.0)))
        add_phase!(r, :I, Ellipsoid(1.0), Dict(:C => TensISO{3}(K_I, 20.0)); fraction = 0.3)
        return idxC(homogenize(r, SelfConsistent(; algorithm = NewtonDefault())))
    end
    function f_fraction(f)
        r = RVE(:M)
        add_matrix!(r, Ellipsoid(1.0), Dict(:C => TensISO{3}(30.0, 10.0)))
        add_phase!(r, :I, Ellipsoid(1.0), Dict(:C => TensISO{3}(60.0, 20.0)); fraction = f)
        return idxC(homogenize(r, SelfConsistent(; algorithm = NewtonDefault())))
    end
    h = 1.0e-5
    d_fd_modulus = (f_modulus(60.0 + h) - f_modulus(60.0 - h)) / (2h)
    d_fd_fraction = (f_fraction(0.3 + h) - f_fraction(0.3 - h)) / (2h)

    @test isapprox(d_incl_modulus, d_fd_modulus; rtol = 1.0e-5)
    @test isapprox(d_fraction, d_fd_fraction; rtol = 1.0e-5)
end

@testset "SelfConsistent / ASC — Symbol shortcuts" begin
    rve = RVE(:M)
    add_matrix!(rve, Ellipsoid(1.0), Dict(:C => TensISO{3}(30.0, 10.0)))
    add_phase!(
        rve, :I, Ellipsoid(1.0), Dict(:C => TensISO{3}(60.0, 20.0));
        fraction = 0.3
    )
    @test homogenize(rve, :sc) ≈ homogenize(rve, SelfConsistent())
    @test homogenize(rve, :SC) ≈ homogenize(rve, SelfConsistent())
    @test homogenize(rve, :self_consistent) ≈ homogenize(rve, SelfConsistent())
    @test homogenize(rve, :asc) ≈ homogenize(rve, AsymmetricSelfConsistent())
end
