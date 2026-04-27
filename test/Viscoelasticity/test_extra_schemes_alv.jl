using Test
using MeanFieldHom
using TensND
using LinearAlgebra

# =============================================================================
#  test_extra_schemes_alv.jl — PCW / ASC / DIFF ALV vs elastic limit and
#  vs the existing ALV schemes (consistency on identity / pure-matrix
#  edge cases).
# =============================================================================

const _to_mandel = MeanFieldHom.Viscoelasticity._tens_to_mandel66

function _setup_2phase_elastic(; k_M = 10.0, μ_M = 4.0,
                                k_I = 20.0, μ_I = 8.0,
                                f_I = 0.2, n_times = 4)
    C_M_t = TensISO{3}(3 * k_M, 2 * μ_M)
    C_I_t = TensISO{3}(3 * k_I, 2 * μ_I)
    times = collect(range(0.0, 1.0; length = n_times))
    return (; C_M_t, C_I_t, times, f_I, f_M = 1 - f_I)
end

function _build_alv(ctx)
    rve = RVE(:M)
    add_matrix!(rve, Ellipsoid(1.0, 1.0, 1.0),
                Dict(:C => heaviside_law(ctx.C_M_t)))
    add_phase!(rve, :I, Ellipsoid(1.0, 1.0, 1.0),
               Dict(:C => heaviside_law(ctx.C_I_t)); fraction = ctx.f_I)
    return rve
end

function _build_el(ctx)
    rve = RVE(:M)
    add_matrix!(rve, Ellipsoid(1.0, 1.0, 1.0), Dict(:C => ctx.C_M_t))
    add_phase!(rve, :I, Ellipsoid(1.0, 1.0, 1.0), Dict(:C => ctx.C_I_t);
               fraction = ctx.f_I)
    return rve
end

function _check_alv_elastic(C_alv::AbstractMatrix, M_ref::AbstractMatrix,
                            n::Int; rtol = 1.0e-12, atol = 1.0e-12)
    for i in 1:n
        rows = (6 * (i - 1) + 1):(6 * i)
        @test isapprox(C_alv[rows, rows], M_ref; rtol = rtol, atol = atol)
        for j in 1:(i - 1)
            cols = (6 * (j - 1) + 1):(6 * j)
            @test maximum(abs, C_alv[rows, cols]) ≤ atol
        end
    end
end

@testset "asymmetric_self_consistent_alv — elastic limit (sphere)" begin
    ctx = _setup_2phase_elastic()
    n = length(ctx.times)
    M_ref = _to_mandel(homogenize(_build_el(ctx), AsymmetricSelfConsistent(), :C))
    C_alv = homogenize_alv(_build_alv(ctx), AsymmetricSelfConsistent(), :C;
                            times = ctx.times)
    _check_alv_elastic(C_alv, M_ref, n; atol = 1.0e-9, rtol = 1.0e-9)
end

@testset "pcw_alv — elastic limit (sphere distribution)" begin
    ctx = _setup_2phase_elastic()
    n = length(ctx.times)
    M_ref = _to_mandel(homogenize(_build_el(ctx), PonteCastanedaWillis(), :C))
    C_alv = homogenize_alv(_build_alv(ctx), PonteCastanedaWillis(), :C;
                            times = ctx.times)
    _check_alv_elastic(C_alv, M_ref, n; atol = 1.0e-12)
end

@testset "differential_alv — elastic limit (sphere)" begin
    ctx = _setup_2phase_elastic()
    n = length(ctx.times)
    sch = DifferentialScheme(; nsteps = 50)
    M_ref = _to_mandel(homogenize(_build_el(ctx), sch, :C))
    C_alv = homogenize_alv(_build_alv(ctx), sch, :C; times = ctx.times)
    _check_alv_elastic(C_alv, M_ref, n; atol = 1.0e-12)
end

@testset "PCW vs Maxwell — equivalence in single-shape case" begin
    # PCW with default UniformDistribution(unit sphere) ≡ Maxwell with
    # spherical distribution shape.
    ctx = _setup_2phase_elastic()
    rve = _build_alv(ctx)
    R_pcw = homogenize_alv(rve, PonteCastanedaWillis(), :C; times = ctx.times)
    R_max = homogenize_alv(rve, Maxwell(), :C; times = ctx.times)
    @test isapprox(R_pcw, R_max; atol = 1.0e-12)
end

@testset "ASC vs SC — same fixed point in elastic limit" begin
    ctx = _setup_2phase_elastic()
    rve = _build_alv(ctx)
    R_sc  = homogenize_alv(rve, SelfConsistent(), :C; times = ctx.times)
    R_asc = homogenize_alv(rve, AsymmetricSelfConsistent(), :C;
                            times = ctx.times)
    @test isapprox(R_sc, R_asc; atol = 1.0e-9, rtol = 1.0e-9)
end

@testset "Differential — independent of nsteps in elastic limit" begin
    ctx = _setup_2phase_elastic()
    rve = _build_alv(ctx)
    R20  = homogenize_alv(rve, DifferentialScheme(; nsteps = 20),  :C; times = ctx.times)
    R100 = homogenize_alv(rve, DifferentialScheme(; nsteps = 100), :C; times = ctx.times)
    # Higher nsteps should converge — finite difference in nsteps is small.
    @test isapprox(R20, R100; atol = 1.0e-3, rtol = 1.0e-3)
end
