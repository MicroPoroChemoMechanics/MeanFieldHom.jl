using Test
using MeanFieldHom
using TensND
using ForwardDiff
using LinearAlgebra

# =============================================================================
#  test_sensitivities_alv.jl — ALV pipeline differentiates correctly via
#  `ForwardDiff` (set_param lens for fractions, closure-captured material
#  parameters for the kernel).  Each AD derivative is compared against a
#  central finite difference at `rtol ≤ 1e-7`.
# =============================================================================

const TIMES = collect(range(0.0, 2.0; length = 8))

function _build_law_M(k_M, μ_M, τ_K = 1.0, τ_μ = 0.5)
    function R_iso(t, tp)
        α = 3 * k_M * (1.0 + 4.0 * exp(-(t - tp) / τ_K))
        β = 2 * μ_M * (0.5 + 1.5 * exp(-(t - tp) / τ_μ))
        return TensISO{3}(α, β)
    end
    return ViscoLaw(R_iso, :relaxation)
end

const _C_INC = TensISO{3}(3 * 10.0, 2 * 4.0)

function _eff_mu_final(rve, scheme)
    R̃ = homogenize_alv(rve, scheme, :C; times = TIMES)
    _, β = iso_params_from_blocks(R̃)
    return β[end, end] / 2
end

function _build_rve_base(f::Real)
    rve = RVE(:M)
    add_matrix!(rve, Ellipsoid(1.0), Dict(:C => _build_law_M(1.0, 1.0)))
    add_phase!(
        rve, :I, Ellipsoid(1.0), Dict(:C => heaviside_law(_C_INC));
        fraction = f
    )
    return rve
end

@testset "ALV sensitivities — d/df via set_param lens" begin
    rve = _build_rve_base(0.2)
    f₀ = 0.2

    function eff_mu_vs_f(f, scheme)
        rve_f = set_param(rve, AmountParameter(:I), f)
        return _eff_mu_final(rve_f, scheme)
    end

    for sch in (Voigt(), Reuss(), Dilute(), MoriTanaka(), Maxwell())
        dμ_AD = ForwardDiff.derivative(f -> eff_mu_vs_f(f, sch), f₀)
        h = 1.0e-5
        dμ_FD = (eff_mu_vs_f(f₀ + h, sch) - eff_mu_vs_f(f₀ - h, sch)) / (2h)
        @test isapprox(dμ_AD, dμ_FD; rtol = 1.0e-7)
    end
end

@testset "ALV sensitivities — d/dμ_M via closure" begin
    function eff_mu_vs_μM(μ_M)
        rve = RVE(:M)
        add_matrix!(rve, Ellipsoid(1.0), Dict(:C => _build_law_M(1.0, μ_M)))
        add_phase!(
            rve, :I, Ellipsoid(1.0), Dict(:C => heaviside_law(_C_INC));
            fraction = 0.2
        )
        return _eff_mu_final(rve, MoriTanaka())
    end

    μM₀ = 1.0
    AD = ForwardDiff.derivative(eff_mu_vs_μM, μM₀)
    h = 1.0e-5
    FD = (eff_mu_vs_μM(μM₀ + h) - eff_mu_vs_μM(μM₀ - h)) / (2h)
    @test isapprox(AD, FD; rtol = 1.0e-6)
end

@testset "ALV sensitivities — gradient over (f, k_M, μ_M)" begin
    function eff_mu_vs_fkμ(p::AbstractVector)
        f, k_M, μ_M = p
        rve = RVE(:M)
        add_matrix!(rve, Ellipsoid(1.0), Dict(:C => _build_law_M(k_M, μ_M)))
        add_phase!(
            rve, :I, Ellipsoid(1.0), Dict(:C => heaviside_law(_C_INC));
            fraction = 0.2
        )
        rve_f = set_param(rve, AmountParameter(:I), f)
        return _eff_mu_final(rve_f, MoriTanaka())
    end

    p₀ = [0.2, 1.0, 1.0]
    ∇AD = ForwardDiff.gradient(eff_mu_vs_fkμ, p₀)
    h = 1.0e-5
    for i in 1:3
        eᵢ = [j == i ? 1.0 : 0.0 for j in 1:3]
        FD = (eff_mu_vs_fkμ(p₀ .+ h .* eᵢ) - eff_mu_vs_fkμ(p₀ .- h .* eᵢ)) / (2h)
        @test isapprox(∇AD[i], FD; rtol = 1.0e-6)
    end
end

@testset "ALV sensitivities — d/dτ_K (relaxation time inside kernel)" begin
    function eff_mu_vs_τK(τ_K)
        rve = RVE(:M)
        add_matrix!(rve, Ellipsoid(1.0), Dict(:C => _build_law_M(1.0, 1.0, τ_K, 0.5)))
        add_phase!(
            rve, :I, Ellipsoid(1.0), Dict(:C => heaviside_law(_C_INC));
            fraction = 0.2
        )
        return _eff_mu_final(rve, MoriTanaka())
    end

    τK₀ = 1.0
    AD = ForwardDiff.derivative(eff_mu_vs_τK, τK₀)
    h = 1.0e-5
    FD = (eff_mu_vs_τK(τK₀ + h) - eff_mu_vs_τK(τK₀ - h)) / (2h)
    @test isapprox(AD, FD; rtol = 1.0e-6)
end

# =============================================================================
#  Phase-2 extension : AD through the iterative / extra ALV schemes
#  (SelfConsistent, AsymmetricSelfConsistent, PonteCastanedaWillis,
#  DifferentialScheme) and through a GEOMETRY parameter (aspect ratio) —
#  previously blocked by hard-coded `Matrix{Float64}` containers.
# =============================================================================

@testset "ALV sensitivities — d/df through SC / ASC / PCW / DIFF" begin
    rve = _build_rve_base(0.2)
    f₀ = 0.2

    function eff_mu_vs_f(f, scheme)
        rve_f = set_param(rve, AmountParameter(:I), f)
        return _eff_mu_final(rve_f, scheme)
    end

    for sch in (
            SelfConsistent(), AsymmetricSelfConsistent(),
            PonteCastanedaWillis(), DifferentialScheme(),
        )
        dμ_AD = ForwardDiff.derivative(f -> eff_mu_vs_f(f, sch), f₀)
        h = 1.0e-5
        dμ_FD = (eff_mu_vs_f(f₀ + h, sch) - eff_mu_vs_f(f₀ - h, sch)) / (2h)
        @test isapprox(dμ_AD, dμ_FD; rtol = 1.0e-5)
    end
end

@testset "ALV sensitivities — d/dω geometry (aspect ratio) MT + SC" begin
    function eff_mu_vs_ω(ω, scheme)
        rve = RVE(:M)
        add_matrix!(rve, Ellipsoid(1.0), Dict(:C => _build_law_M(1.0, 1.0)))
        add_phase!(
            rve, :I, Spheroid(ω), Dict(:C => heaviside_law(_C_INC));
            fraction = 0.2
        )
        return _eff_mu_final(rve, scheme)
    end

    ω₀ = 3.0
    for sch in (MoriTanaka(), SelfConsistent())
        AD = ForwardDiff.derivative(ω -> eff_mu_vs_ω(ω, sch), ω₀)
        h = 1.0e-5
        FD = (eff_mu_vs_ω(ω₀ + h, sch) - eff_mu_vs_ω(ω₀ - h, sch)) / (2h)
        @test isapprox(AD, FD; rtol = 1.0e-5)
    end
end
