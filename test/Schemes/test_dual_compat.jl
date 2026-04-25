# =============================================================================
#  test_dual_compat.jl — cross-cutting ForwardDiff sensitivity test.
#
#  Verifies that every scheme is differentiable through:
#   * the volume fraction of an inclusion phase,
#   * the moduli of a phase (E-like, μ-like).
#
#  Compares the AD derivative against a centred-finite-difference reference
#  to relative tolerance 1e-5.
# =============================================================================

using Test
using MeanFieldHom
using TensND
using ForwardDiff

const RTOL_AD = 1.0e-5

@testset "Schemes — ForwardDiff sensitivity to f (every scheme)" begin
    schemes = [Voigt(), Reuss(), Dilute(), DiluteDual(), MoriTanaka(),
               Maxwell(), PonteCastanedaWillis(),
               SelfConsistent(; abstol = 1.0e-12, maxiters = 200),
               AsymmetricSelfConsistent(; abstol = 1.0e-12, maxiters = 200),
               DifferentialScheme(; nsteps = 50)]

    for sch in schemes
        f_eff(f) = begin
            DT = typeof(f)
            rve = RVE(:M; T = DT)
            add_matrix!(rve, Ellipsoid(1.0), Dict(:C => TensISO{3}(30.0, 10.0)))
            add_phase!(rve, :I, Ellipsoid(1.0), Dict(:C => TensISO{3}(60.0, 20.0));
                       fraction = f)
            return get_array(homogenize(rve, sch))[1, 1, 1, 1]
        end
        df_ad = ForwardDiff.derivative(f_eff, 0.25)
        # Centred FD reference
        h = 1.0e-6
        df_fd = (f_eff(0.25 + h) - f_eff(0.25 - h)) / (2h)
        @test isfinite(df_ad)
        @test isapprox(df_ad, df_fd; rtol = RTOL_AD, atol = 1.0e-7)
    end
end

@testset "Schemes — ForwardDiff sensitivity to a modulus (every scheme)" begin
    schemes = [Voigt(), Reuss(), Dilute(), DiluteDual(), MoriTanaka(),
               Maxwell(), PonteCastanedaWillis(),
               SelfConsistent(; abstol = 1.0e-12, maxiters = 200),
               AsymmetricSelfConsistent(; abstol = 1.0e-12, maxiters = 200),
               DifferentialScheme(; nsteps = 50)]

    for sch in schemes
        f_eff(α) = begin
            rve = RVE(:M)
            add_matrix!(rve, Ellipsoid(1.0), Dict(:C => TensISO{3}(α, 10.0)))
            add_phase!(rve, :I, Ellipsoid(1.0), Dict(:C => TensISO{3}(60.0, 20.0));
                       fraction = 0.3)
            return get_array(homogenize(rve, sch))[1, 1, 1, 1]
        end
        df_ad = ForwardDiff.derivative(f_eff, 30.0)
        h = 1.0e-4
        df_fd = (f_eff(30.0 + h) - f_eff(30.0 - h)) / (2h)
        @test isfinite(df_ad)
        @test isapprox(df_ad, df_fd; rtol = RTOL_AD, atol = 1.0e-7)
    end
end
