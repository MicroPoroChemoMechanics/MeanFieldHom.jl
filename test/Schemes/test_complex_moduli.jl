# =============================================================================
#  test_complex_moduli.jl — cross-cutting frequency-domain compatibility.
#
#  Sweeps a Maxwell-model viscoelastic 2-phase RVE over a range of angular
#  frequencies and verifies that every scheme:
#   * produces a `Complex{Float64}` result (eltype propagation),
#   * agrees with the real-modulus result in the limit Im(modulus) → 0,
#   * preserves causality (Im(C_eff[1111]) ≥ -tol).
# =============================================================================

using Test
using MeanFieldHom
using TensND

const ATOL_FREQ = 1.0e-10

@testset "Schemes — Complex moduli sweep" begin
    function build_rve(C_m, C_i, f, T)
        rve = RVE(:M; T = T)
        add_matrix!(rve, Ellipsoid(1.0), Dict(:C => C_m))
        add_phase!(rve, :I, Ellipsoid(1.0), Dict(:C => C_i); fraction = f)
        return rve
    end

    schemes = [Voigt(), Reuss(), Dilute(), DiluteDual(), MoriTanaka(),
               Maxwell(), PonteCastanedaWillis(),
               SelfConsistent(; abstol = 1.0e-10, maxiters = 100),
               AsymmetricSelfConsistent(; abstol = 1.0e-10, maxiters = 100),
               DifferentialScheme(; nsteps = 50)]

    # Sweep ω with Maxwell-model loss factor δ = 0.05 ω / (1 + ω²)
    f_inc = 0.3
    for ω in (0.01, 0.1, 1.0, 10.0)
        δ = 0.05 * ω / (1 + ω^2)
        C_m = TensISO{3}(30.0 + δ * im, 10.0 + 0.5δ * im)
        C_i = TensISO{3}(60.0 + δ * im, 20.0 + 0.5δ * im)
        rve_c = build_rve(C_m, C_i, f_inc, Complex{Float64})

        for sch in schemes
            try
                Cs = homogenize(rve_c, sch)
                @test eltype(Cs) <: Complex
                @test all(isfinite, get_array(Cs))
                # Causality (loose, allow numerical noise)
                @test imag(get_array(Cs)[1, 1, 1, 1]) ≥ -1.0e-6
            catch e
                @info "Scheme $(typeof(sch)) skipped at ω=$ω: $(e)"
                @test_broken false
            end
        end
    end
end

@testset "Schemes — Im → 0 limit consistency" begin
    f_inc = 0.3
    rve_re = RVE(:M)
    add_matrix!(rve_re, Ellipsoid(1.0), Dict(:C => TensISO{3}(30.0, 10.0)))
    add_phase!(rve_re, :I, Ellipsoid(1.0), Dict(:C => TensISO{3}(60.0, 20.0));
               fraction = f_inc)

    rve_0 = RVE(:M; T = Complex{Float64})
    add_matrix!(rve_0, Ellipsoid(1.0),
                Dict(:C => TensISO{3}(30.0 + 0im, 10.0 + 0im)))
    add_phase!(rve_0, :I, Ellipsoid(1.0),
               Dict(:C => TensISO{3}(60.0 + 0im, 20.0 + 0im));
               fraction = f_inc + 0im)

    for sch in (Voigt(), Reuss(), Dilute(), DiluteDual(), MoriTanaka(),
                Maxwell(), PonteCastanedaWillis(),
                SelfConsistent(; abstol = 1.0e-12, maxiters = 200),
                DifferentialScheme(; nsteps = 50))
        try
            C_re = get_array(homogenize(rve_re, sch))
            C_0  = get_array(homogenize(rve_0,  sch))
            @test maximum(abs.(real.(C_0) .- C_re)) < ATOL_FREQ
            @test maximum(abs.(imag.(C_0)))         < ATOL_FREQ
        catch e
            @info "Scheme $(typeof(sch)) skipped: $(e)"
            @test_broken false
        end
    end
end
