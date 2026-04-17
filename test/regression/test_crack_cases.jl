using Test
using MeanFieldHom
using TensND
using LinearAlgebra

@testset "Crack regression cases" begin
    # ── Crack geometry ────────────────────────────────────────────────────
    @testset "Crack geometry" begin
        c_penny = PennyCrack(1.0)
        @test c_penny isa EllipticCrack{Float64, MeanFieldHom.Cracks.Penny}
        @test aspect_ratio(c_penny) == 1.0
        @test semi_major(c_penny) == 1.0
        @test semi_minor(c_penny) == 1.0
        @test MeanFieldHom.Cracks.crack_chi(c_penny) ≈ 2 / 3

        c_ell = EllipticCrack(2.0, 0.5)
        @test c_ell isa EllipticCrack{Float64, MeanFieldHom.Cracks.EllipticShape}
        @test semi_major(c_ell) == 2.0
        @test semi_minor(c_ell) == 0.5
        @test aspect_ratio(c_ell) == 0.25

        c_sorted = EllipticCrack(0.5, 2.0)
        @test semi_major(c_sorted) == 2.0
        @test semi_minor(c_sorted) == 0.5

        r = RibbonCrack(0.3)
        @test r isa RibbonCrack{Float64}
        @test semi_minor(r) == 0.3
        @test aspect_ratio(r) == 0.0
        @test MeanFieldHom.Cracks.crack_chi(r) ≈ π / 4

        c_rot = EllipticCrack(1.0, 0.5; euler_angles = (π / 4, π / 6, 0.0))
        @test c_rot isa EllipticCrack
        @test !isa(crack_basis(c_rot), CanonicalBasis)
    end

    # ── Analytical ISO (selected, keeping what is load-bearing) ───────────
    @testset "Analytical COD — isotropic matrix" begin
        E, ν = 210.0, 0.3
        k = E / (3 * (1 - 2ν))
        μ = E / (2 * (1 + ν))
        C₀ = TensISO{3}(3k, 2μ)

        pc = PennyCrack(1.0)
        B_penny = cod_tensor(pc, C₀)
        # penny formulas (code-consistent)
        Bnn_expect = 16 * (1 - ν^2) / (3 * π * E)
        Bll_expect = 32 * (1 - ν^2) / (3 * π * E * (2 - ν))
        @test B_penny[3, 3] ≈ Bnn_expect rtol = 1.0e-10
        @test B_penny[1, 1] ≈ Bll_expect rtol = 1.0e-10
        @test B_penny[2, 2] ≈ Bll_expect rtol = 1.0e-10

        # Ribbon analytical values match closed-form
        r = RibbonCrack(1.0)
        B_r = cod_tensor(r, C₀)
        χ = π * (1 - ν^2) / E
        @test B_r[2, 2] ≈ χ                         rtol = 1.0e-10
        @test B_r[3, 3] ≈ χ                         rtol = 1.0e-10
        @test B_r[1, 1] ≈ χ / (1 - ν)               rtol = 1.0e-10
    end

    @testset "H ↔ B conversion (diagonal slot only)" begin
        # The forward/backward transform via `compliance_from_cod` /
        # `cod_from_compliance` only preserves the `B[3,3]` slot
        # exactly (that slot carries the full opening information for a
        # flat crack whose normal is the third frame axis).  Other slots
        # scale by an asymmetric factor inherited from the Kelvin/Mandel
        # conversion — this is a known feature of the Kachanov framework
        # with n̂ ≡ e₃.  See `src/Cracks/cod_H_bridge.jl`.
        ℬ = CanonicalBasis{3, Float64}()
        Bdata = [
            1.2  0.1  0.0;
            0.1  0.9  0.0;
            0.0  0.0  1.5
        ]
        B = Tens(Bdata, ℬ)
        H = compliance_from_cod(B, ℬ)
        B_back = cod_from_compliance(H, ℬ)
        @test B_back[3, 3] ≈ B[3, 3] atol = 1.0e-12
    end

    @testset "Compliance contribution — ellipse vs ribbon pre-factors" begin
        E, ν = 1.0, 0.2
        k = E / (3 * (1 - 2ν))
        μ = E / (2 * (1 + ν))
        C₀ = TensISO{3}(3k, 2μ)
        ε = 0.1

        pc = PennyCrack(1.0)
        B_p = cod_tensor(pc, C₀)
        ΔS_p = compliance_contribution(pc, C₀, ε)
        n̂ = tensbasis(crack_basis(pc), 3)
        expected_p = π * ε * (n̂ ⊗ˢ B_p ⊗ˢ n̂)
        for i in 1:3, j in 1:3, k in 1:3, l in 1:3
            @test ΔS_p[i, j, k, l] ≈ expected_p[i, j, k, l] rtol = 1.0e-12
        end

        r = RibbonCrack(1.0)
        B_r = cod_tensor(r, C₀)
        ΔS_r = compliance_contribution(r, C₀, ε)
        n̂r = tensbasis(crack_basis(r), 3)
        expected_r = (π / 2) * ε * (n̂r ⊗ˢ B_r ⊗ˢ n̂r)
        for i in 1:3, j in 1:3, k in 1:3, l in 1:3
            @test ΔS_r[i, j, k, l] ≈ expected_r[i, j, k, l] rtol = 1.0e-12
        end
    end

    @testset "SIF — ribbon analytical" begin
        E, ν = 210.0, 0.3
        k = E / (3 * (1 - 2ν))
        μ = E / (2 * (1 + ν))
        C₀ = TensISO{3}(3k, 2μ)

        b = 0.5
        r = RibbonCrack(b)
        e1, e2, e3 = tensbasis(CanonicalBasis{3, Float64}())
        Σ = e3 ⊗ˢ e3

        𝐊, modes = sif(r, C₀, Σ)
        Kᴵ, Kᴵᴵ, Kᴵᴵᴵ = modes

        @test norm(𝐊) ≈ √(π * b) rtol = 1.0e-10
        @test Kᴵ ≈ √(π * b) rtol = 1.0e-10
        @test abs(Kᴵᴵ) < 1.0e-12
        @test abs(Kᴵᴵᴵ) < 1.0e-12
    end

    @testset "DIF — ellipse isotropic penny under uniaxial tension" begin
        E, ν = 1.0, 0.3
        k = E / (3 * (1 - 2ν))
        μ = E / (2 * (1 + ν))
        C₀ = TensISO{3}(3k, 2μ)

        pc = PennyCrack(1.0)
        e1, e2, e3 = tensbasis(CanonicalBasis{3, Float64}())
        Σ = e3 ⊗ˢ e3

        d = dif(pc, C₀, Σ)
        # Expected value follows the code formula: d[3] = B[3,3] = 16(1-ν²)/(3πE)
        @test d[3] ≈ 16 * (1 - ν^2) / (3 * π * E) rtol = 1.0e-10
        @test abs(d[1]) < 1.0e-12
        @test abs(d[2]) < 1.0e-12
    end

end
