using Test
using MeanFieldHom
using TensND
using LinearAlgebra

# Cross-check the analytical thermal COD scalar against the Hill 2nd-order
# tensor Taylor limit: lim_{ε→0} ε·a_max·(K₀ − K₀ P(ε) K₀)⁻¹ equals the
# *old* dilute resistivity contribution ΔR at ε=1, which under the new
# convention reads delta_resistivity(crack, R, 1).

function hill_limit_R(crack, K₀; ε_small = 1.0e-8)
    a_max = crack isa EllipticCrack ? max(crack.a, crack.b) : crack.b
    ell = crack isa EllipticCrack ?
        Ellipsoid(crack.a, crack.b, ε_small * a_max, crack_basis(crack)) :
        Ellipsoid(1.0e8 * crack.b, crack.b, ε_small * crack.b, crack_basis(crack))
    P = hill_tensor(ell, K₀)
    P_arr = [P[i, j] for i in 1:3, j in 1:3]
    K_arr = [K₀[i, j] for i in 1:3, j in 1:3]
    Λ = K_arr - K_arr * P_arr * K_arr
    return ε_small * a_max * inv(Λ)
end

@testset "Cracks / conductivity — thermal COD scalar b and R" begin

    @testset "Iso matrix, penny crack" begin
        for k₀ in (1.0, 2.5, 7.3)
            K₀ = TensISO{3}(k₀)
            crack = PennyCrack(1.0)
            b = cod_tensor(crack, K₀)
            @test b ≈ 2 / (π^2 * k₀) rtol = 1.0e-13
            R = compliance_contribution(crack, K₀)
            # R = (3/4) b (n̂ ⊗ n̂)
            @test R[3, 3] ≈ 3 / (2 * π^2 * k₀) rtol = 1.0e-13
            @test R[1, 1] ≈ 0 atol = 1.0e-14
            @test R[2, 2] ≈ 0 atol = 1.0e-14
            # delta_resistivity recovers the dilute correction ΔR
            ΔR = delta_resistivity(crack, R, 1.0)
            @test ΔR[3, 3] ≈ 2 / (π * k₀) rtol = 1.0e-13
        end
    end

    @testset "Iso matrix, non-penny elliptic" begin
        k₀ = 2.0
        K₀ = TensISO{3}(k₀)
        for η in (0.75, 0.5, 0.3)
            a = 1.0
            b_ax = a * η
            crack = EllipticCrack(a, b_ax)
            b = cod_tensor(crack, K₀)
            # Analytical: b = η/(π k₀ ℰ_η)
            k² = 1 - η^2
            ℰ = MeanFieldHom.Elliptic.ell_E(k²)
            @test b ≈ η / (π * k₀ * ℰ) rtol = 1.0e-13
        end
    end

    @testset "Iso matrix, ribbon crack" begin
        for k₀ in (1.0, 2.5)
            K₀ = TensISO{3}(k₀)
            crack = RibbonCrack(0.5)
            b = cod_tensor(crack, K₀)
            @test b ≈ 2 / (π * k₀) rtol = 1.0e-13
            R = compliance_contribution(crack, K₀)
            # R = (2/π) b (n̂ ⊗ n̂)
            @test R[3, 3] ≈ 4 / (π^2 * k₀) rtol = 1.0e-13
            # ΔR = π ε²ᵈ R at ε²ᵈ = 1: ΔR[3,3] = π × 4/(π² k₀) = 4/(π k₀)
            ΔR = delta_resistivity(crack, R, 1.0)
            @test ΔR[3, 3] ≈ 4 / (π * k₀) rtol = 1.0e-13
        end
    end

    @testset "Aligned TI, penny crack (closed-form cross-check)" begin
        k_t, k_n = 1.0, 4.0
        K₀ = TensND.Tens(Matrix(Diagonal([k_t, k_t, k_n])))
        crack = PennyCrack(1.0)
        b = cod_tensor(crack, K₀)
        @test b ≈ 2 / (π^2 * sqrt(k_t * k_n)) rtol = 1.0e-13
        R = compliance_contribution(crack, K₀)
        @test R[3, 3] ≈ 3 / (2 * π^2 * sqrt(k_t * k_n)) rtol = 1.0e-13
    end

    @testset "Aniso matrix, cross-check with Hill limit" begin
        # Full anisotropic K₀ (n̂ NOT an eigenvector)
        K_arr = [3.0 0.5 0.3; 0.5 2.0 0.2; 0.3 0.2 1.5]
        K₀ = TensND.Tens(K_arr)

        crack = PennyCrack(1.0)
        R = compliance_contribution(crack, K₀)
        ΔR = delta_resistivity(crack, R, 1.0)
        R_num = hill_limit_R(crack, K₀; ε_small = 1.0e-9)
        for i in 1:3, j in 1:3
            @test ΔR[i, j] ≈ R_num[i, j] rtol = 1.0e-6 atol = 1.0e-8
        end
    end

    @testset "Aniso ribbon, aligned transverse block" begin
        # Ribbon with diagonal K₀: K_⊥ = diag(k_2, k_3).
        K₀ = TensND.Tens(Matrix(Diagonal([3.0, 4.0, 1.5])))
        crack = RibbonCrack(0.8)
        b = cod_tensor(crack, K₀)
        @test b ≈ 2 / (π * sqrt(4.0 * 1.5)) rtol = 1.0e-13
    end

    @testset "Rotated crack basis — invariance of R magnitude" begin
        k₀ = 1.0
        K₀ = TensISO{3}(k₀)
        # Canonical and rotated pennies should give the same scalar b
        b_canonical = cod_tensor(PennyCrack(1.0), K₀)
        rot = TensND.RotatedBasis(0.3, 0.4, 0.2)
        b_rotated = cod_tensor(EllipticCrack(1.0, 1.0, rot), K₀)
        @test b_canonical ≈ b_rotated rtol = 1.0e-13
    end

end

@testset "Cracks / conductivity — SIF and DIF" begin

    @testset "Ribbon thermal SIF" begin
        K₀ = TensISO{3}(2.5)
        b_ribbon = 0.8
        crack = RibbonCrack(b_ribbon)
        q∞ = TensND.Tens([1.0, 0.0, 3.0])
        KT = sif(crack, K₀, q∞)
        @test KT ≈ sqrt(π * b_ribbon) * 3.0 rtol = 1.0e-13
    end

    @testset "Penny thermal DIF: [T] = b · (n̂·q∞)" begin
        k₀ = 1.5
        K₀ = TensISO{3}(k₀)
        crack = PennyCrack(1.0)
        b = cod_tensor(crack, K₀)
        q∞ = TensND.Tens([2.0, 1.0, 4.0])
        T_ = dif(crack, K₀, q∞)
        @test T_ ≈ b * 4.0 rtol = 1.0e-13
    end

end
