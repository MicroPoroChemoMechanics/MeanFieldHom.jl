using Test
using MeanFieldHom

# =============================================================================
#  test_orientation.jl — familles d'orientation discrétisées
#  (`src/Schemes/orientation.jl`).
#
#  Convention Pichler-Hellmich (2011) / `disc_theta` d'echoes :
#    θ_i = (π/2)(i-1)/(N-1),  bords aux mi-points bornés à [0, π/2],
#    w_i = cos θ_i⁻ − cos θ_i⁺,  Σ w_i = 1.
# =============================================================================

@testset "polar_orientation_bins — convention et normalisation" begin
    for N in (2, 3, 5, 12, 64)
        bins = polar_orientation_bins(N)
        @test length(bins) == N

        θ = [b.θ for b in bins]
        w = [b.weight for b in bins]

        # Angles : bornes incluses, strictement croissants, équirépartis.
        @test θ[1] ≈ 0.0
        @test θ[end] ≈ π / 2
        @test issorted(θ)
        @test all(θ .≥ 0.0)
        @test all(θ .≤ π / 2 + 1.0e-15)
        for i in 1:N
            @test θ[i] ≈ (π / 2) * (i - 1) / (N - 1)
        end

        # Poids : partition de l'hémisphère, Σ w = 1, tous positifs.
        @test sum(w) ≈ 1.0
        @test all(w .> 0.0)
    end
end

@testset "polar_orientation_bins — N = 2, cas limite" begin
    bins = polar_orientation_bins(2)
    @test length(bins) == 2
    @test bins[1].θ ≈ 0.0
    @test bins[2].θ ≈ π / 2

    # Bords : θ₁⁻ = 0 (borné), θ₁⁺ = π/4 ; θ₂⁻ = π/4, θ₂⁺ = π/2 (borné).
    @test bins[1].weight ≈ 1 - cos(π / 4)
    @test bins[2].weight ≈ cos(π / 4)
    @test bins[1].weight + bins[2].weight ≈ 1.0
end

@testset "polar_orientation_bins — type de retour" begin
    bins = polar_orientation_bins(4)
    @test bins isa Vector{@NamedTuple{θ::Float64, weight::Float64}}
    @test bins[1] isa NamedTuple
    @test propertynames(bins[1]) == (:θ, :weight)
end

@testset "polar_orientation_bins — N < 2 rejeté" begin
    @test_throws ArgumentError polar_orientation_bins(1)
    @test_throws ArgumentError polar_orientation_bins(0)
    @test_throws ArgumentError polar_orientation_bins(-3)
end

@testset "polar_orientation_bins — convergence vers la moyenne isotrope" begin
    # Σ wᵢ f(θᵢ) doit converger vers ∫ f(θ) sin θ dθ sur [0, π/2] en O(Δθ²).
    # On teste sur f(θ) = cos²θ, dont la moyenne exacte est 1/3.
    err(N) = abs(sum(b.weight * cos(b.θ)^2 for b in polar_orientation_bins(N)) - 1 / 3)

    e_small = err(16)
    e_big = err(64)
    @test e_big < e_small
    # Ordre 2 : diviser Δθ par 4 doit diviser l'erreur par ~16.
    @test e_big < e_small / 8
    @test e_big < 1.0e-3
end
