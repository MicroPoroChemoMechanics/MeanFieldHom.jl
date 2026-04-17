using Test
using MeanFieldHom
using TensND

@testset "Elasticity — hill_tensor smoke tests" begin
    E, ν = 210.0, 0.3
    k = E / (3 * (1 - 2ν))
    μ = E / (2 * (1 + ν))
    C = TensISO{3}(3k, 2μ)

    # Sphere, ISO
    ell = Ellipsoid(1.0)
    P = hill_tensor(ell, C)
    @test P[1,1,1,1] ≈ P[2,2,2,2] atol=1e-12

    # Auxiliary tensors
    IA = tens_IA(ell)
    @test sum([IA[i,i] for i in 1:3]) ≈ 1.0 atol=1e-12
    UA = tens_UA(ell)
    VA = tens_VA(ell)
    @test UA[1,1,1,1] > 0
    @test VA[1,1,1,1] > 0

    # 2D ellipse
    ell2 = Ellipsoid(1.0, 0.5)
    C2 = TensISO{2}(3k, 2μ)
    P2 = hill_tensor(ell2, C2)
    @test P2[1,1,1,1] > 0
end
