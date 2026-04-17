using Test
using MeanFieldHom
using TensND

@testset "Conductivity — order-2 Hill tensor" begin
    K = TensISO{3}(5.0)
    ell = Ellipsoid(1.0)
    P = hill_tensor(ell, K)
    @test P[1, 1] ≈ 1 / (3 * 5.0) atol = 1.0e-12

    # 2D circle
    K2 = TensISO{2}(3.0)
    ell2 = Ellipsoid(1.0; dim = 2)
    P2 = hill_tensor(ell2, K2)
    @test P2[1, 1] ≈ 1 / (2 * 3.0) atol = 1.0e-12
end
