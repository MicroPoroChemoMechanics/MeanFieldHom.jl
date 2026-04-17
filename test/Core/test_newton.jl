using Test
using MeanFieldHom

@testset "Core — newton potentials" begin
    # Sphere: ΣIᵢ = 4π
    Iv, IIv = MeanFieldHom.Core.newton_potential_3d(1.0, 1.0, 1.0)
    @test sum(Iv) ≈ 4π atol = 1.0e-12

    # 2D circle
    Ia, Ib = MeanFieldHom.Core.newton_potential_2d(1.0, 1.0)
    @test Ia ≈ π atol = 1.0e-12
    @test Ib ≈ π atol = 1.0e-12

    # Moduli extractors
    C_iso = TensND.TensISO{3}(3.0, 2.0)
    E, ν = MeanFieldHom.Core.extract_iso_moduli(C_iso)
    @test E > 0
    @test -1 < ν < 0.5
end
