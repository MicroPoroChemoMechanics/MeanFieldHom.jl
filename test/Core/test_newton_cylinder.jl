using Test
using MeanFieldHom
using ForwardDiff

@testset "Core — newton_potential_3d_cylinder" begin
    # Circular base (b = c): I_b = I_c = 2π, I_bb = I_cc = I_bc = π/b²
    b = 1.5
    Iv, IIv = MeanFieldHom.Core.newton_potential_3d_cylinder(b, b)
    @test Iv[1] == 0.0
    @test Iv[2] ≈ 2π atol = 1.0e-12
    @test Iv[3] ≈ 2π atol = 1.0e-12
    @test sum(Iv) ≈ 4π atol = 1.0e-12
    @test IIv[1] == 0.0
    @test IIv[5] == 0.0
    @test IIv[6] == 0.0
    @test IIv[2] ≈ π / b^2 atol = 1.0e-12
    @test IIv[3] ≈ π / b^2 atol = 1.0e-12
    @test IIv[4] ≈ π / b^2 atol = 1.0e-12

    # Elliptic base (b > c): I_b = 4π c/(b+c), I_c = 4π b/(b+c), I_bc = 4π/(b+c)²
    b2, c2 = 2.0, 1.0
    Iv, IIv = MeanFieldHom.Core.newton_potential_3d_cylinder(b2, c2)
    s = b2 + c2
    @test Iv[1] == 0.0
    @test Iv[2] ≈ 4π * c2 / s atol = 1.0e-12
    @test Iv[3] ≈ 4π * b2 / s atol = 1.0e-12
    @test sum(Iv) ≈ 4π atol = 1.0e-12
    @test IIv[1] == 0.0
    @test IIv[4] ≈ 4π / s^2 atol = 1.0e-12
    @test IIv[2] ≈ 4π / 3 * (1 / b2^2 - 1 / s^2) atol = 1.0e-12
    @test IIv[3] ≈ 4π / 3 * (1 / c2^2 - 1 / s^2) atol = 1.0e-12

    # Numerical consistency with the general 3D Newton potential in the limit
    # a → ∞ — use a = 1e6·max(b,c) as a proxy (larger a triggers numerical
    # overflow in the elliptic-integral formula for a triaxial ellipsoid).
    Iv_inf, _ = MeanFieldHom.Core.newton_potential_3d(1.0e6, b2, c2)
    @test Iv_inf[1] ≈ 0.0 atol = 1.0e-3
    @test Iv_inf[2] ≈ 4π * c2 / s rtol = 1.0e-4
    @test Iv_inf[3] ≈ 4π * b2 / s rtol = 1.0e-4

    # BigFloat — high-precision evaluation
    bb = big"1.5"
    cc = big"0.7"
    Iv_big, _ = MeanFieldHom.Core.newton_potential_3d_cylinder(bb, cc)
    @test sum(Iv_big) ≈ 4 * big(π) atol = big"1e-50"

    # ForwardDiff — derivative w.r.t. transverse semi-axis
    f(x) = MeanFieldHom.Core.newton_potential_3d_cylinder(x, 1.0)[1][2]  # I_b
    df = ForwardDiff.derivative(f, 2.0)
    @test isfinite(df)
    @test !isnan(df)
end
