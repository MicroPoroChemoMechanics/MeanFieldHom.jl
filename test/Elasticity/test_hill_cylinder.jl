using Test
using MeanFieldHom
using TensND
using ForwardDiff

@testset "Elasticity — Hill tensor for infinite cylinders" begin

    K, μ = 175.0, 80.77
    C_iso = TensISO{3}(3K, 2μ)

    @testset "Circular cylinder, iso matrix → TensTI{4}" begin
        cyl = Cylinder(1.5)
        P = hill_tensor(cyl, C_iso)
        @test P isa TensND.TensTI{4}
        # No coupling with axial direction
        @test P[1, 1, 1, 1] ≈ 0.0 atol = 1.0e-14
        @test P[1, 1, 2, 2] ≈ 0.0 atol = 1.0e-14
        @test P[1, 1, 3, 3] ≈ 0.0 atol = 1.0e-14
        # Transverse isotropy : P_{2222} = P_{3333}, P_{1212} = P_{1313}
        @test P[2, 2, 2, 2] ≈ P[3, 3, 3, 3] atol = 1.0e-14
        @test P[1, 2, 1, 2] ≈ P[1, 3, 1, 3] atol = 1.0e-14
        # Limit consistency with a very elongated prolate
        ell_prolate = Ellipsoid(1.0e8, 1.5, 1.5)
        Pref = hill_tensor(ell_prolate, C_iso)
        for (i, j, k, l) in [(2, 2, 2, 2), (3, 3, 3, 3), (2, 2, 3, 3),
                (1, 2, 1, 2), (2, 3, 2, 3)]
            @test P[i, j, k, l] ≈ Pref[i, j, k, l] rtol = 1.0e-6
        end
    end

    @testset "Elliptic cylinder, iso matrix → TensOrtho" begin
        cyl = Cylinder(2.0, 1.0)
        P = hill_tensor(cyl, C_iso)
        @test P isa TensND.TensOrtho
        @test P[1, 1, 1, 1] ≈ 0.0 atol = 1.0e-14
        @test P[1, 1, 2, 2] ≈ 0.0 atol = 1.0e-14
        @test P[1, 1, 3, 3] ≈ 0.0 atol = 1.0e-14
        # Limit consistency
        ell_big = Ellipsoid(1.0e8, 2.0, 1.0)
        Pref = hill_tensor(ell_big, C_iso)
        for (i, j, k, l) in [(2, 2, 2, 2), (3, 3, 3, 3), (2, 2, 3, 3),
                (1, 2, 1, 2), (1, 3, 1, 3), (2, 3, 2, 3)]
            @test P[i, j, k, l] ≈ Pref[i, j, k, l] rtol = 1.0e-10
        end
    end

    @testset "b → c continuity (elliptic → circular)" begin
        # As c → b, the elliptic-cylinder tensor must converge to the
        # circular-cylinder tensor (same physical limit, different dispatch
        # branch).
        b = 1.2
        P_circ = hill_tensor(Cylinder(b), C_iso)
        for ε in (1.0e-3, 1.0e-5, 1.0e-8)
            P_ell = hill_tensor(Cylinder(b, b * (1 - ε)), C_iso)
            @test P_ell[2, 2, 2, 2] ≈ P_circ[2, 2, 2, 2] rtol = max(1.0e-2, 10ε)
            @test P_ell[2, 3, 2, 3] ≈ P_circ[2, 3, 2, 3] rtol = max(1.0e-2, 10ε)
        end
    end

    @testset "Elliptic cylinder, anisotropic matrix — quadrature path" begin
        C_ortho = TensND.TensOrtho(
            210.0, 200.0, 150.0, 80.0, 70.0, 60.0, 90.0, 85.0, 75.0,
            TensND.CanonicalBasis{3, Float64}()
        )
        cyl = Cylinder(2.0, 1.0)
        P = hill_tensor(cyl, C_ortho)
        @test P[1, 1, 1, 1] ≈ 0.0 atol = 1.0e-12
        # Limit consistency with the nested-QuadGK path on a very
        # elongated triaxial (tight 1e-8 comparison).  Real DECUHR at
        # default reltol=1e-6 would not meet this bound.
        ell_big = Ellipsoid(1.0e6, 2.0, 1.0)
        Pref = hill_tensor(ell_big, C_ortho; method = :nestedquadgk)
        for (i, j, k, l) in [(2, 2, 2, 2), (3, 3, 3, 3), (2, 3, 2, 3),
                (1, 2, 1, 2), (1, 3, 1, 3)]
            @test P[i, j, k, l] ≈ Pref[i, j, k, l] rtol = 1.0e-8
        end
    end

    @testset "Ellipsoid → Cylinder / Crack redirection" begin
        @test Ellipsoid(Inf, 2.0, 1.0) isa Cylinder{EllipticCylindrical}
        @test Ellipsoid(1.0, 1.0, Inf) isa Cylinder{CircularCylindrical}
        @test Ellipsoid(Inf, 1.5, 1.5) isa Cylinder{CircularCylindrical}
        @test Ellipsoid(2.0, 1.0, 0.0) isa MeanFieldHom.Cracks.EllipticCrack
        @test Ellipsoid(1.0, 1.0, 0.0) isa MeanFieldHom.Cracks.EllipticCrack
        @test Ellipsoid(Inf, 1.0, 0.0) isa MeanFieldHom.Cracks.RibbonCrack
        @test_throws ArgumentError Ellipsoid(Inf, Inf, 1.0)
        @test_throws ArgumentError Ellipsoid(2.0, 0.0, 0.0)
        @test_throws ArgumentError Ellipsoid(0.0, 0.0, 0.0)
    end

    @testset "ForwardDiff through cylinder analytical path" begin
        f(b) = hill_tensor(Cylinder(b), C_iso)[2, 2, 2, 2]
        df = ForwardDiff.derivative(f, 1.0)
        @test isfinite(df)
        @test !isnan(df)
        g(b) = hill_tensor(Cylinder(b, b / 2), C_iso)[2, 2, 2, 2]
        dg = ForwardDiff.derivative(g, 2.0)
        @test isfinite(dg)
        @test !isnan(dg)
    end

    @testset "Auxiliary tensors on cylinder" begin
        cyl = Cylinder(2.0, 1.0)
        IA = tens_IA(cyl)
        UA = tens_UA(cyl)
        VA = tens_VA(cyl)
        # ΣIᵢ = 1, axial zero
        @test IA[1, 1] ≈ 0.0 atol = 1.0e-14
        @test IA[1, 1] + IA[2, 2] + IA[3, 3] ≈ 1.0 atol = 1.0e-12
        # V^A diagonal equals I^A on the diagonal (by construction)
        @test VA[1, 1, 1, 1] ≈ IA[1, 1] atol = 1.0e-12
        @test VA[2, 2, 2, 2] ≈ IA[2, 2] atol = 1.0e-12
        @test VA[3, 3, 3, 3] ≈ IA[3, 3] atol = 1.0e-12
        # U^A : axial components zero, transverse finite
        @test UA[1, 1, 1, 1] ≈ 0.0 atol = 1.0e-14
        @test UA[1, 1, 2, 2] ≈ 0.0 atol = 1.0e-14

        # Circular case — TensTI{4} structure
        cyl2 = Cylinder(1.3)
        UA2 = tens_UA(cyl2)
        @test UA2 isa TensND.TensTI{4}
    end
end
