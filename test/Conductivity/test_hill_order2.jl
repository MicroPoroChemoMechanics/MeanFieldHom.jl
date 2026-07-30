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

# ── Orientation ──────────────────────────────────────────────────────────────
#
# Regression guard. The isotropic 3-D kernels read the Newton-potential
# components in the *inclusion's own* basis; wrapping them straight into a
# canonical `Tens` used to drop the orientation, so `hill_tensor` of a rotated
# spheroid returned the tensor of the unrotated one — exactly, hence silently.
#
# The internal consistency check `s = P·K₀` in `test_eshelby.jl` cannot catch
# this: both sides drop the rotation identically. What catches it is
# equivariance, and a comparison against the anisotropic kernel, which always
# built the shape tensor with the rotation in it.

_rot(basis) = Matrix(TensND.vecbasis(basis, :cov))
_canon2(t) = Matrix(TensND.components_canon(t))

@testset "Conductivity — the order-2 Hill tensor rotates with the inclusion" begin
    K = TensISO{3}(1.7)
    ang = (0.3, 0.7, 0.2)

    @testset "ellipsoid" begin
        for axes in ((1.0, 1.0, 0.4), (2.5, 1.0, 1.0), (3.0, 2.0, 1.0))
            e0 = Ellipsoid(axes...)
            e1 = Ellipsoid(axes...; euler_angles = ang)
            R = _rot(MeanFieldHom.inclusion_basis(e1))
            P0, P1 = _canon2(hill_tensor(e0, K)), _canon2(hill_tensor(e1, K))
            @test P1 ≈ R * P0 * transpose(R) atol = 1.0e-12
            # And it really is rotated: a non-spherical shape must *change*.
            @test !isapprox(P1, P0; atol = 1.0e-6)
        end
    end

    @testset "cylinder" begin
        c0 = Cylinder(1.0, 0.4)
        c1 = Cylinder(1.0, 0.4; euler_angles = ang)
        R = _rot(MeanFieldHom.inclusion_basis(c1))
        P0, P1 = _canon2(hill_tensor(c0, K)), _canon2(hill_tensor(c1, K))
        @test P1 ≈ R * P0 * transpose(R) atol = 1.0e-12
        @test !isapprox(P1, P0; atol = 1.0e-6)
    end

    @testset "the isotropic kernel agrees with the anisotropic one" begin
        # The sharpest guard available: the anisotropic path is an independent
        # derivation (Giraud's K^(-1/2) change of variable) and was never
        # affected. Reaching for the private kernel is deliberate — the public
        # entry point dispatches `TensISO` to the analytic branch, so there is no
        # other way to compare the two on the same input.
        C = MeanFieldHom.Conductivity
        K = TensISO{3}(1.7)
        for axes in ((1.0, 1.0, 0.4), (3.0, 2.0, 1.0))
            for angles in ((), (0.3, 0.7, 0.2))
                e = Ellipsoid(axes...; euler_angles = angles)
                Pi = _canon2(C._hill_order2_3d_iso(e, K))
                Pa = _canon2(C._hill_order2_3d_aniso(e, K))
                @test Pi ≈ Pa atol = 1.0e-10
            end
        end
    end
end
