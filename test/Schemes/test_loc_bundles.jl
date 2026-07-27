# =============================================================================
#  test_loc_bundles.jl — contract of the bundled localization helpers.
#
#  `Core.loc_and_stiffness` / `Core.loc_and_stress_average` (and the crack
#  analog `Cracks.compliance_and_stiffness_contribution`) exist purely to share
#  the single expensive Hill / COD / recurrence solve between two quantities
#  that Mori-Tanaka and the self-consistent kernels always request together.
#
#  Their contract is therefore stronger than "close enough": the bundle must
#  return values **bitwise identical** to the two separate calls.  These tests
#  assert exactly that with `==` on the flattened component arrays — not `≈` —
#  because that identity is what lets the optimization be verified rather than
#  merely believed.
#
#  Also covers the two branches that had no test coverage before the bundles
#  were introduced:
#    * Mori-Tanaka with a `CrackDensity` phase (`_mt_4`/`_mt_2` `else` branch);
#    * a non-`TensISO` phase property under a non-trivial `symmetrize`
#      (the slow branch of `_phase_stress_strain_average`).
# =============================================================================

using Test
using MeanFieldHom
using TensND
using LinearAlgebra
import ForwardDiff as FD

const MFHC_B = MeanFieldHom.Core

_vals(t) = vec(collect(TensND.get_array(t)))

# Triclinic reference — same KM as test/regression/test_anisotropic.jl.
const _KM_TRI_B = [
    210.0 80.0 75.0 5.0 4.0 3.0;
    80.0 195.0 90.0 -2.0 3.0 -1.0;
    75.0 90.0 220.0 1.0 -2.0 2.0;
    5.0 -2.0 1.0 60.0 2.5 1.5;
    4.0 3.0 -2.0 2.5 65.0 -1.0;
    3.0 -1.0 2.0 1.5 -1.0 55.0
]

@testset "loc_and_stiffness / loc_and_stress_average — bitwise contract" begin
    C_tri = TensND.inv_KM(_KM_TRI_B, CanonicalBasis{3, Float64}())
    C_iso = TensISO{3}(3 * 20.0, 2 * 12.0)
    C_i = TensISO{3}(3 * 80.0, 2 * 50.0)
    C_ti = TensTI{4}(20.0, 30.0, 4.0, 5.0, 8.0, (0.0, 0.0, 1.0))

    fixtures4 = [
        ("sphere/iso", Ellipsoid(1.0), C_i, C_iso),
        ("triaxial/iso", Ellipsoid(3.0, 2.0, 1.0), C_i, C_iso),
        ("oblate/iso", Spheroid(0.2), C_i, C_iso),
        ("triaxial/tri", Ellipsoid(3.0, 2.0, 1.0), C_i, C_tri),
        ("sphere/ti_phase", Ellipsoid(1.0), C_ti, C_iso),
        ("cylinder/iso", Cylinder(1.0), C_i, C_iso),
    ]

    @testset "4th order — $nm" for (nm, geom, C₁, C₀) in fixtures4
        A, N = MFHC_B.loc_and_stiffness(geom, C₁, C₀)
        @test _vals(A) == _vals(strain_strain_loc(geom, C₁, C₀))
        @test _vals(N) == _vals(stiffness_contribution(geom, C₁, C₀))

        A2, B2 = MFHC_B.loc_and_stress_average(geom, C₁, C₀)
        @test _vals(A2) == _vals(strain_strain_loc(geom, C₁, C₀))
        @test _vals(B2) == _vals(stress_strain_loc(geom, C₁, C₀))
    end

    K_iso = TensISO{3}(2.0)
    K_i = TensISO{3}(20.0)
    fixtures2 = [
        ("sphere/iso", Ellipsoid(1.0), K_i, K_iso),
        ("oblate/iso", Spheroid(0.2), K_i, K_iso),
    ]

    @testset "2nd order — $nm" for (nm, geom, K₁, K₀) in fixtures2
        A, N = MFHC_B.loc_and_stiffness(geom, K₁, K₀)
        @test _vals(A) == _vals(gradient_gradient_loc(geom, K₁, K₀))
        @test _vals(N) == _vals(conductivity_contribution(geom, K₁, K₀))

        A2, B2 = MFHC_B.loc_and_stress_average(geom, K₁, K₀)
        @test _vals(A2) == _vals(gradient_gradient_loc(geom, K₁, K₀))
        @test _vals(B2) == _vals(flux_gradient_loc(geom, K₁, K₀))
    end

    @testset "heterogeneous inclusions keep the layer-by-layer assembly" begin
        sph = LayeredSphere(
            (0.5, 1.0),
            (TensISO{3}(3 * 10.0, 2 * 6.0), TensISO{3}(3 * 40.0, 2 * 25.0))
        )
        A, N = MFHC_B.loc_and_stiffness(sph, C_i, C_iso)
        @test _vals(A) == _vals(strain_strain_loc(sph, C_i, C_iso))
        @test _vals(N) == _vals(stiffness_contribution(sph, C_i, C_iso))

        A2, B2 = MFHC_B.loc_and_stress_average(sph, C_i, C_iso)
        @test _vals(A2) == _vals(strain_strain_loc(sph, C_i, C_iso))
        @test _vals(B2) == _vals(stress_strain_loc(sph, C_i, C_iso))
    end

    @testset "the bundle costs strictly less than the two separate calls" begin
        # `@inferred` is deliberately NOT used here: `hill_tensor` resolves its
        # algorithm at run time (`_resolve_algo(Val(method), incl, C₀)`), so the
        # whole localization chain is inferrably `Any` BY DESIGN — that predates
        # these bundles and is what makes `:auto` work.  What must hold is that
        # bundling does not make things worse: the tuple return must not add
        # allocation on top of the work it saves.
        geom = Ellipsoid(3.0, 2.0, 1.0)
        pair() = (
            strain_strain_loc(geom, C_i, C_tri),
            stiffness_contribution(geom, C_i, C_tri),
        )
        bundle() = MFHC_B.loc_and_stiffness(geom, C_i, C_tri)
        pair(); bundle()                      # warm up in this process
        for _ in 1:2
            pair(); bundle()
        end
        a_pair = minimum(@allocated(pair()) for _ in 1:5)
        a_bundle = minimum(@allocated(bundle()) for _ in 1:5)
        @test a_bundle < a_pair
        # One Hill solve instead of two ⇒ close to half.
        @test a_bundle < 0.6 * a_pair
    end

    @testset "ForwardDiff flows through the bundle" begin
        g = FD.derivative(1.0) do x
            A, N = MFHC_B.loc_and_stiffness(
                Ellipsoid(3.0, 2.0, 1.0), x * C_i, C_iso
            )
            TensND.get_array(N)[1, 1, 1, 1]
        end
        h = 1.0e-6
        ref = (
            TensND.get_array(
                stiffness_contribution(Ellipsoid(3.0, 2.0, 1.0), (1 + h) * C_i, C_iso)
            )[1, 1, 1, 1] -
                TensND.get_array(
                stiffness_contribution(Ellipsoid(3.0, 2.0, 1.0), (1 - h) * C_i, C_iso)
            )[1, 1, 1, 1]
        ) / 2h
        @test g ≈ ref rtol = 1.0e-5
    end
end

@testset "compliance_and_stiffness_contribution — bitwise contract" begin
    C_tri = TensND.inv_KM(_KM_TRI_B, CanonicalBasis{3, Float64}())
    C_iso = TensISO{3}(3 * 24.0, 2 * 12.0)

    @testset "elastic — $nm" for (nm, crack, C₀) in [
            ("penny/iso", PennyCrack(1.0), C_iso),
            ("penny/tri", PennyCrack(1.0), C_tri),
            ("ellipse/iso", EllipticCrack(1.0, 0.4), C_iso),
        ]
        H, N = MeanFieldHom.Cracks.compliance_and_stiffness_contribution(crack, C₀)
        @test _vals(H) == _vals(compliance_contribution(crack, C₀))
        @test _vals(N) == _vals(stiffness_contribution(crack, C₀))
    end

    @testset "conductivity" begin
        K₀ = TensISO{3}(2.0)
        crack = PennyCrack(1.0)
        R, N = MeanFieldHom.Cracks.compliance_and_stiffness_contribution(crack, K₀)
        @test _vals(R) == _vals(compliance_contribution(crack, K₀))
        @test _vals(N) == _vals(conductivity_contribution(crack, K₀))
    end
end

# =============================================================================
#  Branches that had no coverage before the bundles were introduced
# =============================================================================

@testset "MoriTanaka — crack phase (previously untested branch)" begin
    C_m = TensISO{3}(30.0, 10.0)
    C33_m = TensND.get_array(C_m)[3, 3, 3, 3]

    rve = RVE(:M)
    add_matrix!(rve, Ellipsoid(1.0), Dict(:C => C_m))
    add_phase!(rve, :CRACK, PennyCrack(1.0), Dict(:C => C_m); density = 0.1)

    Cmt = homogenize(rve, MoriTanaka())
    @test all(isfinite, TensND.get_array(Cmt))
    @test TensND.get_array(Cmt)[3, 3, 3, 3] < C33_m

    # MT → Dilute as the crack density vanishes.
    function build(ε)
        r = RVE(:M)
        add_matrix!(r, Ellipsoid(1.0), Dict(:C => C_m))
        add_phase!(r, :CRACK, PennyCrack(1.0), Dict(:C => C_m); density = ε)
        return r
    end
    Cd = TensND.get_array(homogenize(build(1.0e-4), Dilute()))
    Cm = TensND.get_array(homogenize(build(1.0e-4), MoriTanaka()))
    @test maximum(abs, Cd .- Cm) < 1.0e-4 * maximum(abs, Cd)

    @testset "conductivity crack phase" begin
        K_m = TensISO{3}(2.0)
        r = RVE(:M)
        add_matrix!(r, Ellipsoid(1.0), Dict(:K => K_m))
        add_phase!(r, :CRACK, PennyCrack(1.0), Dict(:K => K_m); density = 0.1)
        Kmt = homogenize(r, MoriTanaka(); property = :K)
        @test all(isfinite, TensND.get_array(Kmt))
        @test TensND.get_array(Kmt)[3, 3] < TensND.get_array(K_m)[3, 3]
    end
end

@testset "non-isotropic phase property under a non-trivial symmetrize" begin
    # The only configuration exercising the slow branch of
    # `_phase_stress_strain_average` (homogeneous inclusion, `sym` non-trivial,
    # `P_i` NOT a `TensISO`): the symmetrization of the PRODUCT `P_i ⊡ A_raw`
    # differs from `P_i ⊡ symmetrize(A_raw)` when `P_i` does not commute with
    # the rotation group being averaged over.
    C_m = TensISO{3}(3 * 20.0, 2 * 12.0)
    C_ti = TensTI{4}(20.0, 30.0, 4.0, 5.0, 8.0, (0.0, 0.0, 1.0))
    θ = π / 4
    axis = (sin(θ), 0.0, cos(θ))

    rve = RVE(:M)
    add_matrix!(rve, Ellipsoid(1.0), Dict(:C => C_m))
    add_phase!(
        rve, :I, Spheroid(5.0; euler_angles = (θ, 0.0, 0.0)), Dict(:C => C_ti);
        fraction = 0.1, symmetrize = TISymmetrize(axis)
    )

    # Only Mori-Tanaka here. `SelfConsistent` on this configuration hits a
    # PRE-EXISTING gap, unrelated to the bundles and reproducible at the commit
    # before them: the running estimate becomes a `TensTI{4,T,8}` (the exact
    # 8-parameter symmetrization result) and the analytical TI-coaxial Hill
    # kernel only has methods for the 5-/6-parameter `TensTI`, so it raises
    # `MethodError: no method matching _hill_3d_ti_coaxial(::Ellipsoid{Spherical},
    # ::TensTI{4,Float64,8})`.  Marked broken rather than silently dropped.
    C_eff = homogenize(rve, MoriTanaka(), :C)
    @test all(isfinite, TensND.get_array(C_eff))
    @test_broken (homogenize(rve, SelfConsistent(), :C); true)

    # Independent reference for the stress average of that phase, built from
    # the raw localization exactly as the slow branch prescribes.
    sym = MeanFieldHom.Schemes.phase_symmetrize(rve, :I)
    P₀_proj = MeanFieldHom.Schemes._project_matrix(C_m, sym)
    geom = rve.phases[:I].geometry
    A_raw = strain_strain_loc(geom, C_ti, P₀_proj)
    CA_ref = MeanFieldHom.Schemes._apply_symmetrize(C_ti ⊡ A_raw, sym)

    A_dil, CA = MeanFieldHom.Schemes._phase_dilute_and_stress_average(rve, :I, :C, C_m)
    @test _vals(CA) == _vals(CA_ref)
    @test _vals(A_dil) == _vals(MeanFieldHom.Schemes._apply_symmetrize(A_raw, sym))

    # And it really is the non-commuting branch: the shortcut form differs.
    @test _vals(CA) != _vals(C_ti ⊡ A_dil)

    @testset "ForwardDiff through the slow branch" begin
        function build(x)
            r = RVE(:M)
            add_matrix!(r, Ellipsoid(1.0), Dict(:C => TensISO{3}(3 * 20.0 * x, 2 * 12.0 * x)))
            add_phase!(
                r, :I, Spheroid(5.0; euler_angles = (θ, 0.0, 0.0)), Dict(:C => C_ti);
                fraction = 0.1, symmetrize = TISymmetrize(axis)
            )
            return r
        end
        f = x -> TensND.get_array(homogenize(build(x), MoriTanaka(), :C))[3, 3, 3, 3]
        @test isfinite(f(1.0))
        # PRE-EXISTING TensND limitation (reproducible before the bundles):
        # differentiating a scheme whose phase property is a `TensTI` promotes
        # the five Walpole parameters to `Dual` while the AXIS stays `Float64`,
        # and the inner constructor
        # `TensTI{order,T,N}(::NTuple{N,T}, ::Tuple{T,T,T})` demands a single
        # `T` for both.  Fix belongs in TensND, not here.
        @test_broken (FD.derivative(f, 1.0); true)
    end
end
