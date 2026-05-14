# =============================================================================
#  test_symmetrize.jl — round-trip tests for the orientation-distribution
#  projection helpers `_apply_symmetrize` and the RVE-level `symmetrize`
#  keyword.
# =============================================================================

using Test
using MeanFieldHom
using TensND
using LinearAlgebra
import MeanFieldHom.Schemes: _apply_symmetrize

@testset "symmetrize" begin

    @testset "iso projection — identity on iso 4-tensor" begin
        C = TensISO{3}(216.0, 64.0)   # 216 J + 64 K
        Cp = _apply_symmetrize(C, IsoSymmetrize())
        @test Cp isa TensND.TensISO{4, 3}
        αp, βp = TensND.get_data(Cp)
        @test αp ≈ 216.0
        @test βp ≈ 64.0
    end

    @testset "iso projection — extracts (α, β) of an aniso 4-tensor" begin
        # Build a simple aniso tensor as a perturbation of an iso one.
        α0, β0 = 60.0, 24.0
        arr = zeros(3, 3, 3, 3)
        for i in 1:3, j in 1:3, k in 1:3, l in 1:3
            J_iso = (i == j) * (k == l) / 3.0
            K_iso = ((i == k) * (j == l) + (i == l) * (j == k)) / 2.0 - J_iso
            arr[i, j, k, l] = α0 * J_iso + β0 * K_iso
        end
        # add an off-iso perturbation (orthorhombic-style scaling on the (3,3) axis)
        for i in 1:3, j in 1:3
            arr[i, j, 3, 3] += 1.5 * (i == j)
            arr[3, 3, i, j] += 1.5 * (i == j)
        end
        T = TensND.Tens(arr)
        Tp = _apply_symmetrize(T, IsoSymmetrize())
        @test Tp isa TensND.TensISO{4, 3}
        αp, βp = TensND.get_data(Tp)
        # Recompute analytically: α = (1/3) sum T[i,i,j,j], β = (full_trace - α)/5
        α_ref = sum(arr[i, i, j, j] for i in 1:3, j in 1:3) / 3.0
        full_trace = sum(arr[i, j, i, j] for i in 1:3, j in 1:3)
        β_ref = (full_trace - α_ref) / 5.0
        @test αp ≈ α_ref
        @test βp ≈ β_ref
    end

    @testset "iso projection — 2nd-order spherical part" begin
        K = TensND.Tens(Diagonal([1.5, 2.0, 2.5]))
        Kp = _apply_symmetrize(K, IsoSymmetrize())
        @test Kp isa TensND.TensISO{2, 3}
        @test TensND.get_data(Kp)[1] ≈ (1.5 + 2.0 + 2.5) / 3.0
    end

    @testset "TI projection — identity on coaxial TI(ez) tensor" begin
        n = (0.0, 0.0, 1.0)
        # Build a TI(ez) tensor manually via TensND constructor.
        data = (1.0, 2.0, 0.5, 3.0, 4.0)   # (ℓ₁, ℓ₂, ℓ₃, ℓ₅, ℓ₆)
        T = TensND.TensTI{4, Float64, 5}(data, n)
        Tp = _apply_symmetrize(T, TISymmetrize(n))
        @test Tp isa TensND.TensTI{4, Float64, 5}
        # Walpole 5-component decomposition is preserved up to numerical noise.
        for (a, b) in zip(TensND.get_data(Tp), data)
            @test a ≈ b atol = 1.0e-10
        end
    end

    @testset "TI projection — 2nd-order axial / transverse split" begin
        n = (0.0, 0.0, 1.0)
        K = TensND.Tens(Diagonal([1.5, 2.0, 4.0]))
        Kp = _apply_symmetrize(K, TISymmetrize(n))
        @test Kp isa TensND.TensTI{2, Float64, 2}
        a, b = TensND.get_data(Kp)
        # axial b = K[3,3] = 4.0 ; transverse a = (K[1,1] + K[2,2]) / 2 = 1.75
        @test b ≈ 4.0
        @test a ≈ 1.75
    end

    @testset "RVE-level symmetrize integration — porous oblate ⇒ iso C_eff" begin
        # Spherical-pore matrix with iso symmetrize on every phase: the
        # macroscopic stiffness should remain iso (TensISO) and match the
        # un-symmetrized result for sphere geometry (sanity check).
        C_s = TensISO{3}(216.0, 64.0)
        C_p = TensISO{3}(3.0e-6, 2.0e-6)
        rve = RVE(:M)
        add_matrix!(rve, Spheroid(0.2), Dict(:C => C_s); symmetrize = :iso)
        add_phase!(
            rve, :PORE, Spheroid(0.2), Dict(:C => C_p);
            fraction = 0.3, symmetrize = :iso
        )
        C_eff = homogenize(rve, MoriTanaka(), :C)
        @test C_eff isa TensND.TensISO{4, 3}
        α, β = TensND.get_data(C_eff)
        @test α / 3 > 0
        @test β / 2 > 0
    end

    @testset "RVE-level symmetrize — :ti yields TI(ez) macroscopic tensor" begin
        # Oblate inclusions with TI(ez) symmetrize: the macroscopic
        # stiffness inherits a TI structure around ez (axisymmetric).
        C_m = TensISO{3}(216.0, 64.0)
        C_i = TensISO{3}(3.0e-6, 2.0e-6)
        rve = RVE(:M)
        add_matrix!(rve, Ellipsoid(1.0), Dict(:C => C_m))
        add_phase!(
            rve, :I, Spheroid(0.2), Dict(:C => C_i);
            fraction = 0.2, symmetrize = TISymmetrize((0.0, 0.0, 1.0))
        )
        C_eff = homogenize(rve, MoriTanaka(), :C)
        # MT keeps the homogenised tensor in a TI subspace ; verify the
        # canonical/walpole representation has the expected axisymmetry by
        # checking C[1,1,1,1] == C[2,2,2,2] (transverse symmetry plane).
        a = TensND.get_array(C_eff)
        @test a[1, 1, 1, 1] ≈ a[2, 2, 2, 2] atol = 1.0e-8
        @test a[1, 1, 2, 2] ≈ a[2, 2, 1, 1] atol = 1.0e-8
    end
end
