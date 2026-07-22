# =============================================================================
#  test_rotational_average.jl — exact rotation-group averages (Core).
#
#  Oracles :
#  1. echoes C++ closed form (`tensor_symmetry.h::transverse_isotropify_
#     around_ez`), transcribed verbatim — the reference implementation the
#     Julia port must match coefficient by coefficient.  echoes' Kelvin-
#     Mandel index convention (4↔23, 5↔31, 6↔12) matches Tensors.jl's.
#  2. Discrete azimuthal quadrature: the average of a 4th-order tensor over
#     N ≥ 5 equally spaced angles about the axis equals the continuous
#     average exactly (the integrand is a trig polynomial of degree ≤ 4).
# =============================================================================

using Test
using MeanFieldHom
using TensND
using LinearAlgebra
using Random
using ForwardDiff
const MC = MeanFieldHom.Core

# ── helpers ─────────────────────────────────────────────────────────────────

function _rand_minor4(rng)
    a = randn(rng, 3, 3, 3, 3)
    b = zeros(3, 3, 3, 3)
    for i in 1:3, j in 1:3, k in 1:3, l in 1:3
        b[i, j, k, l] = (a[i, j, k, l] + a[j, i, k, l] + a[i, j, l, k] + a[j, i, l, k]) / 4
    end
    return b
end

function _rot_axis(n, φ)
    nv = collect(n) ./ norm(collect(n))
    K = [0 -nv[3] nv[2]; nv[3] 0 -nv[1]; -nv[2] nv[1] 0]
    return I(3) .+ sin(φ) .* K .+ (1 - cos(φ)) .* (K * K)
end

function _rotate4(arr, R)
    out = zeros(3, 3, 3, 3)
    for i in 1:3, j in 1:3, k in 1:3, l in 1:3
        s = 0.0
        for a in 1:3, b in 1:3, c in 1:3, d in 1:3
            s += R[i, a] * R[j, b] * R[k, c] * R[l, d] * arr[a, b, c, d]
        end
        out[i, j, k, l] = s
    end
    return out
end

_ti_oracle(arr, n; N = 16) =
    sum(_rotate4(arr, _rot_axis(n, 2π * k / N)) for k in 0:(N - 1)) ./ N

# echoes `transverse_isotropify_around_ez` — verbatim transcription
function _echoes_ti_ez(M)
    t4 = M[6, 6] / 4
    t6 = M[1, 2] / 8 + M[2, 1] / 8 + 3M[1, 1] / 8 + t4 + 3M[2, 2] / 8
    t11 = 3M[1, 2] / 8 + 3M[2, 1] / 8 + M[1, 1] / 8 - t4 + M[2, 2] / 8
    t12 = M[1, 3] + M[2, 3]
    t13 = -M[6, 1] + M[6, 2] + M[1, 6] - M[2, 6]
    t14 = M[3, 1] + M[3, 2]
    t15 = M[4, 4] + M[5, 5]
    t16 = M[4, 5] - M[5, 4]
    TI = zeros(6, 6)
    TI[1, 1] = t6; TI[2, 2] = t6
    TI[1, 2] = t11; TI[2, 1] = t11
    TI[1, 3] = t12 / 2; TI[2, 3] = t12 / 2
    TI[1, 6] = t13 / 4; TI[2, 6] = -t13 / 4
    TI[3, 1] = t14 / 2; TI[3, 2] = t14 / 2
    TI[3, 3] = M[3, 3]
    TI[4, 4] = t15 / 2; TI[5, 5] = t15 / 2
    TI[4, 5] = t16 / 2; TI[5, 4] = -t16 / 2
    TI[6, 1] = -t13 / 4; TI[6, 2] = t13 / 4
    TI[6, 6] = -M[2, 1] / 4 - M[1, 2] / 4 + M[6, 6] / 2 + M[1, 1] / 4 + M[2, 2] / 4
    return TI
end

@testset "rotational averages (Core)" begin
    rng = MersenneTwister(20260722)

    @testset "TI-4 vs echoes closed form (ez), random non-symmetric inputs" begin
        for _ in 1:20
            arr = _rand_minor4(rng)
            M = MC.mandel66_minor(arr)
            avg = MC.transverse_isotropify(TensND.Tens(arr), (0.0, 0.0, 1.0))
            @test avg isa TensND.TensTI{4, Float64, 8}
            TI = MC.mandel66_minor(TensND.get_array(avg))
            @test maximum(abs, TI .- _echoes_ti_ez(M)) < 1.0e-12
        end
    end

    @testset "TI-4 vs discrete quadrature oracle (arbitrary axis)" begin
        nax = (0.36, -0.48, 0.8)
        for _ in 1:5
            arr = _rand_minor4(rng)
            avg = MC.transverse_isotropify(TensND.Tens(arr), nax)
            @test maximum(abs, TensND.get_array(avg) .- _ti_oracle(arr, nax)) < 1.0e-12
        end
    end

    @testset "invariances" begin
        nax = (0.36, -0.48, 0.8)
        arr = _rand_minor4(rng)
        t = TensND.Tens(arr)
        avg = MC.transverse_isotropify(t, nax)
        # idempotence
        avg2 = MC.transverse_isotropify(avg, nax)
        @test maximum(abs, TensND.get_array(avg2) .- TensND.get_array(avg)) < 1.0e-12
        # ISO ∘ TI == ISO
        @test maximum(
            abs,
            TensND.get_array(MC.isotropify(avg)) .- TensND.get_array(MC.isotropify(t))
        ) < 1.0e-12
        # transposition equivariance (catches ℓ₇/ℓ₈ sign errors)
        arrT = permutedims(arr, (3, 4, 1, 2))
        avgT = MC.transverse_isotropify(TensND.Tens(arrT), nax)
        @test maximum(
            abs,
            TensND.get_array(avgT) .- permutedims(TensND.get_array(avg), (3, 4, 1, 2))
        ) < 1.0e-12
        # commutation with a rotation about the axis
        R = _rot_axis(nax, 0.7)
        avgR = MC.transverse_isotropify(TensND.Tens(_rotate4(arr, R)), nax)
        @test maximum(abs, TensND.get_array(avgR) .- TensND.get_array(avg)) < 1.0e-12
    end

    @testset "ISO average == echoes isotropify (α, β)" begin
        arr = _rand_minor4(rng)
        M = MC.mandel66_minor(arr)
        iso = MC.isotropify(TensND.Tens(arr))
        α, β = MC.iso_average_mandel66(M)
        d = TensND.get_data(iso)
        @test d[1] ≈ α
        @test d[2] ≈ β
    end

    @testset "TI-2 : closed form, antisymmetric part, oracle" begin
        nax = (0.36, -0.48, 0.8)
        m2 = randn(rng, 3, 3)
        avg = MC.transverse_isotropify(TensND.Tens(m2), nax)
        @test avg isa TensND.TensTI{2, Float64, 3}
        acc = sum((R -> R * m2 * R')(_rot_axis(nax, 2π * k / 8)) for k in 0:7) ./ 8
        @test maximum(abs, TensND.get_array(avg) .- acc) < 1.0e-12
        # echoes 2nd-order ez closed form
        avgz = MC.transverse_isotropify(TensND.Tens(m2), (0.0, 0.0, 1.0))
        a = TensND.get_array(avgz)
        @test a[1, 1] ≈ (m2[1, 1] + m2[2, 2]) / 2
        @test a[1, 2] ≈ (m2[1, 2] - m2[2, 1]) / 2
        @test a[2, 1] ≈ -(m2[1, 2] - m2[2, 1]) / 2
        @test a[3, 3] ≈ m2[3, 3]
    end

    @testset "ti_average_mandel66 (ALV block form)" begin
        nax = (0.36, -0.48, 0.8)
        arr = _rand_minor4(rng)
        M = MC.mandel66_minor(arr)
        @test maximum(
            abs,
            MC.ti_average_mandel66(M, nax) .- MC.mandel66_minor(_ti_oracle(arr, nax))
        ) < 1.0e-12
    end

    @testset "ForwardDiff through the averages" begin
        nax = (0.36, -0.48, 0.8)
        arr = _rand_minor4(rng)
        f = x -> begin
            ar = arr .* x
            TensND.get_data(MC.transverse_isotropify(TensND.Tens(ar), nax))[7]
        end
        g = ForwardDiff.derivative(f, 2.0)
        h = 1.0e-6
        @test g ≈ (f(2.0 + h) - f(2.0 - h)) / 2h atol = 1.0e-6
    end
end
