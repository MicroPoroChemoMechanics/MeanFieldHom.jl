# =============================================================================
#  test_symmetrize_alv.jl — block-wise orientation averages of ALV Volterra
#  matrices (`_iso_project_blocks`, `_ti_project_blocks`,
#  `_maybe_symmetrize_alv`) against the elastic Core implementations.
# =============================================================================

using Test
using MeanFieldHom
using TensND
using LinearAlgebra
using Random
import MeanFieldHom.Viscoelasticity: _iso_project_blocks, _ti_project_blocks,
    _maybe_symmetrize_alv, _iso_project_mandel66
const MCr = MeanFieldHom.Core

function _rand_minor_mandel(rng)
    a = randn(rng, 3, 3, 3, 3)
    b = zeros(3, 3, 3, 3)
    for i in 1:3, j in 1:3, k in 1:3, l in 1:3
        b[i, j, k, l] = (a[i, j, k, l] + a[j, i, k, l] + a[i, j, l, k] + a[j, i, l, k]) / 4
    end
    return MCr.mandel66_minor(b)
end

@testset "ALV block-wise orientation averages" begin
    rng = MersenneTwister(42)
    nax = (0.36, -0.48, 0.8)
    ez = (0.0, 0.0, 1.0)

    @testset "single block == elastic Core implementation" begin
        M = _rand_minor_mandel(rng)
        # ISO — non-symmetric block: the full Σ M[i,j] must be used, not
        # a symmetric 2·M[1,2] shortcut (regression for Volterra blocks)
        α, β = _iso_project_mandel66(M)
        αc, βc = MCr.iso_average_mandel66(M)
        @test α ≈ αc
        @test β ≈ βc
        # TI about an arbitrary axis
        @test _ti_project_blocks(M, nax) ≈ MCr.ti_average_mandel66(M, nax) atol = 1.0e-12
    end

    @testset "multi-block matrix : block independence + idempotence" begin
        n = 3
        M = zeros(6n, 6n)
        blocks = [_rand_minor_mandel(rng) for _ in 1:n, _ in 1:n]
        for i in 1:n, j in 1:n
            M[(6i - 5):(6i), (6j - 5):(6j)] = blocks[i, j]
        end
        out = _ti_project_blocks(M, nax)
        for i in 1:n, j in 1:n
            @test out[(6i - 5):(6i), (6j - 5):(6j)] ≈
                MCr.ti_average_mandel66(blocks[i, j], nax) atol = 1.0e-12
        end
        # idempotence
        @test _ti_project_blocks(out, nax) ≈ out atol = 1.0e-12
        # ISO ∘ TI == ISO (block-wise)
        @test _iso_project_blocks(out) ≈ _iso_project_blocks(M) atol = 1.0e-12
    end

    @testset "_maybe_symmetrize_alv dispatch incl. TISymmetrize" begin
        M = zeros(12, 12)
        for i in 1:2, j in 1:2
            M[(6i - 5):(6i), (6j - 5):(6j)] = _rand_minor_mandel(rng)
        end
        @test _maybe_symmetrize_alv(M, NoSymmetrize()) === M
        @test _maybe_symmetrize_alv(M, IsoSymmetrize()) ≈ _iso_project_blocks(M)
        @test _maybe_symmetrize_alv(M, TISymmetrize(ez)) ≈ _ti_project_blocks(M, ez)
        @test _maybe_symmetrize_alv(M, TISymmetrize(nax)) ≈ _ti_project_blocks(M, nax)
    end
end
