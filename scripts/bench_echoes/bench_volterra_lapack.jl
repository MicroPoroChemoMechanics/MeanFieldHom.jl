# Compare hand-rolled forward-substitution Volterra ops vs LAPACK
# (LowerTriangular(...) \ B and inv(LowerTriangular(...))).
#
# Goal: decide whether to refactor `volterra_inverse`,
# `volterra_left_divide`, `volterra_divide` to dispatch to BLAS/LAPACK on
# `LowerTriangular` wrappers.
import Pkg
Pkg.activate(@__DIR__; io = devnull)

using LinearAlgebra
using Printf
using BenchmarkTools
using MeanFieldHom

# Build a well-conditioned lower-triangular Volterra-like matrix of size n.
function build_lower(n::Int; T::Type = Float64, scale = T(1.0))
    M = zeros(T, n, n)
    @inbounds for j in 1:n, i in j:n
        # diagonal-dominant lower triangular with decay
        if i == j
            M[i, j] = scale * (1 + 0.1 * i)
        else
            M[i, j] = scale * 0.05 * exp(-T(i - j) / 4)
        end
    end
    return M
end

println("=== volterra_inverse vs inv(LowerTriangular(.)) ===")
@printf "%6s | %14s | %14s | %14s | %8s\n" "n" "hand-rolled" "LAPACK" "diff (rel)" "speedup"
for n in (10, 50, 100, 200, 500, 1000)
    M = build_lower(n)
    invM_hand = volterra_inverse(M; block_size = 1)
    invM_lap = inv(LowerTriangular(M))
    rel = norm(invM_hand .- invM_lap) / norm(invM_hand)

    t_hand = @belapsed volterra_inverse($M; block_size = 1) samples = 5 evals = 1
    t_lap = @belapsed inv(LowerTriangular($M)) samples = 5 evals = 1
    @printf "%6d | %12.3f μs | %12.3f μs | %14.2e | ×%6.2f\n" n (t_hand * 1.0e6) (t_lap * 1.0e6) rel (t_hand / t_lap)
end

println()
println("=== volterra_left_divide(S, M) vs LowerTriangular(S) \\ M ===")
@printf "%6s | %14s | %14s | %14s | %8s\n" "n" "hand-rolled" "LAPACK" "diff (rel)" "speedup"
for n in (10, 50, 100, 200, 500, 1000)
    S = build_lower(n; scale = 1.0)
    M = build_lower(n; scale = 0.7)
    T_hand = volterra_left_divide(S, M; block_size = 1)
    T_lap = LowerTriangular(S) \ M
    rel = norm(T_hand .- T_lap) / norm(T_hand)

    t_hand = @belapsed volterra_left_divide($S, $M; block_size = 1) samples = 5 evals = 1
    t_lap = @belapsed LowerTriangular($S) \ $M samples = 5 evals = 1
    @printf "%6d | %12.3f μs | %12.3f μs | %14.2e | ×%6.2f\n" n (t_hand * 1.0e6) (t_lap * 1.0e6) rel (t_hand / t_lap)
end

println()
println("=== volterra_divide(M, S) vs M / LowerTriangular(S) ===")
@printf "%6s | %14s | %14s | %14s | %8s\n" "n" "hand-rolled" "LAPACK" "diff (rel)" "speedup"
for n in (10, 50, 100, 200, 500, 1000)
    S = build_lower(n; scale = 1.0)
    M = build_lower(n; scale = 0.7)
    T_hand = volterra_divide(M, S; block_size = 1)
    T_lap = M / LowerTriangular(S)
    rel = norm(T_hand .- T_lap) / norm(T_hand)

    t_hand = @belapsed volterra_divide($M, $S; block_size = 1) samples = 5 evals = 1
    t_lap = @belapsed $M / LowerTriangular($S) samples = 5 evals = 1
    @printf "%6d | %12.3f μs | %12.3f μs | %14.2e | ×%6.2f\n" n (t_hand * 1.0e6) (t_lap * 1.0e6) rel (t_hand / t_lap)
end
