import Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."); io = devnull)
using MeanFieldHom, LinearAlgebra, Printf

# Maxwell matrix M_κ and M_μ with different decay times.
T_grid = [0.0, 0.5, 1.0, 1.5, 2.0]   # original bench_layered_alv times
n = length(T_grid)
k0 = 1.0; mu0 = 0.5; eta_k = 0.6; eta_mu = 2.0  # original bench_layered_alv physics
make_R0() = maxwell_iso(k0, mu0, eta_k, eta_mu)
R0 = trapezoidal_matrix(make_R0(), T_grid)
α0, β0 = iso_params_from_blocks(R0)
M_κ_0 = α0 ./ 3
M_μ_0 = β0 ./ 2

# Test commutativity.
println("M_κ_0 * M_μ_0 - M_μ_0 * M_κ_0:")
display(M_κ_0 * M_μ_0 - M_μ_0 * M_κ_0)

# Compute num = 9·I + 4·M_μ_0 and S_b = 3·M_κ_0 + 4·M_μ_0.
num = 9 .* Matrix(I, n, n) .+ 4 .* M_μ_0
S_b = 3 .* M_κ_0 .+ 4 .* M_μ_0

# T_right = num · S_b^{-1} (volterra_divide).
T_right = MeanFieldHom.Viscoelasticity.volterra_divide(num, S_b; block_size = 1)
# T_left = S_b^{-1} · num.
T_left = volterra_inverse(S_b; block_size = 1) * num
# Check.
println("\nT_right (= num · S_b^{-1}):")
display(T_right)
println("\nT_left (= S_b^{-1} · num):")
display(T_left)
@printf "max |T_right - T_left| = %.6e\n" maximum(abs, T_right - T_left)

# Sanity: T_right * S_b should equal num.
@printf "‖T_right * S_b - num‖ = %.6e\n" maximum(abs, T_right * S_b - num)
@printf "‖S_b * T_left - num‖ = %.6e\n" maximum(abs, S_b * T_left - num)
