# Diagnostic: dump Julia's trapezoidal_matrix and α/β for a step kernel.
import Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."); io = devnull)

using MeanFieldHom
using TensND
using LinearAlgebra
using Printf

T_grid = [0.5, 0.7, 1.0, 1.5, 2.0]
n = length(T_grid)

k1 = 5.0 / (3 * (1 - 2*0.3))
mu1 = 5.0 / (2 * (1 + 0.3))
eta1 = 1.0; gamma1 = 1.67
kp = 1e-8 / (3 * (1 - 2*0.2))
mup = 1e-8 / (2 * (1 + 0.2))
t_set = 0.669

C_p_tens = TensISO{3}(3 * kp, 2 * mup)
R1 = maxwell_iso(k1, mu1, eta1, gamma1)

# Build the same step kernel as Python (history-dep).
step_law = ViscoLaw(function (t, tp)
    if t < tp
        return zero(C_p_tens)
    elseif tp ≥ t_set
        return R1.eval_fun(t, tp)
    else
        return C_p_tens
    end
end, :relaxation)

R = trapezoidal_matrix(step_law, T_grid)
α, β = iso_params_from_blocks(R)

println("Julia trapezoidal α (3K) for step kernel:")
display(α)
println("\nJulia trapezoidal β (2μ) for step kernel:")
display(β)

# Print first 6×6 block too.
println("\nFirst 6×6 block (M[1:6, 1:6]):")
display(R[1:6, 1:6])
println("\nM[7:12, 1:6] (block i=2, j=1):")
display(R[7:12, 1:6])
println("\nM[7:12, 7:12] (block i=2, j=2):")
display(R[7:12, 7:12])
