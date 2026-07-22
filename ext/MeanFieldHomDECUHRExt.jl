# =============================================================================
#  MeanFieldHomDECUHRExt.jl
#
#  Weak extension activated when both `DECUHR` and `Integrals` are loaded.
#  It implements the backend seam `MeanFieldHom.Core._decuhr_cubature`, which
#  is a fallback error stub in the core package (see `src/Core/quadrature.jl`).
#
#  With this extension loaded, the `method = :decuhr` path of `hill_tensor`
#  and `cod_tensor` becomes available. Without it, the core stub raises an
#  informative error and users rely on the built-in `method = :nestedquadgk`
#  alternative (QuadGK-based, ForwardDiff-compatible, no extra dependency).
# =============================================================================

module MeanFieldHomDECUHRExt

using MeanFieldHom
using DECUHR: DecuhrAlgorithm
import Integrals

# Real implementation of the DECUHR cubature seam. More specific than the
# core `_decuhr_cubature(args...; kwargs...)` fallback, so it is selected
# whenever this extension is loaded.
function MeanFieldHom.Core._decuhr_cubature(
        integrand, lb::AbstractVector, ub::AbstractVector;
        singul::Int = 1,
        alpha::Float64 = 0.0,
        wrksub::Int = 50_000,
        abstol::Real = 1.0e-8,
        reltol::Real = 1.0e-6,
        maxiters::Int = 100_000,
    )
    prob = Integrals.IntegralProblem(integrand, (lb, ub))
    sol = Integrals.solve(
        prob,
        DecuhrAlgorithm(singul = singul, alpha = alpha, wrksub = wrksub);
        abstol = abstol, reltol = reltol, maxiters = maxiters,
    )
    (
        sol.retcode == Integrals.ReturnCode.Success ||
            sol.retcode == Integrals.ReturnCode.MaxIters
    ) || error("DECUHR failed: retcode = $(sol.retcode)")
    return sol.u
end

end # module MeanFieldHomDECUHRExt
