# =============================================================================
#  MeanFieldHomNonlinearSolveExt.jl
#
#  Weak extension activated when `NonlinearSolve.jl` is loaded together with
#  `MeanFieldHom`. Extends `MeanFieldHom.Schemes._solve_sc` to support every
#  algorithm of the SciML `NonlinearSolve.jl` package
#  (`NewtonRaphson()`, `TrustRegion()`, `Anderson()`, `LimitedMemoryNewtonRaphson()`,
#  …) for the SelfConsistent / AsymmetricSelfConsistent schemes.
#
#  Usage:
#      using MeanFieldHom
#      using NonlinearSolve
#      sc = SelfConsistent(; algorithm = NewtonRaphson(),
#                            abstol = 1e-12, maxiters = 200)
#      C_eff = homogenize(rve, sc)
# =============================================================================

module MeanFieldHomNonlinearSolveExt

using MeanFieldHom
using NonlinearSolve
using TensND

# `NewtonDefault` activates the SciML default Newton-Raphson when this
# extension is loaded.
function MeanFieldHom.Schemes._solve_sc(
        ::MeanFieldHom.Schemes.NewtonDefault,
        step,
        x0::TensND.AbstractTens;
        abstol::Real = 1.0e-10, maxiters::Int = 100, kw...,
    )
    return _solve_sc_nls(NewtonRaphson(), step, x0; abstol, maxiters, kw...)
end

# Generic SciML algorithm dispatch.
function MeanFieldHom.Schemes._solve_sc(
        algo::NonlinearSolve.AbstractNonlinearAlgorithm,
        step,
        x0::TensND.AbstractTens;
        abstol::Real = 1.0e-10, maxiters::Int = 100, kw...,
    )
    return _solve_sc_nls(algo, step, x0; abstol, maxiters, kw...)
end

# ── Internal: vec/unvec wrapper around the tensor `step` ───────────────────

function _solve_sc_nls(algo, step, x0::TensND.AbstractTens{order, dim, T};
                       abstol, maxiters, kw...) where {order, dim, T}
    basis_ref = TensND.get_basis(x0)
    v0 = vec(collect(get_array(x0)))
    function residual!(R, v, _)
        arr = reshape(v, ntuple(_ -> dim, order))
        Cv  = TensND.Tens(arr, basis_ref)
        R  .= vec(collect(get_array(step(Cv)))) .- v
        return R
    end
    prob = NonlinearProblem(residual!, v0, nothing)
    sol  = solve(prob, algo; abstol, maxiters, kw...)
    arr  = reshape(sol.u, ntuple(_ -> dim, order))
    return TensND.Tens(arr, basis_ref)
end

end # module
