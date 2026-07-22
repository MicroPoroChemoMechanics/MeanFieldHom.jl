# =============================================================================
#  MeanFieldHomNonlinearSolveExt.jl
#
#  Weak extension activated when `NonlinearSolve.jl` is loaded together
#  with `MeanFieldHom`.
#
#  STATUS — intentional placeholder (no-op).  The supported
#  self-consistent solvers are the two package built-ins,
#  `MeanFieldHom.Schemes._solve_sc(::AndersonDefault, …)` (Picard +
#  relaxation, Dual-safe) and `_solve_sc(::NewtonDefault, …)`
#  (ForwardDiff-Jacobian Newton with line search) — neither needs this
#  extension.  It is kept as a load hook for a future SciML
#  `NonlinearSolve` algorithm dispatch (trust-region, Anderson with
#  memory > 1); until that lands, loading it changes nothing and
#  `SelfConsistent(algorithm = <SciML alg>)` has no effect beyond the
#  built-ins.  See `docs/src/developer/roadmap.md` (§Open).
# =============================================================================

module MeanFieldHomNonlinearSolveExt

using MeanFieldHom
using NonlinearSolve
using TensND

end # module
