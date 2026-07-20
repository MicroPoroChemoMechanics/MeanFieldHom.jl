# =============================================================================
#  MeanFieldHomNonlinearSolveExt.jl
#
#  Weak extension activated when `NonlinearSolve.jl` is loaded together
#  with `MeanFieldHom`.  The built-in Newton-Raphson SC solver
#  (`MeanFieldHom.Schemes._solve_sc(::NewtonDefault, …)`) ships
#  with the package — no extension needed for the default Newton path.
#  Loading this extension is a no-op at the moment ; it is kept as a
#  hook for future SciML `NonlinearSolve` algorithm dispatch (e.g.
#  trust-region, Anderson with memory > 1).
# =============================================================================

module MeanFieldHomNonlinearSolveExt

using MeanFieldHom
using NonlinearSolve
using TensND

end # module
