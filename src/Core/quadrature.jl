# =============================================================================
#  quadrature.jl
#
#  Backend seam for the optional DECUHR cubature path.
#
#  This file used to also export a `_quadgk` wrapper "all downstream
#  sub-modules should go through".  Nothing ever did — every call site uses
#  `QuadGK.quadgk` directly — so it was removed rather than left as a
#  misleading invitation.  Per-node counting is provided instead by
#  `Core._counted_quadgk` (see `counters.jl`), which the quadrature-heavy
#  kernels do call.
# =============================================================================

"""
    _decuhr_cubature(integrand, lb, ub; singul, alpha, wrksub, abstol, reltol, maxiters)

Backend seam for the optional **DECUHR** cubature path. The real
implementation (via `Integrals.solve(prob, DECUHR.DecuhrAlgorithm(...))`)
lives in the package extension `MeanFieldHomDECUHRExt`, which is loaded only
when both `DECUHR` and `Integrals` are available. Returns the raw solution
vector `sol.u`.

This fallback method is hit when the extension is **not** loaded and raises
an informative error. To use the `:decuhr` method, run `import DECUHR, Integrals`
first, or use the built-in `method = :nestedquadgk` alternative (QuadGK-based,
ForwardDiff-compatible, no extra dependency).
"""
_decuhr_cubature(args...; kwargs...) = error(
    "The `:decuhr` backend requires the DECUHR extension: run " *
        "`import DECUHR, Integrals` first, or use `method = :nestedquadgk` " *
        "(built-in, no extra dependency)."
)
