# =============================================================================
#  quadrature.jl
#
#  Thin wrapper around `QuadGK.quadgk` that normalizes the keyword-argument
#  names (`abstol`, `reltol`, `maxiters`).  All downstream sub-modules should
#  go through this helper so that any future change in backend is localized
#  to this file.
# =============================================================================

"""
    _quadgk(f, a, b; abstol, reltol, maxiters)

Thin wrapper around `QuadGK.quadgk` that exposes the same keyword names
as the rest of `MeanFieldHom` (`abstol`, `reltol`, `maxiters`).  Returns
the `(value, error)` tuple produced by `quadgk`.
"""
@inline function _quadgk(
        f, a, b;
        abstol::Real = 1.0e-8,
        reltol::Real = 1.0e-6,
        maxiters::Int = 1_000_000
    )
    return QuadGK.quadgk(f, a, b; atol = abstol, rtol = reltol, maxevals = maxiters)
end

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
