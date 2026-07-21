# =============================================================================
#  quadrature.jl
#
#  Thin wrapper around `QuadGK.quadgk` that normalises the keyword-argument
#  names (`abstol`, `reltol`, `maxiters`).  All downstream sub-modules should
#  go through this helper so that any future change in backend is localised
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
