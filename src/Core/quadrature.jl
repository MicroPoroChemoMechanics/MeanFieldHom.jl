# =============================================================================
#  quadrature.jl
#
#  Thin wrappers around `QuadGK.quadgk` and `DECUHR.hcubature` that
#  normalise the keyword-argument names (`abstol`, `reltol`, `maxiters`).
#  All downstream sub-modules should go through these helpers so that any
#  future change in backend (e.g. swapping DECUHR for `HCubature.jl`) is
#  localised to this file.
# =============================================================================

"""
    _quadgk(f, a, b; abstol, reltol, maxiters)

Thin wrapper around `QuadGK.quadgk` that exposes the same keyword names
as the rest of `MeanFieldHom` (`abstol`, `reltol`, `maxiters`).  Returns
the `(value, error)` tuple produced by `quadgk`.
"""
@inline function _quadgk(f, a, b;
                         abstol::Real  = 1e-8,
                         reltol::Real  = 1e-6,
                         maxiters::Int = 1_000_000)
    return QuadGK.quadgk(f, a, b; atol=abstol, rtol=reltol, maxevals=maxiters)
end

"""
    _hcubature(f, lo, hi; abstol, reltol, maxiters)

Thin wrapper around `DECUHR.hcubature` — same keyword normalisation as
[`_quadgk`](@ref).  Not used by the kernels shipped with the first
release (they rely on nested `QuadGK` for better ForwardDiff support),
but kept here so that future inclusions can opt in to a cubature
backend without touching the sub-module code.
"""
@inline function _hcubature(f, lo, hi;
                            abstol::Real  = 1e-8,
                            reltol::Real  = 1e-6,
                            maxiters::Int = 1_000_000)
    return DECUHR.hcubature(f, lo, hi;
                            atol=abstol, rtol=reltol, maxevals=maxiters)
end
