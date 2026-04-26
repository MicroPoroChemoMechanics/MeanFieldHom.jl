# =============================================================================
#  sensitivities.jl — autodiff wrappers around `homogenize`.
#
#  Declares the public stubs `derivative`, `gradient`, `jacobian` and
#  `sensitivity`. Concrete methods are provided by the weak extension
#  `MeanFieldHomForwardDiffExt` activated by `using ForwardDiff`; without
#  ForwardDiff loaded, calls fall through to the no-method case and the
#  caller gets a clear error.
#
#  The indirection keeps the package lightweight for users who never need
#  autodiff and matches the pattern already used for NonlinearSolve and
#  SymPy.
# =============================================================================

"""
    derivative(rve, scheme, param::AbstractParameter;
               output = :C, indexer = identity, kw...) -> Number | AbstractTens

Scalar derivative of `homogenize(rve, scheme; property=output)` with respect
to the parameter selected by the lens `param`.

`indexer` is a function that extracts a *scalar* from the effective tensor
(for example `C -> get_data(C)[1] / 3` for the bulk modulus). Without an
`indexer`, the derivative of the full tensor is returned — this is valid as
long as ForwardDiff can propagate through the result (typically when the
output is itself a structured object whose components ForwardDiff handles
natively, via the `_data` Dual fields).

Extra `kw...` (e.g. `Chunk`, `Tag`) are forwarded to `ForwardDiff`.

The method is only available with `using ForwardDiff`; calling it without
the extension loaded raises an explicit error pointing at the extension.

See also [`gradient`](@ref), [`jacobian`](@ref), [`sensitivity`](@ref).
"""
function derivative end

"""
    gradient(rve, scheme, params::AbstractVector{<:AbstractParameter};
             output = :C, indexer, kw...) -> Vector

Gradient of a scalar functional (extracted via `indexer`) of the effective
tensor returned by `homogenize(rve, scheme; property=output)`, with respect
to the vector of lenses `params`.

Implementation: `ForwardDiff.gradient` with automatic chunking
(`ForwardDiff.Chunk(length(params))` via the standard `pickchunksize`
rule).

Without `using ForwardDiff`: error.
"""
function gradient end

"""
    jacobian(rve, scheme, params::AbstractVector{<:AbstractParameter};
             output = :C, indexer = identity, kw...) -> Array

Jacobian (returned as a flat array `length(output_flat) × length(params)`)
of `homogenize(rve, scheme; property=output)` with respect to the vector of
lenses `params`.

`indexer` can be used to reduce the output to a sub-tensor before
flattening. Without an `indexer`, the full effective tensor is flattened
via `get_array` then `vec`.

Without `using ForwardDiff`: error.
"""
function jacobian end

"""
    sensitivity(f, x₀; kind = :auto, kw...)

Generic ForwardDiff wrapper. `f(x)` is a user-supplied closure that builds
or perturbs an RVE, runs `homogenize`, and returns a scalar (or a tensor).

`kind`:
- `:derivative` — scalar `x₀`, scalar output; returns `f'(x₀)`.
- `:gradient`   — vector `x₀`, scalar output; returns `∇f(x₀)`.
- `:jacobian`   — vector `x₀`, tensor output; returns the flattened Jacobian.
- `:auto` (default) — pick one of the above from the types of `x₀` and `f(x₀)`.

This entry point covers cases the `AbstractParameter` lenses cannot
express (compound parametrisations, user-defined inclusion types without
an exposed scalar field, etc.).

Without `using ForwardDiff`: error.
"""
function sensitivity end

# Note: concrete methods are defined in `ext/MeanFieldHomForwardDiffExt.jl`
# and become available only after `using ForwardDiff`. Calling these stubs
# without the extension loaded raises a `MethodError` — load `ForwardDiff`
# to resolve.
