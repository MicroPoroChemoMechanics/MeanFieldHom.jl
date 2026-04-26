# =============================================================================
#  MeanFieldHomForwardDiffExt.jl
#
#  Weak extension activated when `ForwardDiff.jl` is loaded together with
#  `MeanFieldHom`. Provides convenient wrappers around `ForwardDiff` for
#  computing sensitivities of `homogenize(rve, scheme)` with respect to any
#  scalar input parameter, designated by an `AbstractParameter` lens or by a
#  user-supplied closure.
#
#  Public entry points exposed to the main package (defined as stubs in
#  `src/Schemes/sensitivities.jl`) :
#      derivative(rve, scheme, param;       output=:C, indexer=identity, kw...)
#      gradient(rve,   scheme, params;      output=:C, indexer=identity, kw...)
#      jacobian(rve,   scheme, params;      output=:C, indexer=identity, kw...)
#      sensitivity(f, x₀; kind=:auto, kw...)
#
#  Usage:
#      using MeanFieldHom
#      using ForwardDiff
#      ∂kf = derivative(rve, MoriTanaka(), amount(:I);
#                       indexer = C -> get_data(C)[1] / 3)
# =============================================================================

module MeanFieldHomForwardDiffExt

using MeanFieldHom
using ForwardDiff
using TensND

import MeanFieldHom: derivative, gradient, jacobian, sensitivity
import MeanFieldHom: AbstractParameter, RVE, HomogenizationScheme,
    get_param, set_param, homogenize
import MeanFieldHom.Schemes: _set_many

# =============================================================================
#  Helpers
# =============================================================================

# Extract a flat-vector view of a homogenization output: handle AbstractTens,
# AbstractArray, and scalar uniformly.
_flatten(x::Number)        = [x]
_flatten(x::AbstractArray) = vec(x)
_flatten(x::TensND.AbstractTens) = vec(collect(TensND.get_array(x)))

# =============================================================================
#  derivative — scalar input, scalar output via indexer
# =============================================================================

function MeanFieldHom.derivative(rve::RVE, scheme::HomogenizationScheme,
                                 p::AbstractParameter;
                                 output::Symbol = :C, indexer = identity, kw...)
    x₀ = get_param(rve, p)
    f = x -> begin
        rve′ = set_param(rve, p, x)
        out  = homogenize(rve′, scheme; property = output, kw...)
        return indexer(out)
    end
    return ForwardDiff.derivative(f, x₀)
end

# =============================================================================
#  gradient — vector input, scalar output via indexer
# =============================================================================

function MeanFieldHom.gradient(rve::RVE, scheme::HomogenizationScheme,
                               ps::AbstractVector{<:AbstractParameter};
                               output::Symbol = :C, indexer = identity,
                               chunk = nothing, kw...)
    x₀ = [get_param(rve, p) for p in ps]
    f = xs -> begin
        rve′ = _set_many(rve, ps, xs)
        out  = homogenize(rve′, scheme; property = output, kw...)
        return indexer(out)
    end
    if chunk === nothing
        return ForwardDiff.gradient(f, x₀)
    else
        cfg = ForwardDiff.GradientConfig(f, x₀, chunk)
        return ForwardDiff.gradient(f, x₀, cfg)
    end
end

# =============================================================================
#  jacobian — vector input, vector / tensor output
# =============================================================================

function MeanFieldHom.jacobian(rve::RVE, scheme::HomogenizationScheme,
                               ps::AbstractVector{<:AbstractParameter};
                               output::Symbol = :C, indexer = identity,
                               chunk = nothing, kw...)
    x₀ = [get_param(rve, p) for p in ps]
    f = xs -> begin
        rve′ = _set_many(rve, ps, xs)
        out  = homogenize(rve′, scheme; property = output, kw...)
        return _flatten(indexer(out))
    end
    if chunk === nothing
        return ForwardDiff.jacobian(f, x₀)
    else
        cfg = ForwardDiff.JacobianConfig(f, x₀, chunk)
        return ForwardDiff.jacobian(f, x₀, cfg)
    end
end

# Single-parameter convenience: pass one AbstractParameter, return jacobian
# along that single dimension (still a Matrix of size (length(output_flat), 1)).
MeanFieldHom.jacobian(rve::RVE, scheme::HomogenizationScheme,
                      p::AbstractParameter; kw...) =
    MeanFieldHom.jacobian(rve, scheme, [p]; kw...)

# =============================================================================
#  sensitivity — generic ForwardDiff wrapper around a user-supplied closure
# =============================================================================

function MeanFieldHom.sensitivity(f, x₀; kind::Symbol = :auto, kw...)
    actual = kind === :auto ? _autoselect(f, x₀) : kind
    if actual === :derivative
        return ForwardDiff.derivative(f, x₀)
    elseif actual === :gradient
        return ForwardDiff.gradient(f, x₀)
    elseif actual === :jacobian
        # Wrap f so the output is always a flat vector
        wrapper = x -> _flatten(f(x))
        return ForwardDiff.jacobian(wrapper, x₀)
    else
        throw(ArgumentError("sensitivity: unknown kind=:$(actual); expected :derivative, :gradient or :jacobian"))
    end
end

function _autoselect(f, x₀::Number)
    y = f(x₀)
    return y isa Number ? :derivative : :jacobian
end

function _autoselect(f, x₀::AbstractVector)
    y = f(x₀)
    return y isa Number ? :gradient : :jacobian
end

_autoselect(f, x₀) = throw(ArgumentError(
    "sensitivity: cannot autoselect kind for x₀ of type $(typeof(x₀)); pass `kind=:derivative|:gradient|:jacobian`"
))

end # module
