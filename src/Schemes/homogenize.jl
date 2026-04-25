# =============================================================================
#  homogenize.jl — central public entry point for homogenisation schemes.
#
#  Each concrete scheme implements a method
#       _evaluate(rve::RVE, ::ConcreteScheme, ::Val{property}; kw...)
#  in its own file (`voigt.jl`, `dilute.jl`, …). The `Val{property}` form
#  enables compile-time specialisation on the property name (`:C`, `:K`)
#  without paying the Symbol-dispatch cost at runtime.
# =============================================================================

"""
    homogenize(rve::RVE, scheme::HomogenizationScheme; property::Symbol = :C, kw...)
    homogenize(rve::RVE, scheme::Symbol; property::Symbol = :C, kw...)

Compute the effective `property` of `rve` under the chosen `scheme`.

`property` selects which stored phase property is homogenised:

- `:C` (default) — 4th-order elastic stiffness ;
- `:K` — 2nd-order conductivity / diffusivity ;
- arbitrary user-defined symbols are supported as long as every phase
  carries a tensor under that key.

The scheme can be passed either as a *type instance* (full control of
options — `MoriTanaka()`, `SelfConsistent(algorithm = NewtonDefault(), abstol = 1e-12)`,
…) or as a `Symbol` shortcut for the default constructor. The
canonical Symbol aliases are lowercase (`:mt`, `:sc`, `:voigt`, …) to
match the algorithm-method symbols (`:auto`, `:residues`, `:decuhr`,
…) used by the underlying Hill-tensor backends; CamelCase and
upper-case ECHOES-style codes are also accepted (see [`SCHEME_ALIAS`](@ref)).

Extra `kw...` are forwarded to the scheme's `_evaluate` method
(typically `abstol` / `reltol` / `maxiters` for iterative schemes,
`method = :auto | :decuhr | …` for the underlying Hill-tensor backend).
"""
function homogenize(rve::RVE, scheme::HomogenizationScheme;
                    property::Symbol = :C, kw...)
    validate_rve(rve)
    return _evaluate(rve, scheme, Val(property); kw...)
end

"""
    SCHEME_ALIAS :: Dict{Symbol, Type{<:HomogenizationScheme}}

Maps Symbol shortcuts to concrete scheme types for the convenience
overload `homogenize(rve, ::Symbol)`.

The canonical aliases are **all lowercase** (`:voigt`, `:reuss`,
`:dilute`, `:dilute_dual`, `:mori_tanaka`, `:mt`, `:maxwell`, `:pcw`,
`:sc`, `:asc`, `:differential`, `:diff`) for consistency with the
algorithm-method symbols accepted elsewhere in the package
(`:auto`, `:residues`, `:decuhr`, `:nestedquadgk`, `:analytical`).
The CamelCase forms (`:MoriTanaka`, `:Differential`) and the
ECHOES-compatible upper-case codes (`:MT`, `:DIFF`, …) are kept as
extra aliases for ease of porting.
"""
const SCHEME_ALIAS = Dict{Symbol, Type{<:HomogenizationScheme}}(
    # Voigt
    :voigt => Voigt, :Voigt => Voigt, :VOIGT => Voigt, :v => Voigt, :V => Voigt,
    # Reuss
    :reuss => Reuss, :Reuss => Reuss, :REUSS => Reuss, :r => Reuss, :R => Reuss,
    # Dilute
    :dilute => Dilute, :Dilute => Dilute, :dil => Dilute, :Dil => Dilute, :DIL => Dilute,
    # Dilute Dual
    :dilute_dual => DiluteDual, :DiluteDual => DiluteDual,
    :dild => DiluteDual, :DilD => DiluteDual, :DILD => DiluteDual,
    # Mori-Tanaka
    :mori_tanaka => MoriTanaka, :moritanaka => MoriTanaka, :MoriTanaka => MoriTanaka,
    :mt => MoriTanaka, :MT => MoriTanaka,
    # Maxwell
    :maxwell => Maxwell, :Maxwell => Maxwell, :max => Maxwell, :Max => Maxwell, :MAX => Maxwell,
    # Ponte-Castañeda & Willis
    :pcw => PonteCastanedaWillis, :PCW => PonteCastanedaWillis,
    :PonteCastanedaWillis => PonteCastanedaWillis,
    :ponte_castaneda_willis => PonteCastanedaWillis,
    # Self-Consistent
    :self_consistent => SelfConsistent, :SelfConsistent => SelfConsistent,
    :sc => SelfConsistent, :SC => SelfConsistent,
    # Asymmetric Self-Consistent
    :asymmetric_self_consistent => AsymmetricSelfConsistent,
    :AsymmetricSelfConsistent => AsymmetricSelfConsistent,
    :asc => AsymmetricSelfConsistent, :ASC => AsymmetricSelfConsistent,
    # Differential
    :differential => DifferentialScheme, :Differential => DifferentialScheme,
    :diff => DifferentialScheme, :Diff => DifferentialScheme, :DIFF => DifferentialScheme,
)

function homogenize(rve::RVE, scheme::Symbol; kw...)
    haskey(SCHEME_ALIAS, scheme) ||
        throw(ArgumentError("unknown scheme :$(scheme); see MeanFieldHom.Schemes.SCHEME_ALIAS"))
    return homogenize(rve, SCHEME_ALIAS[scheme](); kw...)
end

# =============================================================================
#  Default fallback — explicit "not yet implemented" error so that an
#  unimplemented scheme does not silently dispatch elsewhere.
# =============================================================================

"""
    _evaluate(rve, scheme, ::Val{p}; kw...) -> AbstractTens

Internal entry point that each concrete scheme implements. This generic
fallback throws an explicit `ErrorException` so that a missing
specialisation is reported clearly instead of dispatching to the wrong
method.
"""
function _evaluate(rve::RVE, scheme::HomogenizationScheme, ::Val{p}; kw...) where {p}
    error("homogenize: scheme $(typeof(scheme)) does not yet implement property :$(p)")
end
