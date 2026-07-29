# =============================================================================
#  common.jl — what every finite-element inclusion shares.
#
#  A finite-element evaluation costs one assembly, one factorization and a
#  handful of solves, while iterative schemes (self-consistent, differential)
#  ask for the same tensors again at every iteration with a *slightly*
#  different reference medium.  Memoizing on that reference is what keeps those
#  schemes usable; one-shot schemes (dilute, Mori-Tanaka, Maxwell, PCW) only
#  ever hit one key.
# =============================================================================

"""
    FECache()

Mutable side-store of a finite-element inclusion: the assembled discretization
(built once, on first use) and the response tensors already computed, keyed on
the reference medium.

`assemblies` counts the factorizations actually performed — used by the tests
to prove the memoization works. The field is deliberately kept out of the
struct's numeric fields so that it never interferes with `ForwardDiff`
reconstruction of a geometry parameter.
"""
mutable struct FECache
    setup::Any
    tensors::Dict{Any, Any}
    assemblies::Int
end

FECache() = FECache(nothing, Dict{Any, Any}(), 0)

"""
    fe_assembly_count(incl) -> Int

Number of finite-element assemblies (equivalently, factorizations) actually
performed for `incl` so far. Every distinct reference medium costs one; a
repeat costs none. Useful to check that the memoization of [`FECache`](@ref)
is doing its job.
"""
fe_assembly_count(incl) = _fe_cache(incl).assemblies

"""
    fe_reset!(incl) -> incl

Drop the cached discretization and every memoized response tensor.
"""
function fe_reset!(incl)
    c = _fe_cache(incl)
    c.setup = nothing
    empty!(c.tensors)
    c.assemblies = 0
    return incl
end

"Cache of a finite-element inclusion; every such type defines one method."
_fe_cache(incl) = throw(
    ArgumentError("$(typeof(incl)) is not a finite-element inclusion (no `FECache`)")
)
