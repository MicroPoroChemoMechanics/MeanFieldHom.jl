# Adding a new algorithm

1. Subtype [`MeanFieldHom.Core.AbstractAlgorithm`](@ref) with a plain
   singleton struct.
2. Register a dispatch rule in `Core.dispatch` or in the sub-module
   where the algorithm applies.
3. Add the corresponding `_kernel(inclusion, C₀, ::YourAlgorithm; kw...)`
   method implementations.
4. Document the AD / ForwardDiff compatibility in the algorithm's
   docstring.
