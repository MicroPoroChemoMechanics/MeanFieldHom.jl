# Adding a new inclusion

1. Choose the correct abstract supertype:
   [`AbstractEllipsoidalInclusion`](@ref),
   [`AbstractCrack`](@ref), or
   [`AbstractLayeredInclusion`](@ref).
2. Define your struct and implement the four interface methods:
   [`MeanFieldHom.Core.dimension`](@ref),
   [`MeanFieldHom.Core.element_type`](@ref),
   [`MeanFieldHom.Core.inclusion_basis`](@ref),
   [`MeanFieldHom.Core.shape_trait`](@ref).
3. Add `_kernel(::YourInclusion, C₀, ::Analytical; kw...)` methods for
   every matrix-symmetry class where an analytical formula exists, and
   `_kernel(::YourInclusion, C₀, ::Residue)` / `_kernel(..., ::DECUHR)`
   for numerical fallbacks.
4. If the new inclusion has a specialised TI-aligned path (cf. the
   `AbstractCrack` / `TensTI{4}` rule), register the extra dispatch
   rule at the end of your sub-module file via `import ..Core: _resolve_algo`.
5. Add at least one unit test under `test/<SubModule>/`.
