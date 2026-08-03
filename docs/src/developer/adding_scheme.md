# [Adding a homogenization scheme](@id dev-adding-scheme)

!!! note "Which cell does your scheme serve?"
    `homogenize` accepts any
    [`AbstractHomogenizationCell`](@ref MeanFieldHom.Core.AbstractHomogenizationCell)
    — an [`RVE`](@ref) or a `Laminate`. Schemes stay typed on the **concrete**
    cell they serve (`_evaluate(rve::RVE, ::MyScheme, …)`), so a scheme applied
    to a cell it does not serve falls through to the error-raising fallback
    instead of dispatching somewhere wrong. Type the *scheme* argument too:
    leaving it untyped is ambiguous with that fallback.

`MeanFieldHom.Schemes` ships Voigt/Reuss, dilute (direct and dual),
Mori-Tanaka, Maxwell, Ponte-Castañeda–Willis, self-consistent (symmetric and
asymmetric) and differential, each in an elastic and an ageing-viscoelastic
flavor. A new scheme slots in beside them:

1. Create `src/Schemes/<scheme_name>.jl` and `include` it from
   `src/Schemes/Schemes.jl`.
2. Declare the scheme type (`struct MyScheme <: HomogenizationScheme`) in
   `src/Schemes/scheme_types.jl`, and register its symbol alias in
   `SCHEME_ALIAS` (`src/Schemes/homogenize.jl`) if you want
   `homogenize(rve, :myscheme, :C)` to work.
3. Implement `_evaluate(rve, ::MyScheme, ::Val{p}; kw...)`. Provide **one
   method per tensor order** — dispatch on the order of `matrix_property(rve, p)`
   — so the scheme serves elasticity and transport from one implementation, as
   `_mt_dispatch` does in `mori_tanaka.jl`.
4. **Go through `contribution_helpers.jl`.** Never call
   `stiffness_contribution` and friends directly on a phase geometry: the
   `_phase_*` helpers apply the amount, honor the `symmetrize` setting, hand
   the kernel the correctly pre-projected reference medium, and branch on
   [`is_homogeneous_inclusion`](@ref). Two invariants stated in that file's
   header must hold in your kernel too — all helpers of a given evaluation must
   share the same projected `P₀`, and orientation averaging must never be
   folded into a tensor product.
5. Prefer the bundled seams (`_phase_dilute_and_contribution`,
   `_phase_dilute_and_stress_average`, `_phase_compliance_and_contribution`)
   when you need two objects per phase: they share the single expensive solve.
6. Consider what your kernel *requires* of an inclusion. Needing only the
   contribution tensors keeps the scheme usable with inclusions entered through
   gate C of the [inclusion contract](@ref dev-adding-inclusion); needing
   the concentration tensor restricts it.
7. Export through `src/Schemes/Schemes.jl` and re-export from
   `src/MeanFieldHom.jl`.
8. Add a unit test under `test/Schemes/` and a section in
   `docs/src/manual/schemes.md`.
