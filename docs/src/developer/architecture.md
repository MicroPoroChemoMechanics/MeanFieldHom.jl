# [Architecture](@id dev-architecture)

`MeanFieldHom` is organized around a single principle:

> every high-level entry point dispatches via
> `MeanFieldHom.Core._resolve_algo` on
> `(Val(method), inclusion, C₀)`.

The resolved [`AbstractAlgorithm`](@ref) instance is then passed to an
internal `_kernel` method table maintained by each sub-module.
Sub-modules may *extend* (but not redefine) both `_resolve_algo` and
`_kernel`.

## Cells: what `homogenize` accepts

`RVE` is no longer the only container. Everything `homogenize` accepts is an
[`AbstractHomogenizationCell`](@ref MeanFieldHom.Core.AbstractHomogenizationCell),
declared in `src/Core/cells.jl` alongside the rest of the package's generics:

| Cell | Morphology | Solved by |
| :--- | :--- | :--- |
| `Schemes.RVE` | random, through the Eshelby auxiliary problem | every mean-field scheme |
| `Laminates.Laminate` | periodic stack of parallel layers | `Laminated` (exact), `Voigt`, `Reuss` |

Only the **entry points** are typed on the supertype (`homogenize`,
`get_param`/`set_param`, `derivative`/`gradient`/`jacobian`, and the
`_evaluate` fallback). Every scheme kernel and every `_phase_*` helper stays
typed on `RVE`. That split is deliberate: it is what makes the abstraction
regression-free, and what makes a mis-applied scheme report itself.

`src/Core/cells.jl` also carries the **declarative multiscale** seam. A
property value may be a `Homogenized(cell, scheme)`, resolved lazily by
`resolve_property` — called from exactly two places, `phase_property` and
`layer_property` — and memoized per `(cell, key)` for the duration of one
`homogenize` call by a task-local `ScopedValue`. Two invariants to preserve
when touching it:

- **type inspections must go through the raw accessors**
  (`phase_property_raw`, `cell_container_property`); resolving would run a full
  inner homogenization just to look at a type;
- the cache scope must span the whole `_evaluate`, iterative solvers included,
  or a self-consistent loop would re-solve every nested cell once per
  iteration.


## Sub-module responsibilities

| Sub-module         | Exports                                                                                            |
| ------------------ | -------------------------------------------------------------------------------------------------- |
| `Elliptic`         | type-generic elliptic integrals (`ell_K`, `ell_E`, `ell_F`, `ell_RF`, `ell_RD`)                     |
| `Core`             | abstractions, traits, `_resolve_algo`, Newton potentials, Green kernel helpers, Kelvin dipole field, exact ISO/TI rotation averages, moduli extractors |
| `Elasticity`       | `Ellipsoid`, `Cylinder`, auxiliary tensors, `hill_tensor` + 3D/2D kernels                           |
| `Cracks`           | `EllipticCrack`, `RibbonCrack`, `cod_tensor`, `sif`, `dif`, and the crack methods of the contribution generics |
| `Conductivity`     | additional `_kernel` methods for 2nd-order transport tensors                                        |
| `LayeredSpheres`   | `LayeredSphere`, Hervé-Zaoui recurrences, five interface types, localization fields                 |
| `LayeredSpheroids` | `LayeredSpheroid` (confocal, conduction)                                                            |
| `Schemes`          | `RVE`/`Phase`, `homogenize`, every scheme type, exact symmetrization, ForwardDiff sensitivities     |
| `Viscoelasticity`  | ageing linear viscoelasticity (Volterra pipeline, ALV variant of every scheme)                      |

| `CustomInclusions` | the user-defined inclusion contract: `CustomInclusion`, `check_inclusion_interface` |
| `FiniteElements`   | inclusions solved by finite elements (`FEEllipticCrack`, `FEExcenteredSphere`); the physics lives here, the discretization in a backend extension (`MeanFieldHomFerriteExt`, `MeanFieldHomGridapExt`) |

Two files sit at the **top level** rather than in a sub-module, and are loaded
after every geometry sub-module on purpose: `localization.jl` and
`contribution.jl` implement generics declared in `Core` whose methods need
every sub-module's `_kernel` table to be visible. `CustomInclusions` and
`FiniteElements` are included after *those*, because their fallbacks `invoke`
the generic methods defined there.

## Extension points

Every generic an inclusion may implement is declared — as a bodyless `function`
— in `Core/abstractions.jl`, so that sub-modules **and user code outside the
package** attach their methods to one canonical function. Together with the
open `_resolve_algo` / `_kernel` tables and the neutral
[`AbstractCustomInclusion`](@ref) branch, that is the whole extension surface.
See [Adding a new inclusion](@ref dev-adding-inclusion), [Adding a new algorithm](@ref dev-adding-algorithm) and
[Adding a homogenization scheme](@ref dev-adding-scheme).
