# Architecture

`MeanFieldHom` is organized around a single principle:

> every high-level entry point dispatches via
> `MeanFieldHom.Core._resolve_algo` on
> `(Val(method), inclusion, C₀)`.

The resolved [`AbstractAlgorithm`](@ref) instance is then passed to an
internal `_kernel` method table maintained by each sub-module.
Sub-modules may *extend* (but not redefine) both `_resolve_algo` and
`_kernel`.

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
| `FiniteElements`   | inclusions solved by finite elements (`FEEllipticCrack`, `FEExcenteredSphere`); the solvers live in `MeanFieldHomFerriteExt` |

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
See [Adding a new inclusion](@ref), [Adding a new algorithm](@ref) and
[Adding a homogenization scheme](@ref).
