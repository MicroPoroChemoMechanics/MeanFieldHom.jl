# Architecture

`MeanFieldHom` is organised around a single principle:

> every high-level entry point dispatches via
> [`MeanFieldHom.Core._resolve_algo`](@ref) on
> `(Val(method), inclusion, C₀)`.

The resolved [`AbstractAlgorithm`](@ref) instance is then passed to an
internal `_kernel` method table maintained by each sub-module.
Sub-modules may *extend* (but not redefine) both `_resolve_algo` and
`_kernel`.

## Sub-module responsibilities

| Sub-module       | Exports                                                                             |
| ---------------- | ----------------------------------------------------------------------------------- |
| `Core`           | abstractions, traits, Newton potentials, Green kernel helpers, moduli extractors    |
| `Elasticity`     | `Ellipsoid`, auxiliary tensors, `hill_tensor` + 3D/2D kernels                        |
| `Cracks`         | `EllipticCrack`, `RibbonCrack`, `cod_tensor`, `compliance_contribution`, `sif`, `dif` |
| `Conductivity`   | additional `_kernel` methods for 2nd-order transport tensors                         |
| `Schemes`        | placeholder (no public API yet)                                                      |
