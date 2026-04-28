# API — Sensitivities

Public lenses and autodiff entry points provided by
`MeanFieldHom.Schemes`.  `ForwardDiff` is a **strong dependency** of
`MeanFieldHom` since v0.7.0 — the four `derivative` / `gradient` /
`jacobian` / `sensitivity` functions are available out of the box,
and the built-in [`SelfConsistent`](@ref) Newton-Raphson solver
([`NewtonDefault`](@ref)) uses the same machinery internally.

## Lenses

```@docs
AbstractParameter
AmountParameter
PropertyParameter
GeometryParameter
DistributionShapeParameter
amount
property
geometry
shape_param
get_param
set_param
```

## Autodiff entry points

```@docs
derivative
gradient
jacobian
sensitivity
```
