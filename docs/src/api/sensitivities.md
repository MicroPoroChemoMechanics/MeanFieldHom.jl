# API — Sensitivities

Public lenses and autodiff entry points provided by
`MeanFieldHom.Schemes` (the four `derivative` / `gradient` / `jacobian` /
`sensitivity` functions become callable only after `using ForwardDiff`,
via the `MeanFieldHomForwardDiffExt` weak extension).

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

## Autodiff entry points (require `using ForwardDiff`)

```@docs
derivative
gradient
jacobian
sensitivity
```
