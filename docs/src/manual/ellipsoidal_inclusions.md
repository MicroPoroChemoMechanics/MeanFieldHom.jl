# Ellipsoidal inclusions

```julia
using MeanFieldHom, TensND
E, ν = 210e3, 0.3
λ = E*ν/((1+ν)*(1-2ν)); μ = E/(2*(1+ν))
C₀ = TensISO{3}(3*(λ+2μ/3), 2μ)

# Sphere
hill_tensor(Ellipsoid(1.0), C₀)

# Prolate spheroid
hill_tensor(Ellipsoid(3.0, 1.0, 1.0), C₀)

# 2D ellipse
hill_tensor(Ellipsoid(1.0, 0.5), TensISO{2}(3*(λ+2μ/3), 2μ))
```

## Degenerate limits

When an `Ellipsoid` constructor receives a real semi-axis equal to
`Inf` or `0`, it returns the appropriate dedicated type:

| Call | Returned type | See |
| --- | --- | --- |
| `Ellipsoid(Inf, b, c)` with `b, c > 0` | `Cylinder` | [cylindrical inclusions](cylindrical_inclusions.md) |
| `Ellipsoid(a, b, 0)` with `a, b > 0` | `EllipticCrack` | [cracks](cracks.md) |
| `Ellipsoid(Inf, b, 0)` with `b > 0` | `RibbonCrack` | [cracks](cracks.md) |
| `Ellipsoid(Inf, Inf, c)` | `ArgumentError` (slab, out of scope) | |
| `Ellipsoid(a, 0, 0)` | `ArgumentError` (needle, out of scope) | |

The detection is active only for real element types; with symbolic
types (`SymPy.Sym`, `Symbolics.Num`) call the dedicated constructor
(`Cylinder`, `EllipticCrack`, `RibbonCrack`) explicitly.
