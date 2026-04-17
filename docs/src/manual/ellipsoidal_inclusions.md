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
