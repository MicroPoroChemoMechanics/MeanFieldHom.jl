# Conductivity

```julia
using MeanFieldHom, TensND
K₀ = TensISO{3}(5.0)
P  = hill_tensor(Ellipsoid(1.0), K₀)    # 2nd-order Hill tensor
```
