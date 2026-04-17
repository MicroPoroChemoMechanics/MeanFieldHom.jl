# Cracks

```julia
using MeanFieldHom, TensND
E, ν = 210.0, 0.3
k = E/(3*(1-2ν)); μ = E/(2*(1+ν))
C₀ = TensISO{3}(3k, 2μ)

# Penny-shaped crack
B = cod_tensor(PennyCrack(1.0), C₀)

# Ribbon crack
r = RibbonCrack(0.5)
B_r = cod_tensor(r, C₀)
```
