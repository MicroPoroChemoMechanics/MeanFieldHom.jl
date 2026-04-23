# API — Conductivity

The 2nd-order Hill tensor is obtained through the same entry point
[`hill_tensor`](@ref) — this module only registers additional
`_kernel` methods.

The 2nd-order Eshelby tensor ``\mathbf s = \mathbf P \cdot \mathbf K_0``
is likewise obtained via the dispatching
[`eshelby_tensor`](@ref MeanFieldHom.Core.eshelby_tensor) wrapper on a
2nd-order `K₀`.

```@docs
MeanFieldHom.Conductivity
```
