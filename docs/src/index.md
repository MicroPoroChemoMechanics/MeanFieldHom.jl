# MeanFieldHom.jl

Julia framework for mean-field homogenisation of heterogeneous
materials. Provides Hill polarisation tensors for ellipsoidal
inclusions, crack-opening-displacement tensors and intensity factors
for flat cracks, and second-order Hill tensors for transport problems.
Paves the way for a full homogenisation stack (schemes, RVEs,
viscoelastic laws, user-defined inclusions).

## Sub-modules

- [`MeanFieldHom.Core`](@ref) — abstractions, traits, shared numerics.
- [`MeanFieldHom.Elasticity`](@ref) — Hill polarisation (2D / 3D).
- [`MeanFieldHom.Cracks`](@ref) — COD, SIF, DIF for flat cracks.
- [`MeanFieldHom.Conductivity`](@ref) — 2nd-order Hill tensor.
- [`MeanFieldHom.LayeredSpheres`](@ref) — n-coated-sphere assemblages
  (Hervé–Zaoui, Christensen–Lo) with imperfect interfaces.
- [`MeanFieldHom.Schemes`](@ref) — RVE container, mean-field schemes
  (Voigt, Reuss, Dilute, Mori-Tanaka, Maxwell, PCW, SC, ASC,
  Differential), parameter sensitivities (autodiff).
- [`MeanFieldHom.Viscoelasticity`](@ref) — ageing linear viscoelastic
  homogenisation (`ViscoLaw`, trapezoidal Stieltjes discretisation,
  Volterra algebra, ALV Hill kernel, iso / TI / ortho fast paths,
  cracks, layered spheres, order-2 conductivity / diffusion).

## Quick example

```julia
using MeanFieldHom, TensND
E, ν = 210e3, 0.3
λ = E*ν/((1+ν)*(1-2ν)); μ = E/(2*(1+ν))
C₀ = TensISO{3}(3*(λ+2μ/3), 2μ)
P  = hill_tensor(Ellipsoid(1.0), C₀)
```

See the user manual for more examples and the developer documentation
for guidance on extending the package.
