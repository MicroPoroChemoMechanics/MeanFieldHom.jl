# MeanFieldHom.jl

Julia framework for **mean-field homogenization** of heterogeneous materials:
predict the effective elastic, transport and viscoelastic properties of a
microstructure from the properties, shapes, orientations and volume fractions
of its phases.

Given a reference medium and an inclusion shape, `MeanFieldHom` builds the
Hill polarization tensor ``\mathbb{P}`` (Eshelby's result), derives the
localization tensor for each phase, and assembles a **homogenization scheme**
— dilute, Mori–Tanaka, self-consistent, differential, PCW, or the classical
bounds — into an effective stiffness or conductivity. The same machinery
handles flat cracks (opening-displacement and intensity factors), composite
`n`-layer spheres and confocal spheroids with imperfect interfaces, and
ageing linear viscoelasticity, all under one abstraction hierarchy, a shared
numerical core, and forward-mode automatic differentiation throughout.

`MeanFieldHom` is a pure-Julia reimplementation of the Eshelby/Hill machinery
of the [Echoes](https://github.com/jeanfrancoisbarthelemy/echoes) C++/Python
codebase; see [From Echoes to MeanFieldHom](@ref) for the translation guide
and [Theory — reading path](theory/index.md) for the shared conventions and
bibliography.

## Where to start

| If you want to… | Go to |
| :--- | :--- |
| understand the theory before using the code | [Theory](theory/index.md) — the Eshelby/Hill chain, in the order it is built |
| install the package and run the first example | [Installation](manual/installation.md) |
| learn the API by worked example, topic by topic | [Tutorials](tutorials/index.md) |
| see full micromechanical models of real materials | [Applications](applications/transport.md) — cement paste, concrete, bituminous mixtures |
| look up a function's docstring | [API reference](api/elliptic.md) |
| extend the package (new inclusion, algorithm, scheme) | [Developer guide](developer/architecture.md) |

## Sub-modules

| Module | Responsibility |
| :--- | :--- |
| [`MeanFieldHom.Elliptic`](@ref) | Type-generic Legendre and Carlson elliptic integrals (`ForwardDiff`/`Sym` compatible). |
| [`MeanFieldHom.Core`](@ref) | Abstractions, traits, shared numerics (Green/Newton kernels, Masson residue algorithm, DECUHR seam). |
| [`MeanFieldHom.Elasticity`](@ref) | Hill polarization tensor for ellipsoidal inclusions and infinite cylinders (2D/3D, iso/aniso/TI-coaxial). |
| [`MeanFieldHom.Cracks`](@ref) | Crack-opening-displacement (COD) tensor, compliance contribution, stress/displacement intensity factors. |
| [`MeanFieldHom.Conductivity`](@ref) | Second-order Hill tensor for transport problems (diffusion, conduction, Darcy flow), closed form for any matrix anisotropy. |
| `MeanFieldHom.LayeredSpheres` | `n`-layer composite spheres (Hervé–Zaoui, Christensen–Lo), five interface types, volume-average and pointwise localization. |
| `MeanFieldHom.LayeredSpheroids` | `n`-layer confocal spheroids, conduction, Kapitza / surface-conductive interfaces, series or quadrature evaluation. |
| `MeanFieldHom.Schemes` | RVE container and `homogenize`; Voigt, Reuss, Dilute, Mori–Tanaka, Maxwell, PCW, self-consistent, asymmetric SC, differential; exact vs. best-fit symmetrization; `ForwardDiff` sensitivities. |
| `MeanFieldHom.Viscoelasticity` | Ageing linear viscoelasticity via Volterra operators — every scheme, cracks and layered spheres included. |

## Quick example

```julia
using MeanFieldHom, TensND

# Isotropic matrix, bulk/shear moduli (k, μ) = (30, 10) GPa
C₀ = iso_stiffness(30.0, 10.0)

# Hill polarization tensor for a spherical inclusion
P = hill_tensor(Ellipsoid(1.0), C₀)

# A porous material, 10 % spherical voids, homogenized by Mori-Tanaka
rve = RVE(:M)
add_matrix!(rve, Ellipsoid(1.0), Dict(:C => C₀))
add_phase!(rve, :V, Ellipsoid(1.0), Dict(:C => iso_stiffness(0.01, 0.005)); fraction = 0.1)
k_eff, μ_eff = k_mu(homogenize(rve, MoriTanaka(), :C))
```

Every entry point differentiates through `ForwardDiff` and accepts symbolic
(`SymPy`/`Symbolics`) coefficients out of the box — see
[Derivatives and sensitivities](tutorials/sensitivities.md) and
[Symbolic spheres](tutorials/symbolic_spheres.md).
