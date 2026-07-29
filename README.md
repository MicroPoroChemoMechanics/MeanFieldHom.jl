<p align="center">
  <img src="./docs/src/assets/logo.svg" alt="MeanFieldHom.jl" width="100">
</p>

# MeanFieldHom

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://MicroPoroChemoMechanics.github.io/MeanFieldHom.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://MicroPoroChemoMechanics.github.io/MeanFieldHom.jl/dev/)

[![CI](https://github.com/MicroPoroChemoMechanics/MeanFieldHom.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/MicroPoroChemoMechanics/MeanFieldHom.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/MicroPoroChemoMechanics/MeanFieldHom.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/MicroPoroChemoMechanics/MeanFieldHom.jl)

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/MicroPoroChemoMechanics/MeanFieldHom.jl/blob/main/LICENSE)
[![code style: runic](https://img.shields.io/badge/code_style-%E1%9A%B1%E1%9A%A2%E1%9A%BE%E1%9B%81%E1%9A%B2-pink)](https://github.com/fredrikekre/Runic.jl)

`MeanFieldHom.jl` is a Julia framework for **mean-field homogenization**
of heterogeneous materials: it predicts effective elastic, transport and
viscoelastic properties from the properties, shapes, orientations and
volume fractions of the phases in a microstructure.

It provides Hill polarization tensors for ellipsoidal inclusions and
infinite cylinders (2D/3D, isotropic/anisotropic/TI-coaxial), crack-opening-
displacement tensors with stress and displacement intensity factors for
flat cracks, second-order Hill tensors for transport problems (closed-form
for any matrix anisotropy), composite `n`-layer spheres and confocal
spheroids with imperfect interfaces, ageing linear viscoelasticity, and
the classical mean-field schemes built on top of them (Voigt/Reuss,
dilute, Mori–Tanaka, Maxwell, Ponte Castañeda–Willis, self-consistent,
asymmetric self-consistent, differential) — all under a common
abstraction hierarchy, a shared numerical core, and a central dispatch
mechanism.

The package is geared toward prototyping: forward-mode automatic
differentiation (`ForwardDiff`) and symbolic simplification (`SymPy`,
`Symbolics`) are first-class, not afterthoughts, and every scheme has
`ForwardDiff` sensitivities with respect to fractions, moduli, and
inclusion geometry.

A gallery of full micromechanical models built on the package —
hydrating cement paste, chloride diffusivity, the interfacial transition
zone in concrete, quasi-brittle strength, bituminous mixtures, ageing
creep — lives under [`docs/src/applications/`](docs/src/applications)
and the [Applications](https://MicroPoroChemoMechanics.github.io/MeanFieldHom.jl/stable/applications/transport/)
section of the docs.

## Features

| Sub-module | Responsibility |
| --- | --- |
| `MeanFieldHom.Elliptic` | Type-generic Legendre and Carlson elliptic integrals (`ForwardDiff`/`Sym`-compatible). |
| `MeanFieldHom.Core` | Abstractions, traits, shared numerics (Green / Newton kernels, Masson residue, DECUHR). |
| `MeanFieldHom.Elasticity` | Hill polarization tensor for ellipsoidal inclusions and cylinders (2D / 3D, iso / aniso / TI-coaxial). |
| `MeanFieldHom.Cracks` | COD tensor, compliance contribution, SIF and DIF for elliptic / ribbon cracks. |
| `MeanFieldHom.Conductivity` | 2nd-order Hill tensor for transport problems; closed form for any matrix anisotropy. |
| `MeanFieldHom.LayeredSpheres` | `n`-layer composite spheres, 5 interface types (perfect, spring, membrane, Kapitza, surface-conductive), volume-average and pointwise localization. |
| `MeanFieldHom.LayeredSpheroids` | `n`-layer confocal spheroids, conduction, with Kapitza / surface-conductive interfaces, series or quadrature evaluation. |
| `MeanFieldHom.Schemes` | RVE container and `homogenize`; bounds, dilute, Mori–Tanaka, self-consistent (+ asymmetric), PCW, Maxwell, differential; exact vs. best-fit symmetrization; `ForwardDiff` sensitivities. |
| `MeanFieldHom.Viscoelasticity` | Ageing linear viscoelasticity via Volterra operators, with structured ISO/TI/orthotropic kernel storage — every scheme, cracks and layered spheres included. |

## Installation

`MeanFieldHom.jl` is not (yet) in the General registry; install it directly
from GitHub:

```julia
using Pkg
Pkg.add(url = "https://github.com/MicroPoroChemoMechanics/MeanFieldHom.jl")
```

All dependencies (`TensND.jl`, `OrdinaryDiffEq.jl`, `Elliptic.jl`,
`Polynomials.jl`, `PolynomialRoots.jl`, `QuadGK.jl`, `Tensors.jl`, …) are
resolved from the Julia General registry.

Three package extensions activate on weak dependencies, each optional:

- [`DECUHR.jl`](https://github.com/MicroPoroChemoMechanics/DECUHR.jl) +
  `Integrals.jl` — adaptive cubature backend (`method = :decuhr`); the
  built-in `method = :nestedquadgk` (QuadGK-based, `ForwardDiff`-compatible)
  covers the same cases without it.
- `NonlinearSolve.jl` — lets the self-consistent schemes solve their fixed
  point with any SciML algorithm (`NewtonRaphson`, `TrustRegion`, …)
  instead of the built-in Anderson/Newton iteration, with exact
  `ForwardDiff` sensitivities through the fixed point either way.
- `SymPy.jl` — symbolic closed forms for the elliptic integrals.

Type-generic elliptic integrals themselves are always bundled, as the
`MeanFieldHom.Elliptic` submodule.

## Quick start

```julia
using MeanFieldHom, TensND

# Isotropic matrix, bulk/shear moduli (k, μ) = (30, 10) GPa
C₀ = iso_stiffness(30.0, 10.0)

# Hill polarization for a sphere
P = hill_tensor(Ellipsoid(1.0), C₀)

# Crack opening displacement for a penny-shaped crack
B = cod_tensor(PennyCrack(1.0), C₀)

# Conductivity — second-order Hill tensor
K₀ = TensISO{3}(5.0)
P_cond = hill_tensor(Ellipsoid(1.0), K₀)

# A porous material, 10 % spherical voids, homogenized by Mori-Tanaka
rve = RVE(:M)
add_matrix!(rve, Ellipsoid(1.0), Dict(:C => C₀))
add_phase!(rve, :V, Ellipsoid(1.0), Dict(:C => iso_stiffness(0.01, 0.005)); fraction = 0.1)
k_eff, μ_eff = k_mu(homogenize(rve, MoriTanaka(), :C))
```

Every entry point accepts `method = :auto | :residues | :decuhr` and
the keyword tuple `(abstol, reltol, maxiters)`; every entry point also
differentiates through `ForwardDiff` and accepts symbolic (`SymPy`/
`Symbolics`) coefficients. See the in-line docstrings and the
[Tutorials](https://MicroPoroChemoMechanics.github.io/MeanFieldHom.jl/stable/tutorials/)
for details.

## Tests

```shell
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Documentation

Built with Documenter.jl and deployed at the badges above. Six sections,
roughly in reading order:

| Section | Content |
| --- | --- |
| [Theory](https://MicroPoroChemoMechanics.github.io/MeanFieldHom.jl/stable/theory/) | the Eshelby/Hill chain — polarization tensor → localization → schemes — and its specializations (cracks, layered inclusions, viscoelasticity). |
| [Manual](https://MicroPoroChemoMechanics.github.io/MeanFieldHom.jl/stable/manual/installation/) | installation and a topic-by-topic reference for each inclusion family. |
| [Tutorials](https://MicroPoroChemoMechanics.github.io/MeanFieldHom.jl/stable/tutorials/) | worked examples: bounds and schemes, layered spheres/spheroids, viscoelasticity, sensitivities, symbolic computation. |
| [Applications](https://MicroPoroChemoMechanics.github.io/MeanFieldHom.jl/stable/applications/transport/) | full micromechanical models — cement paste, ITZ concrete, bituminous mixtures, strength, ageing creep. |
| [Developer guide](https://MicroPoroChemoMechanics.github.io/MeanFieldHom.jl/stable/developer/architecture/) | architecture, dispatch, and how to add an inclusion / algorithm / scheme. |
| [API reference](https://MicroPoroChemoMechanics.github.io/MeanFieldHom.jl/stable/api/elliptic/) | every public docstring, grouped by sub-module. |

Build locally:

```shell
julia --project=docs -e 'using Pkg; Pkg.develop(path = "."); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

## License

MIT License — see [LICENSE](LICENSE) for details.

## Citation

See [CITATION.cff](CITATION.cff) for citation details.

**BibTeX entry:**

```bibtex
@software{meanfieldhom_jl,
  author = {Barthélémy, Jean-François},
  title  = {MeanFieldHom.jl: Mean-field homogenization of heterogeneous materials},
  url    = {https://github.com/MicroPoroChemoMechanics/MeanFieldHom.jl},
  year   = {2026}
}
```

## Credits and Acknowledgements

Developed by [Jean-François Barthélémy](https://github.com/jfbarthelemy),
researcher at [Cerema](https://www.cerema.fr/en) in the research team
[UMR MCD](https://mcd.univ-gustave-eiffel.fr/).

Parts of this codebase were developed with the assistance of Anthropic's
*Claude Code*, under the author's review and numerical validation.
