<p>
  <img src="./docs/src/assets/logo.svg" width="100">
</p>

# MeanFieldHom

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://MicroPoroChemoMechanics.github.io/MeanFieldHom.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://MicroPoroChemoMechanics.github.io/MeanFieldHom.jl/dev/)

[![CI](https://github.com/MicroPoroChemoMechanics/MeanFieldHom.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/MicroPoroChemoMechanics/MeanFieldHom.jl/actions/workflows/CI.yml)

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/MicroPoroChemoMechanics/MeanFieldHom.jl/blob/main/LICENSE)
[![code style: runic](https://img.shields.io/badge/code_style-%E1%9A%B1%E1%9A%A2%E1%9A%BE%E1%9B%81%E1%9A%B2-pink)](https://github.com/fredrikekre/Runic.jl)

`MeanFieldHom.jl` is a Julia framework for **mean-field homogenization**
of heterogeneous materials. It provides Hill polarisation tensors for
ellipsoidal inclusions, crack-opening-displacement tensors, stress and
displacement intensity factors for flat cracks, and second-order Hill
tensors for transport problems — all under a common abstraction
hierarchy, a shared numerical core, and a central dispatch mechanism.

The package is geared toward prototyping, symbolic simplification
(`SymPy`, `Symbolics`) and forward-mode automatic differentiation
(`ForwardDiff`).

## Features

| Sub-module                 | Responsibility                                                                          |
| -------------------------- | --------------------------------------------------------------------------------------- |
| `MeanFieldHom.Core`        | Abstractions, traits, shared numerics (Green / Newton kernels, Masson residue, DECUHR). |
| `MeanFieldHom.Elasticity`  | Hill polarisation tensor for ellipsoidal inclusions (2D / 3D, iso / aniso).             |
| `MeanFieldHom.Cracks`      | COD tensor, compliance contribution, SIF and DIF for elliptic / ribbon cracks.          |
| `MeanFieldHom.Conductivity`| 2nd-order Hill tensor for transport problems.                                           |
| `MeanFieldHom.Schemes`     | Placeholder for mean-field schemes (dilute, Mori–Tanaka, self-consistent, PCW, …).      |

## Installation

`MeanFieldHom.jl` is a **private** package (public release planned). Clone it
and instantiate its environment:

```shell
git clone https://github.com/MicroPoroChemoMechanics/MeanFieldHom.jl
julia --project=MeanFieldHom.jl -e 'using Pkg; Pkg.instantiate()'
```

Its MPCM dependencies
[`DECUHR.jl`](https://github.com/MicroPoroChemoMechanics/DECUHR.jl) (adaptive
cubature backend) and `TensND.jl` (structured tensors) are pinned to their
GitHub repositories via `[sources]`; all other dependencies (`Integrals.jl`,
`OrdinaryDiffEq.jl`, `Elliptic.jl`, `Polynomials.jl`, `PolynomialRoots.jl`,
`QuadGK.jl`, `Tensors.jl`, …) come from the Julia General registry.

Type-generic elliptic integrals are bundled as the
`MeanFieldHom.Elliptic` submodule.

## Quick start

```julia
using MeanFieldHom, TensND

# Isotropic matrix, E = 210 GPa, ν = 0.3
E, ν = 210e3, 0.3
λ = E*ν/((1+ν)*(1-2ν));  μ = E/(2*(1+ν))
C₀ = TensISO{3}(3*(λ+2μ/3), 2μ)

# Hill polarisation for a sphere
P = hill_tensor(Ellipsoid(1.0), C₀)

# Crack opening displacement for a penny-shaped crack
B = cod_tensor(PennyCrack(1.0), C₀)

# Conductivity — second-order Hill tensor
K₀ = TensISO{3}(5.0)
P_cond = hill_tensor(Ellipsoid(1.0), K₀)
```

Every entry point accepts `method = :auto | :residues | :decuhr` and
the keyword tuple `(abstol, reltol, maxiters)`; see the in-line
docstrings for details.

## Tests

```shell
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Documentation

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

The package accompanies:

> Barthélémy, J.-F. (2026). *Stress Intensity Factors in Anisotropic
> Media.*

## Credits and Acknowledgements

Developed by [Jean-François Barthélémy](https://github.com/jfbarthelemy),
researcher at [Cerema](https://www.cerema.fr/en) in the research team
[UMR MCD](https://mcd.univ-gustave-eiffel.fr/).

Parts of this codebase were developed with the assistance of Anthropic's
*Claude Code*, under the author's review and numerical validation.
