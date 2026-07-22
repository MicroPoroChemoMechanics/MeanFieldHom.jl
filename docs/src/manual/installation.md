# Installation

`MeanFieldHom` is not (yet) in the General registry; install it directly from
GitHub:

```julia
using Pkg
Pkg.add(url = "https://github.com/MicroPoroChemoMechanics/MeanFieldHom.jl")
```

Its dependencies (`TensND.jl`, `OrdinaryDiffEq.jl`, `Elliptic.jl`, `QuadGK.jl`,
`Polynomials.jl`, `PolynomialRoots.jl`, `Tensors.jl`, …) come from the Julia
General registry. The `DECUHR` cubature backend (`import DECUHR, Integrals`)
and `SymPy` symbolic closed forms are optional package extensions. Type-generic
elliptic integrals are bundled internally as the
[`MeanFieldHom.Elliptic`](@ref MeanFieldHom.Elliptic) submodule.

For development from a clone of the repository, instantiate the project
before first use:

```shell
cd /path/to/MeanFieldHom.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Run the test suite with:

```shell
julia --project=. -e 'using Pkg; Pkg.test()'
```
