# [Installation](@id man-installation)

`MeanFieldHom` is registered in Julia's General registry.

```julia
julia> import Pkg; Pkg.add("MeanFieldHom")
```

or, from the Pkg REPL mode (`]`):

```julia
pkg> add MeanFieldHom
```

No additional registry is required: its dependencies (`TensND.jl`,
`OrdinaryDiffEq.jl`, `Elliptic.jl`, `QuadGK.jl`, `Polynomials.jl`,
`PolynomialRoots.jl`, `Tensors.jl`, …) come from General as well. The `DECUHR`
cubature backend (`import DECUHR, Integrals`) and `SymPy` symbolic closed forms
are optional package extensions. Type-generic elliptic integrals are bundled
internally as the
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
