# Installation

`MeanFieldHom` depends on one Julia package hosted under the
`MicMacTools` organisation (`DECUHR.jl`) and on the registered
`TensND.jl` (and a handful of registered scientific-computing
dependencies such as `Elliptic.jl`, `QuadGK.jl`, `Polynomials.jl`,
`PolynomialRoots.jl`, `Tensors.jl`). Type-generic elliptic integrals
are bundled internally as the [`MeanFieldHom.Elliptic`](@ref
MeanFieldHom.Elliptic) submodule, with an optional `SymPy` weak
extension for symbolic closed forms. Instantiate the project before
first use:

```shell
cd /path/to/MeanFieldHom.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Run the test suite with:

```shell
julia --project=. -e 'using Pkg; Pkg.test()'
```
