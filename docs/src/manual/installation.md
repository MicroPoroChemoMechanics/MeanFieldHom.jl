# Installation

`MeanFieldHom` is a **private** package (public release planned), so it is not
registered in any public registry. Its MPCM dependencies `DECUHR.jl` (adaptive
cubature backend) and `TensND.jl` (structured tensors) are pinned to their
GitHub repositories via `[sources]`; the remaining dependencies (`Integrals.jl`,
`OrdinaryDiffEq.jl`, `Elliptic.jl`, `QuadGK.jl`, `Polynomials.jl`,
`PolynomialRoots.jl`, `Tensors.jl`) come from the Julia General registry.
Type-generic elliptic integrals are bundled internally as the
[`MeanFieldHom.Elliptic`](@ref MeanFieldHom.Elliptic) submodule, with an
optional `SymPy` weak extension for symbolic closed forms.

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
