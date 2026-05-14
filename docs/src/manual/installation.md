# Installation

`MeanFieldHom` is released through the dedicated
[MPCM-Registry](https://codeberg.org/MicroPoroChemoMechanics/MPCM-Registry)
on Codeberg. Its sister package `DECUHR.jl` (adaptive cubature backend)
lives in the same registry; `TensND.jl` (structured tensors) is also
in MPCM-Registry, and a handful of registered scientific-computing
dependencies (`Elliptic.jl`, `QuadGK.jl`, `Polynomials.jl`,
`PolynomialRoots.jl`, `Tensors.jl`) come from the Julia General
Registry. Type-generic elliptic integrals are bundled internally as
the [`MeanFieldHom.Elliptic`](@ref MeanFieldHom.Elliptic) submodule,
with an optional `SymPy` weak extension for symbolic closed forms.

Add the MPCM-Registry once, then install the package:

```julia
julia> using Pkg
pkg> registry add https://codeberg.org/MicroPoroChemoMechanics/MPCM-Registry
pkg> add MeanFieldHom
```

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
