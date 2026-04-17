# Installation

`MeanFieldHom` depends on three Julia packages shipped alongside it in
the ECHOES repository (`DECUHR.jl`, `GenericElliptic.jl`, and a
registered `TensND.jl`).  Instantiate the project before first use:

```shell
cd interface/julia/MeanFieldHom.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Run the test suite with:

```shell
julia --project=. -e 'using Pkg; Pkg.test()'
```
