# `scripts/bench_echoes/` — MeanFieldHom.jl vs Echoes (C++) benchmark

Side-by-side numerical + timing comparison between the Julia implementation
and the Echoes C++ reference called from Python via `PyCall.jl`.

## Prerequisites

Echoes must already be importable from a Python interpreter on this
machine (the Boost.Python build works; no recompilation needed here).

```python
# Must succeed in the target interpreter:
>>> import echoes
```

The first time you run this benchmark, point PyCall at that interpreter
and rebuild:

```julia
julia> ENV["PYTHON"] = raw"C:\path\to\python.exe"   # or /usr/bin/python3 on Linux
julia> using Pkg; Pkg.activate("."); Pkg.build("PyCall")
```

## Run

From the package root `MeanFieldHom.jl/`:

```bash
julia --project=scripts/bench_echoes scripts/bench_echoes/benchmark.jl
```

Two summary tables are printed at the end:

- Hill tensor **P** — 5 cases (iso sphere / prolate / oblate, cubic prolate,
  triclinic triaxial).
- Crack compliance **ΔS** (ε = 1) — 3 cases (penny in iso / cubic /
  triclinic matrix).

Each Echoes reference is computed with both `RESIDUES` and `NUMINT3D`
(DECUHR) algorithms.  `RESIDUES` occasionally fails (Masson polynomial
degenerates on perfectly isotropic stiffness) — such cases are reported
as `FAIL` in the summary, and `NUMINT3D` is the robust fallback.

Timings are collected with `BenchmarkTools.@belapsed`.

## Notes

- This subfolder has its **own Project.toml** — it does not pollute the
  main `MeanFieldHom.jl` environment with PyCall / BenchmarkTools.
- `:residues` is used on both sides for the cubic and triclinic cases
  (same mathematical algorithm → expected agreement at machine
  precision).
- For the isotropic cases, the Julia side uses the analytic path
  (`hill_tensor` with `:auto`) and the Python side uses `RESIDUES` +
  `NUMINT3D` — a small discrepancy (~1e-6 to 1e-8) is expected.
