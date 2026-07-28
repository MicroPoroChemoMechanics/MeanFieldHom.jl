# Performance benchmarks

`scripts/bench/` holds a benchmark suite designed to make performance claims
*measurements* rather than tables of numbers: a registry of cases, three
independent measurement channels, a committed baseline, a calibrated noise
floor and a numerical gate. This page describes what it measures and reports
the results it has produced. For the reasoning behind the hot paths themselves,
see [Performance notes](performance_notes.md).

## Why not `@elapsed`

The suite replaces a hand-rolled convention that could not support the work,
for three reasons:

1. **It cannot resolve the microsecond tier.** `@elapsed f()` measures exactly
   one call; the acoustic-tensor inverse and a `TensISO ⊡ TensISO` run in
   10–200 ns, below clock granularity. `BenchmarkTools` auto-tunes the number
   of evaluations per sample so such a kernel is measured against a window
   orders of magnitude above the timer's resolution.
2. **It reported allocations from a different call than the one it timed** —
   `@allocated f()` then `@elapsed f()`, two separate invocations, paired in
   the output. That pairing is meaningless for an allocation investigation.
3. **It reported only the minimum**, while one failure mode to guard against
   is precisely a min/total confusion.

`@benchmark` alone is not enough either: it wants `$`-interpolated arguments to
avoid measuring global lookup, but interpolating a pre-built closure **hoists
the capture out of the measured region** — which is exactly what would hide a
closure-boxing regression.

## The three channels

| Channel | How | Answers |
| :--- | :--- | :--- |
| `time` | `BenchmarkTools.@benchmark`, min + median, auto-tuned `evals` | how fast |
| `alloc` | hand-rolled loop of `@allocated` on the **verbatim** thunk, no `$`, ≥2 warm-ups + N samples in the *same* process, **min and max** | how much garbage, and is it stable |
| `counters` | instrumented `Ref{Int}` counters read around a *separate* clean call | how much **work** |

Reporting allocation max as well as min is the cheap upgrade that catches
instability: `max ≠ min` signals a type instability and immediately invalidates
any "one measurement on a fresh process" reading. The `counters` channel
distinguishes *faster* from *did less work* — a drop in adaptive-quadrature
node count is a change of behaviour (and of accuracy), not a speed-up.

## The noise floor and the numerical gate

`--repeat-suite=2` runs the whole registry twice and computes, over the
**control cases only**, `noise = p90(|Δt|/t)`. A case is reported as MOVED only
beyond `max(3·noise, 3 %)`; a *control* case moving that far invalidates the
whole run. Controls are chosen so that no planned change can touch them:
analytical Hill branches, the elliptic integrals, Voigt/Reuss, and the schemes
that call only one of the de-duplicated helpers.

Every case also carries a `checksum` closure **evaluated on the same call that
is timed**, so non-regression is asserted on the real code path:

* `--gate=bitwise` — sha256 of the canonical `%.17g` rendering, a true
  bit-identity test whose stored values are exact and re-diffable;
* `--gate=1e-14` — scale-relative bound
  `max|new−ref| ≤ tol · max(1, max|ref|)`. The relative form is required:
  component magnitudes span ``\mathbb C \sim O(200)``,
  ``\mathbb A_{\varepsilon\varepsilon} \sim O(1)`` and a stiff-matrix
  ``\mathbb H \sim O(10^{-3})``.

Iterative schemes are pinned (`abstol = reltol = 1e-10`, `maxiters = 300`,
`select_best = true`) so the iteration count is part of the contract via the
`sc_iterations` counter.

## The registry

67 cases in six groups:

| Group | Cases | What it probes |
| :--- | ---: | :--- |
| `kernels/` | 11 | Hill and COD tensors per back-end (residues, nested QuadGK, DECUHR), per geometry, plus the per-node acoustic inverse |
| `control/` | 11 | analytical Hill branches, elliptic integrals, and the six schemes that touch only one shared helper — the calibration set |
| `schemes/` | 12 | Mori-Tanaka and self-consistent variants: isotropic, porous, anisotropic matrix, cracked, orientation-binned, conductivity |
| `sens/` | 3 | `ForwardDiff` sensitivities, including `Dual` through nested QuadGK |
| `alv/` | 7 | ageing linear viscoelasticity: trapezoidal assembly, Volterra inversion, full `homogenize_alv` at 50 and 100 steps |
| `tensnd/` | 23 | `TensND` primitives: `getindex`, `get_array`, `collect`, `⊡`, `inv`, symmetry projections, Mandel conversions |

`scripts/bench/baseline.json` is the only versioned report; `results/` is
gitignored.

## Instrumentation is free

The work counters would normally perturb the innermost loop. Measured against
an uninstrumented `git worktree`, same workloads, same machine:

| Case | Time (before → after) | Allocations |
| :--- | :--- | :--- |
| `hill_tensor` triaxial/triclinic, `:nestedquadgk` | 78.601 → 78.572 ms | 103 363 248 B → identical |
| `hill_tensor` triaxial/triclinic, `:residues` | 4.278 → 4.296 ms | 7 377 152 B → identical |
| `cod_tensor` penny/triclinic, `:residues` | 2.729 → 2.716 ms | 3 577 360 → 3 577 280 B (−0.002 %) |
| `hill_tensor` iso sphere (analytical) | 0 → 0 | 0 → 0 |

The timed pass asserts the counters are off, so it can never see the wrapper.

## Results of the optimisation campaign

Run `--gate=1e-14 --repeat-suite=2`, 67 cases, idle machine: **24 moved, 0 gate
failures, noise floor (controls, p90 of |Δt|/t) = 0.8 %**.

| Case | Time | Allocations |
| :--- | ---: | ---: |
| `tensnd/getindex.ortho` | −99.7 % | −95.8 % |
| `tensnd/dcontract.ortho_ortho` | −99.3 % | −94.9 % |
| `tensnd/dcontract.iso_ortho` | −99.3 % | −95.6 % |
| `tensnd/collect.ortho` | −94.4 % | −98.6 % |
| `tensnd/inv_KM.gen` | −87.0 % | +0.0 % |
| `kernels/cod.nqgk.ellipse03.tri` | −85.3 % | −99.4 % |
| `kernels/hill.decuhr.tri.321` | −62.0 % | −80.7 % |
| `schemes/mt.aniso_matrix` | −50.8 % | −50.0 % |
| `schemes/mt.porous.oblate.isosym` | −50.2 % | −14.2 % |
| `schemes/mt.crack.penny.tri` | −49.7 % | −50.0 % |
| `schemes/mt.crack.penny` | −49.0 % | −35.3 % |
| `tensnd/get_array.ortho` | −44.4 % | +0.0 % |
| `kernels/hill2.aniso` | −35.2 % | −18.5 % |
| `kernels/hill.dual.nqgk.tri` | −21.5 % | −0.0 % |
| `schemes/mt.theta_binned_ti.n20` | −17.3 % | −7.2 % |

The four largest wins came from four different causes: sharing the single
expensive Hill/COD solve between the two quantities the schemes always request
together (the `mt.*` row block, with the work counters confirming
40 → 20 Hill solves on the binned case); returning `StaticArrays` from the
innermost COD loop instead of writing through heap 3×3 temporaries; replacing
16 scalar quadratures by one vector-valued one in the anisotropic 2D Hill; and
closed-form `TensOrtho` indexing and contraction in `TensND`.

**Correctness.** 63 of the 67 cases stayed **bit-identical** (`0.0e+00`). The
four that moved did so by pure floating-point reassociation —
`cod.nqgk.ellipse03.tri` at 6.9e-16, `hill.decuhr.tri.321` at 1.7e-18,
`dcontract.iso_ortho` at 4.0e-17, `dcontract.ortho_ortho` at 1.9e-17 — all far
below the 1e-14 tolerance. Those four still fail a strict `--gate=bitwise` run
against the committed baseline, which predates the campaign; that is expected,
not a regression.

**What went up.** One deterministic item: `schemes/mt.conductivity.iso2`,
+22.2 % allocation — the price of the bundled-helper tuple, on the cheapest
case of the suite (528 B in total).

!!! warning "A control case may legitimately move"
    `control/alv.voigt.n50` came out at −5.6 %, i.e. *faster*; the harness
    flags any control deviation without looking at the sign. Verified rather
    than assumed: reproducible over five fresh processes (−4.8 to −6.8 %),
    allocations identical to the byte, checksum bit-identical, work counters
    unchanged. The control set was chosen assuming the changes would not touch
    shared primitives, but `inv_KM`, `tensor_or_array` and basis comparison are
    global — so a control can move for a real reason.

## A worked example: heterogeneous RVE amounts

A later change made `RVE.amounts` heterogeneous, so that a `Dual` or complex
volume fraction no longer needs declaring. It is a good illustration of the
method, because the first measurement pointed at the wrong thing twice.

**Against the committed baseline**, the change looked catastrophic on the cheap
schemes — `voigt.iso2` +75 %, `dilute_dual.iso2` +99 %, `mt.iso2.sphere`
+78 % — and innocuous elsewhere. But the baseline predates the optimisation
campaign, so those numbers mixed the change with everything the campaign had
already gained. A paired A/B against a worktree at the current commit put
`mt.iso2.sphere` back inside the noise: **+59 % was an artefact of the
reference, not an effect of the change.**

**The allocation channel then contradicted the obvious explanation.**
Allocations were *flat or lower*, so the cost was not the boxing one would
assume. Direct measurement on both trees, identified by `pathof`, localised it:

| | previous commit | modified |
| :--- | ---: | ---: |
| `amount_value(rve.amounts[:I])` | 24.7 ns / 16 B | 26.9 ns / 16 B |
| `matrix_volume_fraction` | 12.4 ns / 0 B | 61 ns / 48 B |

The per-phase access was **already** dynamic before the change — the dict's
value type had always been abstract. What had been lost was the type
*refinement* that made `a isa VolumeFraction` land on a concrete
`VolumeFraction{Float64}`, plus a recomputed sum in
`matrix_volume_fraction`. Caching the matrix fraction in a field and putting
the per-phase product behind a barrier that dispatches on the amount's
concrete type brought every case back to or below the previous commit, with
**all 67 checksums bit-identical** and allocations down on 21 cases
(−19 152 B).

Applying the same barrier to the crack `delta_*` paths, where the amount is a
trailing argument and the wrapper needs varargs, *cost* ~1.4 KB per call on
`mt.crack.penny` for no gain — those cases run at tens of microseconds. It was
reverted. The lesson is the harness's whole point: **three channels, a paired
reference, and a willingness to drop the half of a change that does not
measure.**

## A pitfall the harness exists to catch

Earlier in the project, replacing a `Matrix` + `inv` by twelve local scalars
inside a nested quadrature closure made the per-node computation ~10× cheaper
in isolation and **increased total allocation by 24 %** (10.22 MB → 12.69 MB),
because a deeply nested closure boxed differently depending on the number of
its own locals. The fix was to extract the computation into a top-level
function with no captures (→ 8.34 MB). A single `@allocated` on a fresh
process reported the *opposite* of the truth. See
[Performance notes](performance_notes.md) for the general form of this trap.

## Running

```shell
julia --project=scripts/bench -e 'using Pkg; Pkg.instantiate()'

# capture a reference (refuses a dirty tree unless --force)
JULIA_NUM_THREADS=1 julia --project=scripts/bench scripts/bench/bench_suite.jl \
    --record-baseline --label=my-baseline --repeat-suite=2

# after a change
JULIA_NUM_THREADS=1 julia --project=scripts/bench scripts/bench/bench_suite.jl \
    --label=my-change --out=scripts/bench/results/mine.json \
    --baseline=scripts/bench/baseline.json --gate=bitwise --repeat-suite=2

# correctness only, no timing (fast)
julia --project=scripts/bench scripts/bench/bench_suite.jl --verify-only
```

Exit codes: `0` clean, `1` gate failure, `2` control regression, `3` setup
failure.

!!! tip "Comparing against a commit, not against the committed baseline"
    The baseline is a snapshot and drifts away from `main`. To attribute a
    change to *your* edit, measure both sides yourself: `git worktree add
    ../pkg.ref HEAD`, run the suite in each tree back to back, and diff the two
    JSON reports. The relative `[sources]` paths in
    `scripts/bench/Project.toml` resolve correctly as long as the worktree is a
    sibling of the package directory.
