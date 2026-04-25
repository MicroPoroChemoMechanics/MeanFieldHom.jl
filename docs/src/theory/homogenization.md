# Homogenisation schemes

`MeanFieldHom.Schemes` provides ten classical mean-field homogenisation
schemes for computing the *effective* property tensor of a multi-phase
heterogeneous medium given (i) the phase geometries, (ii) the phase
properties, and (iii) the phase volume fractions or crack densities.
This page summarises the mathematics behind each scheme and the design
choices that govern their implementation.

## Notation

The Representative Volume Element (RVE) consists of:

- a *matrix* phase of property tensor ``\mathbb C_0`` (or ``\mathbf K_0``
  for the 2nd-order conductivity problem),
- one or more *inclusion* phases of property tensors
  ``\mathbb C_i``, geometries ``\mathcal G_i``, and amounts
  ``f_i`` (volume fraction) or ``\varepsilon_i`` (crack density).

For each inclusion the **dilute strain concentration tensor**
``\mathbb A_\mathrm{dil}^{(i)}`` and the **size-independent stiffness
contribution** ``\mathbb N_i = (\mathbb C_i - \mathbb C_0):
\mathbb A_\mathrm{dil}^{(i)}`` are the natural building blocks
([Kachanov & Sevostianov 2018](@cite kachanov2018)). The dual
**compliance contribution** ``\mathbb H_i = (\mathbb S_i - \mathbb S_0):
\mathbb A_\sigma^{(i)}`` is more natural for cracks (whose stiffness
contribution is the rank-1 limit of a divergent eigenvalue).

## Bounds

| Scheme | Formula |
| --- | --- |
| **Voigt** | ``\langle \mathbb C \rangle = \sum_i f_i \mathbb C_i`` (upper bound, [Hill 1965](@cite hill1965)) |
| **Reuss** | ``\langle \mathbb S \rangle^{-1}`` (lower bound) |

Cracks are ignored in both bounds: their volume contribution vanishes in
the penny limit (`c → 0`) while their density stays finite.

## One-shot schemes (require a matrix)

| Scheme | Effective stiffness |
| --- | --- |
| **Dilute** | ``\mathbb C_0 + \sum_i f_i \mathbb N_i`` (first order in `f`) |
| **DiluteDual** | ``\big(\mathbb S_0 + \sum_i f_i \mathbb H_i\big)^{-1}`` |
| **Mori-Tanaka** | ``\mathbb C_0 + \big(\sum_i f_i \mathbb N_i\big) : \big(f_m\,\mathbb I + \sum_i f_i \mathbb A_\mathrm{dil}^{(i)}\big)^{-1}`` ([Mori-Tanaka 1973](@cite mori1973), [Christensen 1990](@cite christensen1990)) |
| **Maxwell** | ``\mathbb C_0 + \Sigma : (\mathbb I - \mathbb P_d : \Sigma)^{-1}`` with `P_d` the Hill tensor of the *outer distribution shape* |
| **PCW** | identical algebraic form, distribution-shape-aware ensemble interpretation ([Ponte-Castañeda & Willis 1995](@cite ponte1995)) |

The **distribution shape** is stored at the RVE level (default: unit
sphere ⇒ Mori-Tanaka limit). Any `AbstractInclusion` can be used; the
hierarchy [`AbstractDistributionShape`](@ref) leaves room for a future
`PairwiseDistribution` extension following [Willis 1982](@cite willis1982).

## Iterative schemes

| Scheme | Iteration |
| --- | --- |
| **SelfConsistent** ([McLaughlin 1977](@cite mclaughlin1977)) | ``\mathbb C^{(n+1)} = \big(\sum_i f_i \mathbb C_i : \mathbb A_\mathrm{dil}^{(i)}(\mathbb C^{(n)})\big) : \big(\sum_i f_i \mathbb A_\mathrm{dil}^{(i)}(\mathbb C^{(n)})\big)^{-1}`` |
| **AsymmetricSelfConsistent** | switches between stiffness- and compliance-form iteration based on the matrix-vs-Voigt-bound contrast |

The default solver is a damped Picard fixed point (Anderson with memory
1, Dual-safe). Loading `NonlinearSolve.jl` activates the
`MeanFieldHomNonlinearSolveExt` extension, which accepts every SciML
non-linear algorithm (`NewtonRaphson()`, `TrustRegion()`,
`Anderson()`, …) via the `algorithm` keyword of [`SelfConsistent`](@ref).

## Differential scheme

The **DifferentialScheme** integrates the Norris ODE
([Norris 1985](@cite norris1985))
```math
\frac{\mathrm d \mathbb C}{\mathrm d f_i} = (\mathbb C_i - \mathbb C):\mathbb A_\mathrm{dil}^{(i)}(\mathbb C)
```
along a user-selectable trajectory. At step `k` the per-phase increment
``\Phi_k`` is the solution of a small linear system that compensates for
the matrix already homogenised at step ``k-1``:

```math
\big(\mathbb I - \mathrm{diag}(\mathrm{path}_j(k-1)\,f_j)\big)\,\Phi_k
  = (\mathrm{path}_i(k) - \mathrm{path}_i(k-1))\,f_i .
```

Three trajectories are shipped:

- **Proportional** (default) — every phase grows at the same relative
  rate `k/N`; target fractions reached simultaneously.
- **Sequential(order)** — phases are introduced one after another; each
  ramps from 0 to its target fraction over its allotted slice of steps.
- **CustomPath(dict)** — explicit per-phase monotone trajectory with
  `path_i(0) = 0` and `path_i(N) = 1`.

For multi-phase RVEs the three trajectories agree in the dilute limit
(`f → 0`) and diverge by an amount proportional to `f` at finite
fractions — a *physical* feature of the differential scheme, not a
numerical artefact. Cracks (`CrackDensity`) are integrated separately
since they don't compete for matrix volume.

## Number-type compatibility

Every scheme is mandated to support:

- `Float64` — default;
- `ForwardDiff.Dual` — sensitivity analysis through fractions, moduli,
  geometric parameters;
- `Complex{Float64}` — frequency-domain viscoelasticity (parity with
  the C++ ECHOES library, templated on `T = double | complex<double>`);
- `SymPy.Sym`, `Symbolics.Num`, `BigFloat` — best-effort, with explicit
  documentation of any limitation (the iterative SC solvers are not
  symbolic-friendly because the linear-system Jacobian must be numeric).
