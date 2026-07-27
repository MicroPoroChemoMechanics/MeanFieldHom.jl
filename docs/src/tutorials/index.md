# Tutorials

Mean-field homogenization replaces a heterogeneous microstructure — a matrix
carrying inclusions, pores, or cracks — by an equivalent homogeneous medium with
the same overall response. `MeanFieldHom` computes that response from a handful
of ingredients: phase properties, geometries, volume fractions, and a **scheme**
encoding an assumption about how the phases interact.

These pages teach that workflow with worked examples, the underlying equations,
and plots you can reproduce line by line. They are grouped by theme; within each
group, simplest first.

## Fundamentals

Start here. The **porous material** — a solid riddled with pores — recurs
throughout and gets two pages, because it is both the simplest non-trivial
microstructure and the one where the choice of scheme matters most.

1. [A first homogenization](first_estimate.md) — build your first RVE and
   compare the dilute and Mori–Tanaka estimates.
2. [Bounds and classical schemes](bounds_and_schemes.md) — Voigt/Reuss bounds,
   self-consistent, and where every scheme sits between them.
3. [Porous materials and the self-consistent trap](porous_materials.md) — why
   soft pores break the naive self-consistent iteration.
4. [Porous benchmark: all schemes](porous_benchmark.md) — every scheme side by
   side on the canonical porosity sweep, spheres and oblate pores alike.
5. [The differential scheme and path dependence](differential_paths.md) —
   incremental homogenization and why mixing order matters.

## Inclusions and geometries

Beyond the single homogeneous ellipsoid: flat cracks, and **composite**
inclusions that carry no Hill tensor at all and plug into the schemes through
their volume-averaged concentration tensors instead.

- [Cracks and crack density](cracks.md) — from volume fraction to crack density,
  and the crack-opening-displacement tensor.
- [Layered spheres](generated/layered_sphere.md) — Hervé–Zaoui `n`-layer
  localization and layer averages.
- [Layered spheroids: geometry and effective conductivity](generated/layered_spheroid_effective.md)
  — the confocal `n`-layer spheroid, its API, the effective conductivity of a
  composite reinforced by such particles, the equivalent particle, and the
  numerical accuracy of the harmonic series.
- [Imperfect interfaces: what they do to the local fields](generated/layered_spheroid_interfaces.md)
  — pointwise temperature and flux, streamlines, an animated sweep over the
  interface conductance, and an interactive 3-D view.
- [Highly conducting interfaces](generated/layered_spheroid_hc.md) — the
  equivalent conductivity of an HC-coated particle against aspect ratio.

## Beyond elasticity

- [Viscoelastic composites](viscoelasticity.md) — complex moduli in the
  frequency domain, and a first taste of ageing creep.

## Differentiation and solvers

- [Derivatives and sensitivities](sensitivities.md) — differentiate any
  homogenization result with `ForwardDiff`, no finite differences.
- [From derivatives to a strength criterion](strength_criteria.md) — turn those
  derivatives into a macroscopic strength criterion for a porous solid.
- [Nonlinear solvers for the self-consistent fixed point](nonlinear_solvers.md) —
  solve the fixed point with `NonlinearSolve.jl` instead of Picard iteration, and
  check that sensitivities agree regardless of which solver found it.
- [Symbolic spheres](symbolic_spheres.md) — derive the Eshelby/Hill tensors and
  the dilute, Mori–Tanaka and self-consistent estimates in closed form, then
  substitute numbers at the very end.

## Interoperability and tools

- [From Echoes to MeanFieldHom](from_echoes.md) — translate your Echoes
  (C++/Python) scripts, and cross-check with Echoes live via PyCall.
- [Symmetrization](generated/symmetrization.md) — exact rotation-group average
  versus best-fit projection, and when the two differ.

!!! note "Pages generated from runnable scripts"
    Some pages above live under `tutorials/generated/`: they are produced from
    the demos in `scripts/` by [Literate.jl](https://github.com/fredrikekre/Literate.jl)
    at build time. Each also ships as a Jupyter notebook and as a cleaned
    standalone `.jl` file, and every one of them runs on its own with
    `julia --project=. scripts/NN_name.jl`. Nothing else distinguishes them from
    the hand-written pages.

## Prerequisites

These tutorials assume you can install and load the package (see
[Installation](../manual/installation.md)) and focus on *why* each scheme and
computation is used, not on restating the API. For a terse reference — the full
scheme list, aliases, and keyword arguments — see
[Homogenization schemes](../manual/schemes.md) and the other
[Manual](../manual/ellipsoidal_inclusions.md) pages. For the underlying
derivations, see [Theory](../theory/index.md). Once comfortable here, the
[Applications](../applications/transport.md) section walks through full case
studies drawn from published work — cement paste, transport properties, strength,
bituminous binders, ageing creep.
