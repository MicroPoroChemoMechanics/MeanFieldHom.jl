# Tutorials

Mean-field homogenization replaces a heterogeneous microstructure — a
matrix carrying inclusions, pores, or cracks — by an equivalent
homogeneous medium with the same overall (elastic, conductive,
viscoelastic, or strength) response. `MeanFieldHom` computes that
response from a handful of ingredients: phase properties, geometries,
volume fractions, and a **scheme** that encodes an assumption about how
the phases interact.

These pages teach that workflow from the ground up, **simplest to most
complex**, with worked examples, the underlying equations, and plots
you can reproduce line by line. The **porous material** — a solid
riddled with pores — recurs throughout and gets two dedicated pages,
because it is both the simplest non-trivial microstructure and the one
where the choice of scheme matters most.

## Reading path

1. [A first homogenization](01_first_estimate.md) — build your first
   RVE and compare the dilute and Mori–Tanaka estimates.
2. [Bounds and classical schemes](02_bounds_and_schemes.md) — Voigt/Reuss
   bounds, self-consistent, and where every scheme sits between them.
3. [Porous materials and the self-consistent trap](03_porous_materials.md)
   — why soft pores break the naive self-consistent iteration.
4. [Porous benchmark: all schemes](04_porous_benchmark.md) — every
   scheme side by side on the canonical porosity sweep, spheres and
   oblate pores alike.
5. [The differential scheme and path dependence](05_differential_paths.md)
   — incremental homogenization and why mixing order matters.
6. [Cracks and crack density](06_cracks.md) — from volume fraction to
   crack density, and the crack-opening-displacement tensor.
7. [Viscoelastic composites](07_viscoelasticity.md) — complex moduli in
   the frequency domain, and a first taste of ageing creep.
8. [Derivatives and sensitivities](08_sensitivities.md) — differentiate
   any homogenization result with `ForwardDiff`, no finite differences.
9. [From derivatives to a strength criterion](09_strength_criteria.md)
   — turn those derivatives into a macroscopic strength criterion for a
   porous solid.
10. [From Echoes to MeanFieldHom](10_from_echoes.md) — translate your
    Echoes (C++/Python) scripts, and cross-check with Echoes live via
    PyCall.

## Prerequisites

These tutorials assume you can install and load the package (see
[Installation](../manual/installation.md)) and focus on *why* each
scheme and computation is used, not on restating the API. For a terse
reference — the full scheme list, aliases, and keyword arguments — see
[Homogenization schemes](../manual/schemes.md) and the other
[Manual](../manual/ellipsoidal_inclusions.md) pages. For the underlying
derivations, see [Theory](../theory/overview.md). Once comfortable
here, the [Applications](../applications/transport.md) section walks
through full case studies (cement paste, transport properties,
strength, bituminous binders, ageing creep).
