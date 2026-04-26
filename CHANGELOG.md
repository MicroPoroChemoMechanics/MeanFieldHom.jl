# Changelog

## v0.4.0 — Friendly autodiff sensitivities, RVE-level symmetrize, Hill-symmetric SC

A small but expressive API exposing `ForwardDiff`-based derivatives of any
homogenisation result with respect to any scalar input parameter — physical
(stiffness coefficient, conductivity), geometric (radii, semi-axes, crack
opening, distribution-shape envelope) or volume-fraction / crack-density —
*and* for arbitrary scalar fields of inclusion types defined later by the
user. The autodiff path unlocks geometric and user-type sensitivities that
were not previously practical, and the multi-scale chain rule is taken care
of automatically by composing several `homogenize` calls inside a single
closure.

The release also ships an RVE-level orientation-distribution projection
(`symmetrize`), a corrected Hill-symmetric self-consistent scheme that
percolates exactly at φ=0.5 for spherical pores, a `Spheroid` convenience
constructor, Dual-stable SC convergence, and a `select_best` mode that
mirrors the C++ reference's behaviour at percolation thresholds.

### Additions

- **Lens hierarchy** `AbstractParameter` with four concrete kinds
  (`AmountParameter`, `PropertyParameter`, `GeometryParameter`,
  `DistributionShapeParameter`) plus user-friendly helpers (`amount`,
  `property`, `geometry`, `shape_param`).
- **`get_param(rve, p)` / `set_param(rve, p, value)`** — read / immutable
  update of the scalar designated by a lens, with automatic eltype
  promotion to integrate `ForwardDiff.Dual` cleanly.
- **Public autodiff entry points** `derivative`, `gradient`, `jacobian`
  and the closure fallback `sensitivity(f, x₀)`. They become available
  after `using ForwardDiff` (weak extension `MeanFieldHomForwardDiffExt`).
- **Generic geometry-field reflection** `_replace_geom_field` based on
  `@generated` reconstruction with uniform sibling-field eltype
  promotion. User-defined inclusions whose constructor follows the
  parametric Julia auto-generated pattern (`MyType{T,B}(args...)`) are
  differentiable through their scalar fields without any library change.
- **Symbol selectors for property tensors** mapping named coefficients
  (`:bulk`, `:shear`, `:transverse`, `:axial`, `:ℓ₁`..`:ℓ₆`) to the
  positional indices of `get_data(tensor)` for `TensISO{2}`,
  `TensISO{4,3}`, `TensTI{2}` and `TensTI{4}`.
- **`MeanFieldHomForwardDiffExt`** weak extension activating the public
  API on `using ForwardDiff`. ForwardDiff is registered in `[weakdeps]`
  alongside NonlinearSolve and SymPy; no new hard dependency.
- **RVE-level orientation symmetrize** via the `symmetrize` keyword on
  `add_matrix!` / `add_phase!`. Three options:
  - `:none` (default): no projection.
  - `:iso`: Reynolds average over `SO(3)` ⇒ isotropic contribution
    (`TensISO`).
  - `:ti` / `TISymmetrize(axis)`: Reynolds average over rotations around
    `axis` ⇒ transversely-isotropic contribution (`TensTI(axis)`).
  Implemented for tensor orders 2 and 4. The TI projection currently
  routes the matrix through an iso projection during the
  localisation-tensor computation (workaround for non-coaxial inclusion
  families); see [`src/Schemes/symmetrize.jl`](src/Schemes/symmetrize.jl)
  for the rationale.
- **`Spheroid(ω; euler_angles)`** convenience constructor on top of
  `Ellipsoid`, mirroring the `spheroidal(omega)` helper of the C++
  reference: `ω = c/a` with one polar semi-axis equal to `ω` and two
  equatorial ones equal to `1`. Eshelby/Hill computations are
  scale-invariant so only the aspect ratio matters.
- **`select_best` keyword on the SC fixed-point solver** — when `true`,
  the solver tracks the best iterate seen during Picard iteration
  (smallest residual on the value field) and returns it at the end.
  Useful for high-contrast iterations that oscillate around the fixed
  point near percolation thresholds; matches the C++ reference's
  `select_best=True` mode.

### Fixes

- **Hill-symmetric self-consistent**: every phase now contributes a
  non-trivial dilute concentration `A_α = inv(I + P(C_α − C_eff))`
  computed in the iterating effective medium, including the matrix
  phase. The previous SC step treated the matrix as having `A = I`
  (Mori-Tanaka-style), which gave the upper SC branch only and
  misplaced the porous-sphere percolation threshold. With the fix,
  porous spheres percolate exactly at φ = 0.5.
- **Dual-stable SC convergence criterion** — the Picard convergence
  test now requires both the value AND every partial of the residual to
  fall below `abstol`. Without that, the value can converge while the
  partials carry residual error of order `‖∂step/∂x‖ × abstol`,
  producing numerically wrong sensitivities through the SC fixed point.
- **TI symmetrize Walpole normalisation**: the `_apply_symmetrize` for
  `TISymmetrize` now divides the W₅ and W₆ projection coefficients by
  `‖W_k‖² = 2`, matching the basis-decomposition convention of
  `TensND.TensTI{4}`. Round-trip on a coaxial TI(ez) tensor is now
  exact.

### Documentation

New manual page `manual/sensitivities.md` (motivation, lens API, closure
fallback, user-inclusion tutorial, multi-scale chain-rule example, and a
section on the `symmetrize` keyword) and auto-API page
`api/sensitivities.md`. Both wired into `docs/make.jl`.

### Scripts

- `scripts/26_sensitivities.jl` — tour of the API (lenses + gradient +
  jacobian + cross-check vs the Christensen 1990 closed form for
  `∂k_MT/∂f`, agreement to ~1e-16).
- `scripts/27_user_inclusion_sensitivity.jl` — extensibility demo on a
  user-defined inclusion type `MyBlob{T,B}` with two numeric fields
  (`radius`, `eccentricity`).
- `scripts/28_multiscale_strength.jl` — three-scale upscaling of
  cement-paste / mortar elasticity and quasi-brittle compression
  strength following Pichler & Hellmich 2011 (SC + MT + MT). The single
  iso hydrate phase + global-μ autodiff approximation matches the
  effective moduli (k, μ, E) of the reference Python implementation to
  rtol ≈ 1e-3 across the (wc, α) grid.
- `scripts/29_porous_schemes.jl` — porous benchmark across all ten
  schemes (sphere and oblate ω = 0.2 with iso symmetrize). After the
  Hill-symmetric SC fix, spherical-pore SC percolates exactly at φ=0.5.
- `scripts/bench_echoes/benchmark_porous.jl`, `benchmark_pichler.jl`
  — PyCall cross-validation against the C++ reference; ten schemes ×
  two cases (sphere / oblate) for porous, six wc curves × twelve α
  points for Pichler. The moduli match the reference to rtol_mod ≈ 1e-3
  across both benchmarks.

### Tests

Three new cross-cutting test files:

- `test_parameters.jl` (round-trip, type-promotion, no-mutation
  invariants on every lens kind),
- `test_sensitivities.jl` (FD vs autodiff cross-check on every scheme,
  closed-form Christensen 1990 match to `~1e-12`, closure fallback,
  `MyBlob` user-inclusion demonstration),
- `test_symmetrize.jl` (round-trip iso/TI projections on 2nd- and
  4th-order tensors, TI(ez) coaxial preservation, integration with
  `homogenize`).

Total: 3421 tests pass.

### Breaking changes

- **SC results differ for systems near percolation** because of the
  Hill-symmetric SC fix. The pre-v0.4 SC step treated the matrix as
  Mori-Tanaka-style (A = I) and therefore selected the upper branch
  unconditionally. The new step is the textbook Hill / Budiansky 1965
  symmetric SC. Users who relied on the old upper-branch behaviour for
  porous-sphere systems should switch to `MoriTanaka` (which is also
  not broken by the fix).
- **`homogenize` API** keeps the `homogenize(rve, scheme; property=:C)`
  kwarg form for backward compatibility but the recommended signature
  is now `homogenize(rve, scheme, property::Symbol)` with `property`
  required and positional.

### Notes

- `Complex{T}` autodiff is not supported (ForwardDiff does not mix Dual
  and Complex cleanly). Symbolic differentiation goes through SymPy on
  the closed-form schemes directly (already supported).
- The `AsymmetricSelfConsistent` scheme follows the symmetric-SC fixed
  point in compliance space when matrix is stiff. The C++ reference's
  ASC uses a different formulation (compliance-form Mori-Tanaka with
  iterating reference) which converges to a different branch for
  porous oblate systems away from percolation; this is documented in
  `scripts/bench_echoes/benchmark_porous.jl`.

## v0.3.0 — RVE container + 10 homogenisation schemes

New `MeanFieldHom.Schemes` sub-module: a Representative Volume Element
container plus the ten classical mean-field homogenisation schemes ported
from C++ ECHOES, with a few Julia-idiomatic improvements.

### Additions

- **`RVE`** container with `add_matrix!`, `add_phase!`, helpers
  (`matrix_phase`, `inclusion_phase_names`, `phase_property`,
  `volume_fraction`, `crack_density`, `matrix_volume_fraction`,
  `validate_rve`).  Volume fractions are stored at the RVE level rather
  than on the inclusions — a single inclusion remains usable for
  localisation-tensor calculations without any RVE machinery.
- **`AbstractAmount`** hierarchy with `VolumeFraction` (solid inclusions)
  and `CrackDensity` (flat cracks); the matrix amount is implicit
  (`1 − Σ f_inc`) and crack densities are excluded from that sum.
- **`AbstractDistributionShape`** hierarchy with `UniformDistribution`
  (single outer envelope, current behaviour); leaves an extension hook
  for a future `PairwiseDistribution` (Willis 1982) without breaking
  the public API.
- **Ten homogenisation schemes**: `Voigt`, `Reuss`, `Dilute`,
  `DiluteDual`, `MoriTanaka`, `Maxwell`, `PonteCastanedaWillis`,
  `SelfConsistent`, `AsymmetricSelfConsistent`, `DifferentialScheme`.
- **`homogenize(rve, scheme; property=:C)`** central entry point.  The
  scheme can be a type instance (`MoriTanaka()`,
  `SelfConsistent(algorithm=NewtonRaphson(), abstol=1e-12)`) or a
  `Symbol` shortcut. Canonical Symbol aliases are lowercase
  (`:mt`, `:sc`, `:diff`, …) for consistency with the algorithm-method
  symbols (`:auto`, `:residues`, `:decuhr`); CamelCase and ECHOES
  upper-case codes (`:MT`, `:DIFF`, …) are kept as backwards-compatible
  aliases.
- **Differential trajectories**: `Proportional` (default), `Sequential`
  (phase-by-phase), `CustomPath` (per-phase explicit trajectory) — all
  validated for monotonicity and boundary conditions.
- **SciML weak extension** `MeanFieldHomNonlinearSolveExt` (activated by
  `using NonlinearSolve`) makes every algorithm of `NonlinearSolve.jl`
  available to `SelfConsistent` / `AsymmetricSelfConsistent` via the
  `algorithm` keyword.
- **Conductivity (`property = :K`)** is supported by every scheme
  through 2nd-order tensor algebra (gradient-gradient localisation,
  resistivity contributions for cracks).

### Number-type compatibility

Every new scheme is fully `ForwardDiff.Dual` and `Complex{Float64}`
compatible (frequency-domain viscoelasticity); symbolic `Sym` / `Num`
work on the closed-form schemes (Voigt, Reuss, Dilute, DiluteDual,
Mori-Tanaka, Maxwell, PCW). The asymmetric SC heuristic uses the
Inf-norm rather than the SVD-based 2-norm so it works seamlessly under
`Dual`.

### Documentation

New theory page `theory/homogenization.md`, manual page
`manual/schemes.md`, API page `api/schemes.md`. Bibliography augmented
with `mori1973`, `christensen1990`, `mclaughlin1977`, `norris1985`,
`ponte1995`, `willis1982`. New scripts `scripts/20_voigt_reuss_bounds.jl`
through `scripts/25_echoes_crosscheck.jl`. The latter cross-validates
Mori-Tanaka against the [Christensen 1990](@cite christensen1990) closed
form (exact match to 6 sig. figs. on bulk and shear at five fractions).

### Tests

Around 270 new tests covering construction, numerical bounds, closed
forms, Dual sensitivity (every scheme), Complex moduli sweep, Symbol
shortcuts, and CustomPath validation. Total 3312 tests passing.

## v0.2.0 — alignment with TensND 0.2 (breaking)

Follow-up to TensND 0.2's API unification. MeanFieldHom is iso-functional —
all outputs are unchanged — but every mention of a TensND symbol now uses
the new snake_case + UPPERCASE-acronym convention.

### Breaking changes

- `TensND.TensWalpole` references (type annotations, dispatch rules,
  constructor calls) now use `TensND.TensTI{4}`.  The struct layout is
  identical so numerical behaviour is unchanged.
- Accessor renames propagated from TensND: `getbasis` → `get_basis`,
  `tensbasis` → `tens_basis`, `invKM` → `inv_KM`, `getdata` → `get_data`,
  `getarray` → `get_array`, `getvar` → `get_var`, `getdim` → `get_dim`,
  `getorder` → `get_order`.
- Predicate renames: `isISO` → `is_ISO`, `isTI` → `is_TI`,
  `isOrtho` → `is_ORTHO`.
- Tensor factory renames in scripts and docs: `tensId2` → `tens_Id2`,
  `tensJ4` → `tens_J4`, `tensTI` → `tens_TI`, etc.

### Additions

None — functional surface unchanged.

### Migration guide

If you have your own code depending on MeanFieldHom dispatch, apply the
same renames as listed in TensND's v0.2 changelog. All MeanFieldHom tests
(2865) pass without behavioural change after migration.
