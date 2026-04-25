# Changelog

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
