# Changelog

## v0.6.0 — TI ALV fast path, order-2 ALV, BLAS Volterra, ALV cracks roadmap

**TI Walpole-basis fast path** for ALV homogenisation : when every
phase 4-tensor and the matrix kernel are TI 4-tensors with the
**common canonical axis n = e₃** (every 6×6 Mandel block matches the
Walpole structure), `homogenize_alv` now routes through new
`*_alv_ti(ℓ_…)` primitives that operate on **6** `n × n` Volterra
matrices `(ℓ₁, ℓ₂, ℓ₃, ℓ₄, ℓ₅, ℓ₆)` instead of the generic `(6n × 6n)`
block.  The Walpole 2×2 part `[[ℓ₁, ℓ₃]; [ℓ₄, ℓ₂]]` is packed as a
`(2n)×(2n)` block-Volterra matrix and inverted via
`volterra_inverse(_; block_size = 2)`; the two scalar shears
`(ℓ₅, ℓ₆)` go through the LAPACK scalar fast path below.

ISO inputs are subsumed automatically (iso ⊂ TI), so a TI-matrix +
iso-inclusion combination — common in layered concrete /
fibre-reinforced ALV — gets the fast path "for free".  Storage is 6 ·
n² doubles per phase (vs 36 · n² generic), and the inverse cost drops
from `O((6n)²)` to `O((2n)²) + 2·O(n²) ≈ ×3` cheaper than the generic
6n×6n path.

The TI fast path is integrated into all six schemes (Voigt / Reuss /
Dilute / DiluteDual / Mori-Tanaka / Maxwell), via a new
`_try_ti_tuples` helper analogous to the existing `_try_iso_pairs`
detection.

**Order-2 ALV** : new sub-module covering vector-tensor ageing
viscoelasticity (thermal / electrical conductivity, diffusivity,
permittivity).  Operators are stored as `(3n × 3n)` lower-block-
triangular matrices with 3×3 blocks.  Mirrors the order-4 API:

  * `homogenize_alv_order2(rve, scheme, prop; times)` — public entry
    point, dispatching on Voigt / Reuss / Dilute / DiluteDual /
    Mori-Tanaka / Maxwell schemes.
  * `hill_kernel_order2(ell, K_0_law, times)` — Hill polarisation for
    iso ALV matrix + ellipsoidal inclusion (time-space decoupling
    `P̃[block(i,j)] = α₀^{-vol}[i,j] · 𝐈^A`, with `𝐈^A` from the
    existing elastic `tens_IA(ell)`).
  * `iso_order2_params_from_blocks` / `iso_order2_blocks_from_params`
    — per-component parameter extraction (single scalar α for iso
    order-2).
  * `voigt_alv_order2`, `reuss_alv_order2`, `dilute_alv_order2`,
    `dilute_dual_alv_order2`, `mori_tanaka_alv_order2`,
    `maxwell_alv_order2`.

`trapezoidal_matrix(law, times)` now accepts both order-4 (4-tensor /
6×6 Mandel) and order-2 (`TensND.AbstractTens{2,3}` / 3×3 matrix)
sample types, dispatching to the appropriate (B·n)×(B·n) layout.

The order-2 elastic-limit test verifies that ALV reduces to the
existing elastic conductivity `homogenize` to machine precision for
both spherical and spheroidal inclusions.  `scripts/bench_echoes/
bench_order2_alv.jl` provides a template for a Julia–ECHOES
crosscheck on the `fluage_echoes_maxwell_ordre2.py` setup.

**Volterra BLAS / LAPACK fast path** : `volterra_inverse`,
`volterra_left_divide` and `volterra_divide` now dispatch to LAPACK
`trtri` / `trsm` via the `LowerTriangular(...)` wrapper for
`BlasFloat` element types and `n ≥ 64` (the crossover where BLAS
overhead amortises).  Measured speedups vs the hand-rolled forward
substitution: **×9.7** at `n = 500`, **×14.6** at `n = 1000`.  The
hand-rolled fallback is preserved for small grids and for
non-BlasFloat element types (`BigFloat`, `Sym`, `ForwardDiff.Dual`).

**ALV cracks roadmap** : new file
`src/Viscoelasticity/cracks_alv.jl` documents the planned
`cod_kernel_alv` / `compliance_contribution_alv` API, the time-space
decoupling formulas for pure penny / interface-stiffness cracks in
iso ALV matrices, and the integration points with the existing
`CrackDensity` amount in `Schemes`.  Implementation is scheduled for
v0.6.1.

## v0.5.3 — Non-uniform time grid + multi-layer ALV stability + iso fast path

**Bug fix (membrane interface convention)** : when v0.5.2 introduced
the C++-convention σ-form shear M-matrix, the elastic-limit unit test
for `MembraneInterface` started failing because the membrane jump
expressions and the M-matrix used different angular-component
normalisations.  Resolution: keep the M-matrix in the
**Christensen–Lo / SymPy convention** (matching the elastic
state-space recurrence in `LayeredSpheres`) and derive the analytic
`M^{-1}` from the C++ closed form via the diagonal conjugation

```text
M_C++ = D_row · M_Christ–Lo · D_col,
   D_row = diag(1/2, 1, 1/2, 1),  D_col = diag(1, 1, 2, 2)
⇒ M_Christ–Lo^{-1} = D_col · M_C++^{-1} · D_row
```

This gives the **best of both worlds**: numerical stability of the
closed-form (only `(3κ+4μ)^{-vol}` and `μ^{-vol}` `n×n` Volterra
inverses are needed) under the elastic-compatible mode normalisation
where the Christensen–Lo membrane jumps stay correct.  The
`Membrane interface (elastic limit)` test is back to `≤ 1e-10`
tolerance with no convention asterisk.

**Performance — iso-symmetry fast path** : when every phase 4-tensor
and the matrix kernel are iso (`TensISO{4,3}`-valued),
`homogenize_alv` automatically routes through new
`*_alv_iso(αβ_…)` scheme primitives that operate on two scalar
`n × n` Volterra matrices `(α = 3K, β = 2μ)` instead of the generic
`(6n × 6n)` Mandel block matrix.  Theoretical speedups :

  - matrix-matrix product : ~108× cheaper (`216 n³` → `2 n³`)
  - matrix inverse        : ~18× cheaper
  - storage               : 18× smaller

Measured speedup on the script-37 setup (5 phases + Maxwell matrix +
spherical inclusions, MT scheme): **×2.5–7** depending on `n_times`.
The detection happens once per phase via a cheap iso-form pattern
check (`_is_iso_block`); on failure the generic 6n×6n path is
selected automatically — no API change.  New internal helpers
(`iso_schemes_alv.jl`) :
`voigt_alv_iso`, `reuss_alv_iso`, `dilute_alv_iso`,
`dilute_dual_alv_iso`, `mori_tanaka_alv_iso`, `maxwell_alv_iso`,
plus `dilute_concentration_alv_iso`, `dilute_contribution_alv_iso`
used inside `_inclusion_alv_quantities`.

**Bug fix (multi-layer ALV stability)** : on a non-uniform time grid
(e.g. `logspace`) with a matrix relaxation kernel that has multiple
time constants and/or layers with extreme modulus contrast (pores,
step-activated `ViscoLaw`s), the layered-sphere ALV recurrence
diverged from the ECHOES Python reference by 1e-3 to 1e-2.  Two
compounding root causes :

1. **Right vs left Volterra divide.**  The Hervé–Zaoui closed-form
   interface transition `T = M_b^{-1} · M_a` requires the Volterra
   inverse on the **left** of the numerator.  Our `volterra_divide`
   implemented `num · S^{-vol}` (right) ; the two products are equal
   only when `[num, S] = 0`.  Lower-triangular Volterra trapezoidal
   matrices commute pairwise iff they are Toeplitz (uniform time
   grid + same kernel structure) — non-uniform grids broke this
   assumption silently.

2. **Generic 4n×4n inversion of the shear M-matrix is FP-unstable
   for soft phases.**  Even the block forward-substitution
   `volterra_inverse(_; block_size = 4)` collapses when
   `det(M[t,t]) → 0` (pore-like or step-activated layers).
   ECHOES C++ uses a **closed-form analytic 4×4 inverse** whose only
   `n × n` Volterra inverses are `(3κ + 4μ)^{-vol}` and `μ^{-vol}` —
   both regular for any non-vacuum modulus.

**Fix** :

- Added `volterra_left_divide(S, M; block_size = 1|6)` (forward
  substitution on `S · T = M`, rows i = j..n) and switched every
  closed-form transition (perfect / spring / membrane, both bulk and
  shear) to use it.
- Added `_shear_M_inverse_alv(r, M_κ, M_μ, n)` returning the
  closed-form `M(r; κ, μ)^{-1}` in time-major 4n×4n layout, mirroring
  C++ `inclusion_sphere_nlayers.h::set_visco_inv_matrix_dev` and then
  conjugated to the Christensen–Lo convention.  Used by
  `_shear_layer_transfer_alv` (intra-layer transfer) and
  `_shear_amp_blocks_alv` (state → amplitude extraction).
- Reverted the v0.5.2 τ-scaling for the shear M-matrix — the
  closed-form inverse is naturally written in σ-form and FP
  stability now comes from the closed form rather than the rescaling.

**Impact** : `script 37 :layers` now produces smooth, monotonic
creep curves matching the Python reference figure visually
(bounded between elastic limit and matrix Maxwell), in place of
the previous unbounded / oscillating output.  Bench results :

- `bench_layered_alv.jl` (N=2 stiff elastic + Maxwell matrix,
  uniform grid) : 1e-16 (unchanged).
- `bench_layered_alv_step.jl` (N=3 step-activated layers + pore +
  Maxwell matrix, non-uniform grid) : bulk α 1e-15, shear β 1e-14.
- `bench_step_n2.jl` (N=2 step layers, no pore) : 1e-14.
- `bench_layered_alv_nopore.jl` (N=4 elastic, varied moduli) : 1e-15.

## v0.5.1 — Multi-layer sphere shear localisation bug fix

**Bug fix** : on a non-uniform time grid (e.g. `logspace`) with a
matrix relaxation kernel that has multiple time constants and/or
layers with extreme modulus contrast (pores, step-activated
`ViscoLaw`s), the layered-sphere ALV recurrence diverged from the
ECHOES Python reference by 1e-3 to 1e-2.  Two compounding root
causes :

1. **Right vs left Volterra divide.**  The Hervé–Zaoui closed-form
   interface transition `T = M_b^{-1} · M_a` requires the Volterra
   inverse on the **left** of the numerator.  Our `volterra_divide`
   implemented `num · S^{-vol}` (right) ; the two products are equal
   only when `[num, S] = 0`.  Lower-triangular Volterra trapezoidal
   matrices commute pairwise iff they are Toeplitz (uniform time
   grid + same kernel structure) — non-uniform grids broke this
   assumption silently.

2. **Generic 4n×4n inversion of the shear M-matrix is FP-unstable
   for soft phases.**  Even the block forward-substitution
   `volterra_inverse(_; block_size = 4)` collapses when
   `det(M[t,t]) → 0` (pore-like or step-activated layers).
   ECHOES C++ uses a **closed-form analytic 4×4 inverse** whose only
   `n × n` Volterra inverses are `(3κ + 4μ)^{-vol}` and `μ^{-vol}` —
   both regular for any non-vacuum modulus.

**Fix** :

- Added `volterra_left_divide(S, M; block_size = 1|6)` (forward
  substitution on `S · T = M`, rows i = j..n) and switched every
  closed-form transition (perfect / spring / membrane, both bulk and
  shear) to use it.
- Added `_shear_M_inverse_alv(r, M_κ, M_μ, n)` returning the
  closed-form `M(r; κ, μ)^{-1}` in time-major 4n×4n layout, mirroring
  C++ `inclusion_sphere_nlayers.h::set_visco_inv_matrix_dev`.  Used
  by `_shear_layer_transfer_alv` (intra-layer transfer) and
  `_shear_amp_blocks_alv` (state → amplitude extraction).
- Reverted the v0.5.2 τ-scaling for the shear M-matrix — the
  closed-form inverse is naturally written in σ-form and FP
  stability now comes from the closed form rather than the rescaling.
- Rewrote `_shear_M_matrix_alv` to use the C++ ECHOES mode
  normalisation (mode 1 contributes `U = a · r`, not Christensen–Lo's
  `2a · r`) so it stays consistent with the analytic `M^{-1}`
  formula.  Mode-2 dev contribution factor
  `F_k = (21/5) μ^{-vol} (3κ + μ) (r_b⁵ − r_a⁵)/(r_b³ − r_a³)` was
  unchanged ; it cancels the mode-2 amplitude scaling implicitly.

**Impact** : `script 37 :layers` now produces smooth, monotonic
creep curves matching the Python reference figure visually
(bounded between elastic limit and matrix Maxwell), in place of
the previous unbounded / oscillating output.  Bench results :

- `bench_layered_alv.jl` (N=2 stiff elastic + Maxwell matrix,
  uniform grid) : 1e-16 (unchanged).
- `bench_layered_alv_step.jl` (N=3 step-activated layers + pore +
  Maxwell matrix, non-uniform grid) : bulk α 1e-15, shear β 1e-14.
- `bench_step_n2.jl` (N=2 step layers, no pore) : 1e-14.
- `bench_layered_alv_nopore.jl` (N=4 elastic, varied moduli) : 1e-15.

## v0.5.1 — Multi-layer sphere shear localisation bug fix

**Bug fix** : `LayeredSpheres._shear_localization` and the
companion ALV `shear_localization_alv` previously returned only the
mode-1 amplitude `a_k`, which is the correct per-layer dev β only
when `b_k` (mode-2 amplitude) vanishes — true for `N = 1` (single
sphere) and for the degenerate `N = 2` cases tested in
`test_christensen.jl` (shell ≡ matrix or core ≡ shell).  For
genuinely multi-layer composite spheres with distinct core, shell
and matrix moduli, `b_k` is non-zero and contributes to the
volume-averaged deviatoric strain via the mode-2 r³ profile.

The corrected per-layer dev localisation is

```text
β_k = a_k + b_k · (21/5) (3κ_k + μ_k)/μ_k · (r_k⁵ − r_{k-1}⁵)/(r_k³ − r_{k-1}³)
```

(modes 3 and 4 contribute zero to the layer-volume-averaged
deviatoric strain by angular orthogonality).  This matches ECHOES
C++ `inclusion_sphere_nlayers.h::get_visco_layer_average_strain_Strain`
to machine precision.

**Impact** : the ALV `:layers` topology of `script 37` now matches
the Python `fluage_echoes_solid.py` reference for N ≥ 2.  The
elastic `stiffness_contribution(LayeredSphere, C₀)` and the ALV
`stiffness_contribution_alv(LayeredSphere, C₀_law, times)` produce
the correct effective dilute / MT moduli.

**Validation** : new cross-check benchmark
`scripts/bench_echoes/bench_layered_alv.{py,jl}` and a regression
test `shear_localization_alv — N=2 cross-check vs ECHOES Python` in
`test/Viscoelasticity/test_layered_alv.jl` pin the ALV per-layer
α(t,t') and β(t,t') Volterra blocks to ECHOES Python at machine
precision (1e-16 on the diagonal, 1e-6 on the off-diagonal blocks).

## v0.5.0 — Ageing linear viscoelasticity (ALV) module

A new `MeanFieldHom.Viscoelasticity` sub-module brings time-domain
viscoelastic homogenisation to the package, mirroring the capabilities
of the C++ ECHOES `viscoelasticity/visco_law.h` and
`homogenization_maxwell.h`.  Reference: Sanahuja IJSS 2013 ;
Barthélémy-Giraud-Lavergne-Sanahuja IJSS 2016 ;
Barthélémy-Giraud-Sanahuja-Sevostianov IJES 2019 ; ECHOES manual
chapter 7 and appendix `viscoelastic_hill_kernel.qmd`.

### Highlights

- **`ViscoLaw`** : abstract relaxation `R(t,t')` or creep `J(t,t')`
  kernel, scalar- or 4-tensor-valued, with built-in convenience
  constructors `maxwell_relaxation`, `kelvin_creep`, `maxwell_iso`,
  `kelvin_iso`, `heaviside_law`.
- **`trapezoidal_matrix`** : Sanahuja-2013 trapezoidal discretisation of
  the Stieltjes integral on a time grid `times`, returning a dense
  `Matrix{T}` of size `(B·n) × (B·n)` in lower-block-triangular form
  (`B = 6` for 4-tensor in Mandel convention, `B = 1` for scalar
  kernels).
- **`volterra_inverse`** : block-triangular forward-substitution that
  takes a discrete relaxation matrix to its discrete creep matrix in
  `O(B³ n²)` flops.
- **`hill_kernel`** : discrete ALV Hill polarisation tensor for an
  ellipsoidal inclusion in an isotropic ALV matrix, using the
  time-space decoupling formula of the manual appendix : reuses the
  elastic auxiliary tensors `tens_UA`, `tens_VA` and combines them with
  two scalar Volterra inverses (longitudinal and shear moduli).
  Machine-precision agreement with the elastic Hill tensor in the
  Heaviside (elastic) limit.
- **`homogenize_alv(rve, scheme, prop; times)`** : public entry point
  that builds the discrete operators for every phase, computes the
  ALV Hill kernel, and dispatches to the appropriate scheme function.
  Implemented schemes : `Voigt`, `Reuss`, `Dilute`, `DiluteDual`,
  `MoriTanaka`, `Maxwell`.  Each one's output coincides with the
  corresponding elastic homogenisation in the Heaviside limit (verified
  to machine precision in the test suite).
- **`Phase.properties` relaxed to `Dict{Symbol, Any}`** : a phase can now
  carry either an elastic `AbstractTens` or a `ViscoLaw` under the
  same key (`:C`).  No regression in the elastic test suite (3421/3421).
- **Scripts** : `scripts/33_visco_law_basics.jl` (kernels + plot),
  `scripts/37_fluage_echoes_solid.jl` (Sanahuja-style ageing creep of a
  solidifying composite, whole-pores topology, mirroring
  `tests/python/creep/fluage_echoes_solid.py` after [@sanahuja2013] and
  chapter 9 §"Ageing creep of solidifying cementitious materials").
- **Tests** : `test/Viscoelasticity/` adds 525 new tests across
  `test_visco_law.jl`, `test_trapezoidal.jl`, `test_volterra_inverse.jl`,
  `test_hill_alv_iso.jl`, `test_schemes_alv.jl`.  Total package test
  count : 3946/3946 PASS.

### Self-Consistent ALV (added in 0.5.0)

- **`self_consistent_alv(rve, prop; times, abstol, reltol, maxiters,
  damping, verbose, select_best)`** — symmetric SC fixed-point iteration
  on the `(6n × 6n)` block matrix.  Each iteration recomputes the
  per-phase Hill kernels using the running estimate's iso parameters,
  computes the dilute concentration tensors, and forms the next
  iterate.  Convergence on the Frobenius norm.
- Plumbed into the dispatcher via `homogenize_alv(rve, SelfConsistent(),
  :C; times = T)`.
- Tests against the elastic SC limit pass at machine precision.

### N-layer sphere ALV — full bulk + shear recurrence (added in 0.5.0)

- **`bulk_localization_alv(sphere::LayeredSphere, C0_law, times)`** —
  per-layer bulk localisation matrices `α_k(t,t')` of size `n × n`
  (one per layer).  Extends the elastic Hervé-Zaoui bulk recurrence
  ([`LayeredSpheres/bulk_recurrence.jl`]) by replacing every scalar
  modulus with its `n × n` trapezoidal Volterra matrix, building
  `(2n × 2n)` block transfer matrices.
- **`shear_localization_alv(sphere::LayeredSphere, C0_law, times)`** —
  per-layer deviatoric (Y₂-harmonic) localisation matrices `β_k(t,t')`
  of size `n × n`.  Builds the Hervé-Zaoui 1993 4×4 fundamental matrix
  in **time-major** layout (`(4n × 4n)` block-lower-triangular with
  4×4 diagonal blocks), inverts it via `volterra_inverse(_;
  block_size = 4)`, propagates two probe states through the layers,
  and selects the linear combination matching unit far-field
  `(a_{N+1}, b_{N+1}) = (I_n, 0)` via a final `(2n × 2n)`
  `block_size = 2` Volterra solve.  Verified to machine precision
  against `LayeredSpheres._shear_localization` in the Heaviside
  (elastic) limit.
- **`bulk_state_seq_alv`** / `_shear_state_seq_alv` — forward
  propagation of the discrete state vectors through every layer.
- **ALV interface transfers** (`_bulk_interface_T_alv`,
  `_shear_interface_T_alv`) cover the same set of imperfect interface
  models as the elastic counterpart : `PerfectInterface`,
  `SpringInterface(kn, kt)` (primal — displacement jump driven by
  `kn`/`kt`) and `MembraneInterface(κs, μs)` (dual — traction jump
  from surface elasticity).  Each interface parameter may be a plain
  scalar (constant in time, the elastic limit) **or** a `ViscoLaw`
  scalar kernel — in the latter case the jump itself is ageing and
  the corresponding `n × n` block is the parameter's trapezoidal
  matrix.  The `(4n × 4n)` shear block-diagonal (4×4 sense) cleanly
  reduces to the elastic 4×4 jump for scalar parameters.
- **Composite-sphere assembly**:
  `strain_strain_loc_alv(sphere, C0_law, times)` builds the volume-
  averaged strain-strain localisation tensor
  `A_avg = ⟨α⟩ 𝕁 + ⟨β⟩ 𝕂` (`(6n × 6n)`), and
  `stiffness_contribution_alv(sphere, C0_law, times)` builds the
  size-independent stiffness contribution
  `N = 3 Σ_k f_k (M_κ_k − M_κ_0) ∘ α_k 𝕁
       + 2 Σ_k f_k (M_μ_k − M_μ_0) ∘ β_k 𝕂` (`(6n × 6n)`).
- **`homogenize_alv` extended to `LayeredSphere` phases**: the
  per-inclusion quantities (`A_dut`, `N_dut`) are computed via the
  layered-sphere recurrence instead of the Hill kernel + dilute
  pipeline.  Verified to machine precision against the elastic
  reference for the Dilute scheme and against the elastic MT for the
  `t = t' = 0` block of a Maxwell relaxation kernel.
- **`scripts/37_fluage_echoes_solid.jl`** now exposes a `MODEL`
  constant (`:whole_pores` / `:layers`) selecting the topology;
  the `:layers` branch reproduces the Python `sphere_nlayers(...)`
  setup of `tests/python/creep/fluage_echoes_solid.py` exactly.

### Deferred to follow-up

- Cracks in ALV (extrapolating from the elastic `cod_tensor` /
  `compliance_contribution` infrastructure).
- Self-Consistent ALV with `LayeredSphere` phases (the current SC
  ALV iteration handles only `Ellipsoid`-geometry inclusions).
- Anisotropic ALV Hill kernel (numerical surface integral with Volterra
  inverse of the 3×3 acoustic tensor at each integration point).

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
