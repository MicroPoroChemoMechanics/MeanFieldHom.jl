# Cross-validation against Echoes

`MeanFieldHom` is a port of the C++ [`Echoes`](@cite echoes) code, so almost
every quantity it computes has an independent reference implementation. This
page collects what is actually compared, at which tolerance, and — just as
important — **where the two codes deliberately disagree**.

The cross-checks live in `scripts/bench_echoes/`. They call Echoes live through
`PyCall`; see [From Echoes to MeanFieldHom](../tutorials/from_echoes.md) for
the API translation table and the `PyCall` setup recipe. A subset of the
viscoelastic benchmarks also ships committed `*_python.json` dumps, so the
Julia side can be checked **without** a working Echoes install.

## The five production cross-checks

| Script | Echoes reference | Validated quantities | Tolerance |
| :--- | :--- | :--- | :--- |
| `benchmark.jl` | Hill ``\mathbb P`` / crack ``\Delta\mathbb S`` (residues, DECUHR) | ``\mathbb P``, ``\Delta\mathbb S`` on 5 + 3 cases | ~1e-6 (iso), machine (aniso) |
| `benchmark_nlayers.jl` | n-layer sphere reference | n-layer averages / local fields | 1e-6 |
| `benchmark_porous.jl` | porous benchmark | porous SC / MT moduli | 1e-6 |
| `benchmark_pichler.jl` | Pichler et al. mortar model [pichler2011](@cite) | ``k, \mu, E`` (mortar) + strength ``f_c`` | moduli 1 %, ``f_c`` 2 % |
| `benchmark_hill_derivative.jl` | Hill-derivative reference | ``\partial\mathbb P/\partial\mathbb C`` | ISO ~1e-15, ORTHO ~1e-6 |

Status of the last full run (recorded in `scripts/bench/DIAGNOSTIC.md`):
`benchmark_pichler` 24/24, `benchmark_hill_derivative` 17/17,
`benchmark_nlayers` local stresses to 5.8e-16, `benchmark_porous` 134/140
(the six known `DifferentialScheme` cases below).

## Agreement measured, by quantity

| Quantity | Agreement | Where |
| :--- | :--- | :--- |
| Hill ``\mathbb P``, anisotropic matrix | machine precision (same residue algorithm both sides) | `benchmark.jl` §1 |
| Hill ``\mathbb P``, isotropic matrix | 1e-6 … 1e-8 (analytic path vs quadrature) | `benchmark.jl` §1 |
| Hill ``\mathbb P``, 2nd order (conductivity) | 1e-6 | `benchmark.jl` §3 |
| ``\partial\mathbb P/\partial\mathbb C`` | ISO ~1e-15, ORTHO ~1e-6 | `benchmark_hill_derivative.jl` |
| Crack ``\mathbb H``, penny | machine precision | `benchmark.jl` §2 |
| Porous moduli — Voigt, Reuss, Dilute, DiluteDual, Maxwell, PCW, MT | ``\le`` 2.5e-14 | `benchmark_porous.jl` |
| Porous moduli — MT, published sweep | ``\le`` 6.4e-8 over ``\varphi \in [0.1, 0.9]`` | [From Echoes](../tutorials/from_echoes.md) |
| Porous moduli — SC | 6.4e-8 at ``\varphi = 0.1``, 5.7e-7 at ``\varphi = 0.3``; both codes collapse toward the numerical floor past ``\varphi \approx 0.4`` | [From Echoes](../tutorials/from_echoes.md) |
| n-layer sphere, bulk ``\alpha_k`` | 4.6e-13 over 30 random 2–8-layer stacks | `benchmark_nlayers.jl` §1 |
| n-layer sphere, shear ``\beta_k`` | 4.4e-14 over the same 30 stacks | `benchmark_nlayers.jl` §1 |
| n-layer sphere, local ``\sigma_{rr}``, ``\sigma_{\theta\theta}`` | 5.8e-16 over 80 points | `benchmark_nlayers.jl` §4 |
| Rotational averages (TI / ISO), C++ transcribed verbatim | < 1e-12 over 20 random draws | `test/Core/test_rotational_average.jl` |
| ALV per-layer ``\alpha(t,t')``, ``\beta(t,t')`` | 1e-16 diagonal, 1e-6 off-diagonal | `test/Viscoelasticity/test_layered_alv.jl` |
| ALV layered benchmarks (4 configurations) | 1e-16 / 1e-15 / 1e-14 / 1e-15 | `scripts/bench_echoes/bench_layered_alv*.jl` |
| Dual-interface (`DUALDISC`) concentrations | 1e-8 | `test/LayeredSpheres/test_interfaces.jl` |
| Dual-interface effective diffusivity | 1e-5 | `test/LayeredSpheres/test_conductivity.jl` |
| Cracked ALV creep, PCW | ``\le`` 1e-4 | `scripts/60_alv_cracks_interface.jl` |
| Cracked ALV creep, SC + interface stiffness | 1.4e-4 | see CHANGELOG v0.5.x |
| Three-scale mortar ``k, \mu, E`` | within 1 % (typically 1e-3) | `benchmark_pichler.jl` |
| Ageing solidifying creep, both topologies | better than 1 % | [Ageing creep](../applications/ageing_creep.md) |

## Deliberate divergences

These are **not** numerical error, and a future contributor should not "fix"
them. Each is a convention or a formulation choice, verified rather than
assumed.

### The elliptic-crack ``\mathbb H`` normalisation

`MeanFieldHom` normalises the crack compliance uniformly by the **minor**
semi-axis ``b``, Echoes by the **major** semi-axis ``a`` for the ellipse:

```math
\mathbb H_{\mathrm{MFH}} = \tfrac{3}{4}\,\hat n \otimes^s \mathbf B \otimes^s \hat n,
\qquad
\mathbb H_{\mathrm{Echoes}} = \tfrac{3\eta}{4}\,\hat n \otimes^s \mathbf B \otimes^s \hat n,
```

so ``\mathbb H_{\mathrm{Echoes}} = \eta\,\mathbb H_{\mathrm{MFH}}`` for the
ellipse and the two agree at ``\eta = 1`` (penny). Measured against Echoes at
``\eta = 0.7, 0.5, 0.3, 0.1``: the ratio is ``\eta`` to four decimals, and on
the `MeanFieldHom` side the ``3/4`` is ``\eta``-independent to machine
precision. Details in [COD tensors](../theory/cod_tensors.md); pinned by
`test/regression/test_crack_cases.jl`.

!!! warning "The penny case hides this"
    At ``\eta = 1`` both conventions coincide. A crack test that only exercises
    the penny limit cannot detect a normalisation error — see
    [Testing conventions](testing_conventions.md).

### Mori-Tanaka on a cracked RVE

`MeanFieldHom` uses the ``\mathbb B \cdot \mathbb A^{-1}`` form, Echoes an
additive closure. The two differ at finite crack density: at ``d = 0.30``,
traction-free, MFH gives ``\varepsilon_{xx}(t\to\infty) \approx 0.481`` against
Echoes' ``0.559``. At low density both MT variants agree with PCW to
``\text{rtol} \le`` 1e-3. Discussion in
[Viscoelasticity manual](../manual/viscoelasticity.md).

### Asymmetric self-consistent, porous oblate

The C++ reference's ASC is a compliance-form Mori-Tanaka with an iterating
reference, which converges to a **different branch** for porous oblate systems
away from percolation. The divergence is structural, not a tolerance issue.

### DifferentialScheme, porous, ``\varphi \ge 0.5``

Six of the 140 `benchmark_porous.jl` cases fall outside tolerance, all
`DifferentialScheme`, with the relative error growing from 2.6e-3 at
``\varphi = 0.50`` to 6.2e-2 at ``\varphi = 0.95``. This is a **pre-existing,
reproducible** gap against Echoes, unchanged digit for digit across the
optimisation campaign — not a regression.

### Strength ``f_c``: the 2 % gap is by construction

Water and air phases are regularised to a small positive stiffness instead of
exact zero. That regularisation is the source of the ~2 % gap on ``f_c``; the
elastic moduli are unaffected at the 1 % level. See
[Strength criteria](../applications/strength.md).

## A divergence that wasn't

Worth recording, because the failure mode is easy to repeat. The n-layer
sphere's shear localization ``\beta_k`` used to be listed as *not* comparable
to Echoes: a direct comparison disagreed by 1–50 % on genuine multi-layer
stacks, the gap was attributed to a layer-indexing convention in
`echoes.layer_eE`, and the check was replaced by §3's analytical limits.

All of that was wrong. ``\beta_k`` was being computed as the bare mode-1
amplitude ``a_k``, dropping the mode-2 contribution
``b_k \cdot F_k`` — the ``r^3`` displacement profile integrates to a non-zero
deviatoric strain across a shell of finite thickness. The omission cancels
exactly in the degenerate configurations that §3 tests (vanishing core,
core ≡ shell, single layer), which is why the implementation could
simultaneously "match the analytical limits" and "disagree with Echoes". The
missing term was added later as part of an unrelated fix, but the conclusion
about Echoes was never re-tested. With the comparison restored, ``\beta_k``
agrees to **4.4e-14** — better than ``\alpha_k``.

Two cheap checks would have caught it immediately:

- ``\alpha_k`` and ``\beta_k`` are read from the **same** `layer_eE(k)`
  matrix, so a layer-index error could not have spared ``\alpha_k``, which
  matched to 5e-13 throughout. The stated explanation was internally
  inconsistent.
- The ALV per-layer ``\beta(t,t')`` Volterra blocks were already pinned to
  Echoes at 1e-16 on the diagonal. Since ALV reduces to elasticity in the
  non-ageing limit, an elastic ``\beta_k`` wrong by tens of percent was
  impossible.

!!! warning "Degenerate test cases hide missing modes"
    A localization term that vanishes in every symmetric or degenerate
    configuration is invisible to a test suite built from those
    configurations. This is the same trap as the penny-crack limit hiding the
    ``\eta`` normalisation above — see
    [Testing conventions](testing_conventions.md).

## Running the cross-checks

The Echoes benchmarks need the weak-dependency extensions loaded, since they
compare back-ends:

```julia
import DECUHR, Integrals      # activates MeanFieldHomDECUHRExt
include("scripts/bench_echoes/benchmark.jl")
```

Without a `PyCall`/Echoes install, the viscoelastic layered benchmarks still
run against their committed JSON dumps
(`scripts/bench_echoes/*_python.json`), and every agreement figure pinned in
the table above is also asserted by the ordinary test suite.
