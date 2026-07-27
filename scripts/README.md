# `scripts/` — MeanFieldHom.jl demos & validation

Numbered demonstration / validation scripts, grouped in blocks by theme.
Each is self-contained (`Pkg.activate(joinpath(@__DIR__, ".."))`) and, where
relevant, states the reference benchmark it reproduces.

Shared code lives in [`common/`](common/) — currently the Pichler-Hellmich
three-scale model (`common/pichler_model.jl`), used by both the demo script
`41_multiscale_strength.jl` and the cross-check
`bench_echoes/benchmark_pichler.jl`.

## Numbering blocks

| Block | Theme |
|---|---|
| 01–09 | Tensor / Hill / Eshelby toolbox |
| 10–19 | Cracks & COD (16–19 reserved for future conductive / resistive conduction cracks) |
| 20–29 | Elastic homogenization schemes |
| 30–39 | Layered n-layer sphere / spheroid |
| 40–49 | Strength & multiscale (Pichler-Hellmich) |
| 50–59 | Viscoelasticity & ALV |
| 60–69 | ALV cracks, interfaces & cross-validations |
| 70+   | Symmetrization showcases |

## Coverage map

`—` = no direct reference benchmark (native demonstration).

### 01–09 Tensor toolbox
| Script | reference / topic | Notes |
|---|---|---|
| `01_auxiliary_tensors.jl` | — | geometric tensors `tens_IA/UA/VA` |
| `02_hill_elasticity.jl` | `eshelby`/`hill` API | Hill P, elasticity |
| `03_hill_conductivity.jl` | 2nd-order `hill` | conductivity Hill |
| `04_forwarddiff.jl` | — | AD through Hill tensors |
| `05_symbolic.jl` | — | SymPy genericity |
| `06_cylinder.jl` | cylinder Hill | transverse-plane quadrature |
| `07_hill_ti_coaxial.jl` | `hill(...,TI)` | Barthélémy 2020 TI-coaxial closed form |
| `08_hill_derivatives.jl` | `hill_derivative` | ∂P/∂C by ForwardDiff (ISO, TI), validated vs finite differences |

### 10–19 Cracks & COD
| Script | reference / topic | Notes |
|---|---|---|
| `10_cod_isotropic.jl` | `crack_compliance` (iso) | COD / H tensor |
| `11_cod_TI.jl` | `crack_compliance` (TI) | Hoenig / Kanaun-Levin |
| `12_cod_aniso_residue.jl` | `crack_compliance(...,RESIDUES)` | general anisotropy |
| `13_cod_ribbon.jl` | ribbon crack | 2D ribbon COD |
| `14_sif_computation.jl` | — | stress/displacement intensity factors |
| `15_cracks_iso_interface.jl` | iso cracks + spring interface | Sevostianov spring interface |

### 20–29 Elastic schemes
| Script | reference / topic | Notes |
|---|---|---|
| `20_voigt_reuss_bounds.jl` | VOIGT/REUSS | bounds |
| `21_dilute_vs_mori_tanaka.jl` | DIL/MT | dilute vs MT |
| `22_self_consistent_porous.jl` | SC | porous SC percolation |
| `23_differential_trajectories.jl` | DIFF | Norris DEM trajectories |
| `24_differential_loading_paths.jl` | DIFF | path-dependence demo |
| `25_echoes_crosscheck.jl` | Christensen 1990 | cross-check |
| `26_sensitivities.jl` | `homogenize_derivative` | AD sensitivities tour |
| `27_user_inclusion_sensitivity.jl` | — | user-defined inclusion + AD |
| `28_porous_schemes.jl` | porous benchmark | porous scheme comparison |
| `29_symbolic_schemes.jl` | — | SymPy/Symbolics closed forms: Eshelby/Hill, dilute, MT, porous/rigid limits, hand-derived self-consistent |

### 30–39 Layered n-layer sphere / spheroid
| Script | reference / topic | Notes |
|---|---|---|
| `30_average_nlayers.jl` | n-layer sphere | volume-average concentration (sphere) |
| `31_local_nlayers.jl` | n-layer sphere | pointwise localization fields (sphere) |
| `32_spheroid_effective_conductivity.jl` | Kushch 2015 setting | **published tutorial** — confocal geometry and API, Kapitza sweep, exact equivalent particle (`𝐤ᵉᑫ = 𝐁_Ω·𝐀_Ω⁻¹`), series-truncation convergence and quadrature vs. BigFloat |
| `35_spheroid_interfaces.jl` | local fields (spheroid) | **published tutorial** — what an interface does to the local fields: temperature map, streamlines (bilateral seeding), GIF over β, interactive 3D |
| `37_spheroid_hc_conductivity.jl` | Kushch 2015, HC interface | **published tutorial** — highly conducting (surface-conductive) interfaces vs. aspect ratio |

### 40–49 Strength & multiscale
| Script | reference / topic | Notes |
|---|---|---|
| `40_porous_strength_criterion.jl` | — | porous strength criterion |
| `41_multiscale_strength.jl` | Pichler et al. (CCR 2011) | full 3-scale + strength (ω=1e4). Cross-checked in `bench_echoes/benchmark_pichler.jl` (moduli 1 %, fc 2 %) |
| `42_cementpaste_iso.jl` | Pichler et al. (CCR 2011), ISO | elasticity-only ISO variant (**ω=100**, αmax·(1−1e-3)) |
| `43_secant_elastoplasticity.jl` | Suquet (1997) / Ponte Castañeda (1991); Gurson (1977) | **published tutorial** — modified secant method on a porous plastic solid: n-shell composite sphere + SC + `ForwardDiff` second moments; ported from echoes `echoes_tests/elastoplasticity_porous.py` |

### 50–59 Viscoelasticity & ALV
| Script | reference / topic | Notes |
|---|---|---|
| `50_visco_law_basics.jl` | `visco_law` | Maxwell/Kelvin kernels |
| `51_frequency_sweep_viscoelastic.jl` | complex moduli | frequency sweep |
| `52_rabotnov_mittag_leffler.jl` | Rabotnov / Mittag-Leffler | Rabotnov closed form |
| `53_ageing_creep_solid.jl` | solidifying creep | ALV creep |
| `54_ageing_creep_ellipsoid2.jl` | ellipsoid-2 creep | ALV creep |
| `55_ageing_creep_dirichlet_chains.jl` | Granger creep | ageing creep (Granger–Bažant 1995 law) |
| `56_ageing_creep_order2.jl` | order-2 creep | order-2 ALV |
| `57_ageing_creep_cracks.jl` | crack creep | ALV crack creep |
| `58_alv_kernel_types.jl` | — | structured ALV kernel types |
| `59_alv_sensitivities.jl` | — | AD through the ALV pipeline |

### 60+ ALV cracks, cross-validations / symmetrization
| Script | reference / topic | Notes |
|---|---|---|
| `60_alv_cracks_interface.jl` | crack + interface creep | finite interface stiffness |
| `61_freq_vs_time.jl` | Sanahuja (2013) trapezoidal Volterra | **published tutorial** — complex-modulus route vs. `homogenize_alv`, cross-checked through a forward Laplace-Carson transform; O(Δt²) agreement. Ported from echoes `creep/comparison_freq_time.py` |
| `70_symmetrization_showcase.jl` | `symmetrize` / `.paramsym` | **exact rotation average vs best-fit projection** on a non-major-symmetric concentration tensor |

## Conventions worth knowing

- **Exact vs best-fit symmetrization.** Inside scheme kernels the orientation
  average is EXACT (`transverse_isotropify` → `TensTI{4,T,8}`, non-major-
  symmetric content preserved). `best_fit_ti` (→ `TensTI{4,T,5}`) is the
  echoes `.paramsym(sym=TI)` reporting projection — never used in kernels.
  `70_symmetrization_showcase.jl` demonstrates the difference.
- **Water/air TINY = 1e-3.** The Pichler scripts regularize the exactly-zero
  echoes water/air stiffness with a small positive `TINY`, which selects the
  physical (percolating) Self-Consistent branch. Expect a matching small
  offset from echoes near α→0. The ISO variant (`42_cementpaste_iso.jl`) uses
  the exact echoes convention where it is robust.
- **Needle aspect ratio.** The full CCR2011 model uses ω = 1e4; the companion
  iso variant uses ω = 100 (both faithful to their echoes originals).

## Not yet ported
Biaxial strength envelope (Pichler et al., CCR 2013) and the multi-model
`E(w/c)` comparison — future ports.

## Literate.jl convention (pilot, 2026-07-24)

A script converted to this contract stays runnable exactly as before
(`julia scripts/NN_*.jl`) **and** becomes a source for
[Literate.jl](https://github.com/fredrikekre/Literate.jl), which generates a
Documenter markdown page, a Jupyter notebook, and a cleaned standalone
script from the same file (`julia --project=docs docs/literate.jl`).

**Publication policy** — a script is only *published as a tutorial page*
(added to `PUBLISHED_SCRIPTS` in `docs/literate.jl` and to the `pages`
tree of `docs/make.jl`) if no other tutorial or application already
covers its topic. Scripts that duplicate one (the majority — see
`Assets/plans/MFH_LITERATE_SCRIPTS.md` for the full classification) keep
the plain banner style and are never regenerated into a competing page.

There is no separate "Gallery" section any more: a page generated from a
script and one written by hand are both tutorials, and the reader has no
reason to care which is which. `PUBLISHED_SCRIPTS` maps each script to
its **page name**, so scripts keep their numeric prefixes (a running
order) while pages carry thematic names — inserting a tutorial never
forces a renumbering. That mapping is what Literate's `name` option is
for.

Converting a script to the contract, whether or not it ends up promoted:

- **Title & prose.** Replace the `# ===...===` banner with a Literate
  header: `# # Title`, then prose paragraphs as plain `# ` lines, math as
  ```` # ```math ... ``` ````, section dividers as `# ## §N Title`.
- **`Pkg.activate`.** Suffix both the `import Pkg` and `Pkg.activate(...)`
  lines with `#jl` — kept in the standalone script and the generated
  "cleaned script", stripped from the generated markdown/notebook (which
  run inside the `docs` environment, where `MeanFieldHom` is already
  available via `[sources] path=".."`).
- **Figures.** End the plotting code with the plot object as a bare,
  unmarked final expression (captured inline by `@example`/notebook
  execution). Suffix `figdir`/`mkdir`/`savefig`/`display`/the "Saved:"
  `@printf` with `#jl` — the standalone run still writes the PNG to
  `scripts/figures/`, the doc page shows the figure inline instead.
- **Determinism.** Any script using `Random` needs `Random.seed!(<const>)`
  near the top, so the generated doc page (and notebook) render identical
  numbers on every rebuild.
- **Don't combine `#md #nb` on one line** — Literate's marker matching only
  recognizes a single trailing tag; a line with two strips it from *all*
  three outputs. `gr()` needs no marker at all — leave it plain, exactly as
  the hand-written tutorials already do.

See `docs/literate.jl` for the generator entry point and
`Assets/plans/MFH_LITERATE_SCRIPTS.md` for the gap-filler vs. duplicate
classification of all 41 scripts.
