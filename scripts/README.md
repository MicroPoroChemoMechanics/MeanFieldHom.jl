# `scripts/` — MeanFieldHom.jl demos & echoes cross-checks

Numbered demonstration / validation scripts, grouped in blocks by theme.
Each is self-contained (`Pkg.activate(joinpath(@__DIR__, ".."))`) and, where
relevant, states the echoes (C++/Python) counterpart it mirrors.

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
| 30–39 | Layered n-layer sphere |
| 40–49 | Strength & multiscale (Pichler-Hellmich) |
| 50–59 | Viscoelasticity & ALV |
| 60–69 | ALV cracks / interfaces |
| 70+   | Symmetrization showcases |

## Coverage map (script ↔ echoes)

Paths under `echoes_cpp/tests/python/` unless noted; `—` = no direct echoes
counterpart (native demonstration).

### 01–09 Tensor toolbox
| Script | echoes counterpart | Notes |
|---|---|---|
| `01_auxiliary_tensors.jl` | — | geometric tensors `tens_IA/UA/VA` |
| `02_hill_elasticity.jl` | `eshelby`/`hill` API | Hill P, elasticity |
| `03_hill_conductivity.jl` | 2nd-order `hill` | conductivity Hill |
| `04_forwarddiff.jl` | — | AD through Hill tensors |
| `05_symbolic.jl` | — | SymPy genericity |
| `06_cylinder.jl` | cylinder Hill | transverse-plane quadrature |
| `07_hill_ti_coaxial.jl` | `hill(...,TI)` | Barthélémy 2020 TI-coaxial closed form |
| `08_hill_derivatives.jl` | `hill_derivative` / `derive_eshelby.py` | ∂P/∂C by ForwardDiff (ISO, TI), validated vs finite differences |

### 10–19 Cracks & COD
| Script | echoes counterpart | Notes |
|---|---|---|
| `10_cod_isotropic.jl` | `crack_compliance` (iso) | COD / H tensor |
| `11_cod_TI.jl` | `crack_compliance` (TI) | Hoenig / Kanaun-Levin |
| `12_cod_aniso_residue.jl` | `crack_compliance(...,RESIDUES)` | general anisotropy |
| `13_cod_ribbon.jl` | ribbon crack | 2D ribbon COD |
| `14_sif_computation.jl` | — | stress/displacement intensity factors |
| `15_cracks_iso_interface.jl` | `cracksiso.py` | Sevostianov spring interface |

### 20–29 Elastic schemes
| Script | echoes counterpart | Notes |
|---|---|---|
| `20_voigt_reuss_bounds.jl` | VOIGT/REUSS | bounds |
| `21_dilute_vs_mori_tanaka.jl` | DIL/MT | dilute vs MT |
| `22_self_consistent_porous.jl` | SC | porous SC percolation |
| `23_differential_trajectories.jl` | DIFF | Norris DEM trajectories |
| `24_differential_loading_paths.jl` | DIFF | path-dependence demo |
| `25_echoes_crosscheck.jl` | Christensen 1990 | cross-check |
| `26_sensitivities.jl` | `homogenize_derivative` | AD sensitivities tour |
| `27_user_inclusion_sensitivity.jl` | — | user-defined inclusion + AD |
| `28_porous_schemes.jl` | `echoes_tests/porous.py` | porous scheme comparison |
| `29_symbolic_schemes.jl` | — | SymPy/Symbolics closed forms: Eshelby/Hill, dilute, MT, porous/rigid limits, hand-derived self-consistent |

### 30–39 Layered n-layer sphere
| Script | echoes counterpart | Notes |
|---|---|---|
| `30_average_nlayers.jl` | `spheroid_nlayers/` | volume-average concentration |
| `31_local_nlayers.jl` | `spheroid_nlayers/` | pointwise localization fields |

### 40–49 Strength & multiscale
| Script | echoes counterpart | Notes |
|---|---|---|
| `40_porous_strength_criterion.jl` | — | porous strength criterion |
| `41_multiscale_strength.jl` | `cementpaste_mortar_Pichler_CCR2011.py` | full 3-scale + strength (ω=1e4). Cross-checked in `bench_echoes/benchmark_pichler.jl` (moduli 1 %, fc 2 %) |
| `42_cementpaste_iso.jl` | `cementpaste_mortar_iso_Pichler_CCR2011.py` | elasticity-only ISO variant (**ω=100**, αmax·(1−1e-3)) |

### 50–59 Viscoelasticity & ALV
| Script | echoes counterpart | Notes |
|---|---|---|
| `50_visco_law_basics.jl` | `visco_law` | Maxwell/Kelvin kernels |
| `51_frequency_sweep_viscoelastic.jl` | complex moduli | frequency sweep |
| `52_rabotnov_mittag_leffler.jl` | `mittag_leffler/` | Rabotnov closed form |
| `53_ageing_creep_solid.jl` | `creep/` solid layers | ALV creep |
| `54_ageing_creep_ellipsoid2.jl` | `creep/` ellipsoid2 | ALV creep |
| `55_ageing_creep_dirichlet_chains.jl` | `creep/` Granger | ageing creep (Granger–Bažant 1995 law) |
| `56_ageing_creep_order2.jl` | `creep/` order-2 | order-2 ALV |
| `57_ageing_creep_cracks.jl` | `creep/` cracks | ALV crack creep |
| `58_alv_kernel_types.jl` | — | structured ALV kernel types |
| `59_alv_sensitivities.jl` | — | AD through the ALV pipeline |

### 60+ ALV cracks / symmetrization
| Script | echoes counterpart | Notes |
|---|---|---|
| `60_alv_cracks_interface.jl` | `creep/` cracks + interface | finite interface stiffness |
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

## Not yet ported (echoes side)
`cementpaste_mortar_Pichler_biax_CCR2013.py` (biaxial strength envelope) and
`cementpaste.py` (multi-model `E(w/c)` comparison) — future ports.
