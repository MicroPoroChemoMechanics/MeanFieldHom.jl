# Tutorials

Mean-field homogenization replaces a heterogeneous microstructure — a matrix
carrying inclusions, pores, or cracks — by an equivalent homogeneous medium with
the same overall response. `MeanFieldHom` computes that response from phase
properties, geometries, volume fractions, and a **scheme** encoding an
assumption about how the phases interact.

Simplest first within each group. Pages under `generated/` are produced from the
runnable demos in `scripts/` by
[Literate.jl](https://github.com/fredrikekre/Literate.jl).

## Fundamentals

The **porous material** gets two pages: the simplest non-trivial microstructure,
and the one where the choice of scheme matters most.

| Page | What it shows |
| :--- | :--- |
| [A first homogenization](first_estimate.md) | build an RVE; dilute vs Mori–Tanaka |
| [Bounds and classical schemes](bounds_and_schemes.md) | Voigt/Reuss bounds, self-consistent, where each scheme sits |
| [Porous materials and the self-consistent trap](porous_materials.md) | why soft pores break the naive SC iteration |
| [Porous benchmark: all schemes](porous_benchmark.md) | every scheme on the canonical porosity sweep, spheres and oblate pores |
| [The differential scheme and path dependence](differential_paths.md) | incremental homogenization; why mixing order matters |
| [Comparing loading-path trajectories](differential_loading_paths.md) | the same target fractions, four trajectories, watched τ by τ |

## Inclusions, geometries and orientation

| Page | What it shows |
| :--- | :--- |
| [Cracks and crack density](cracks.md) | volume fraction → crack density; the COD tensor |
| [Layered spheres](generated/layered_sphere.md) | Hervé–Zaoui `n`-layer localization and layer averages |
| [Layered spheroids: geometry and effective conductivity](generated/layered_spheroid_effective.md) | the confocal `n`-layer spheroid, the equivalent particle, harmonic-series accuracy |
| [Imperfect interfaces: what they do to the local fields](generated/layered_spheroid_interfaces.md) | pointwise temperature and flux, streamlines, conductance sweep, 3-D view |
| [Highly conducting interfaces](generated/layered_spheroid_hc.md) | equivalent conductivity of an HC-coated particle vs aspect ratio |
| [Symmetrization](generated/symmetrization.md) | exact rotation-group average vs best-fit projection |

Composite inclusions carry no Hill tensor at all: they enter the schemes through
their volume-averaged concentration tensors instead.

## Beyond elasticity

| Page | What it shows |
| :--- | :--- |
| [Viscoelastic composites](viscoelasticity.md) | complex moduli in the frequency domain; a first taste of ageing creep |
| [Frequency or time?](generated/freq_vs_time.md) | the complex-modulus and time-domain ALV routes, cross-checked on the same non-ageing composite |

## Differentiation and solvers

| Page | What it shows |
| :--- | :--- |
| [Derivatives and sensitivities](sensitivities.md) | differentiate any result with `ForwardDiff`, no finite differences |
| [From derivatives to a strength criterion](strength_criteria.md) | those derivatives as a macroscopic strength criterion |
| [Nonlinear solvers for the self-consistent fixed point](nonlinear_solvers.md) | `NonlinearSolve.jl` instead of Picard, and sensitivities that agree either way |
| [Nonlinear homogenization by the secant method](generated/secant_elastoplasticity.md) | elastic–perfectly plastic porous solid, closed by second moments |

## Interoperability and tools

| Page | What it shows |
| :--- | :--- |
| [From Echoes to MeanFieldHom](from_echoes.md) | translate Echoes (C++/Python) scripts; cross-check live via PyCall |
| [Symbolic spheres](symbolic_spheres.md) | the same tensor algebra on `SymPy` / `Symbolics` expressions: Eshelby/Hill tensors and the closed-form estimates |
