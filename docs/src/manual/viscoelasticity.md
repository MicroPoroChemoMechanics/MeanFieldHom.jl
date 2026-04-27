# Viscoelastic homogenisation — user manual

The ALV pipeline reuses the [`RVE`](@ref) machinery of the elastic
side. Switching to viscoelasticity is essentially a matter of
replacing each phase property by a [`ViscoLaw`](@ref) and passing
`times` to [`homogenize_alv`](@ref).

## 1. Defining a constitutive law

A [`ViscoLaw`](@ref) wraps a 2-argument kernel function and a mode
flag (`:relaxation` for `R(t,t')`, `:creep` for `J(t,t')`). The
return value can be:

* a scalar / complex number (1D ALV problems),
* a `TensISO{4,3}` / `TensTI{4}` / `TensOrtho` / `Matrix{6×6}` (3D
  4-tensor in Mandel form),
* a `TensISO{2,3}` / `Matrix{3×3}` (3D 2-tensor for conductivity).

```julia
using MeanFieldHom, TensND

# Maxwell isotropic relaxation : R(t,t') = (3K∞ + (3K₀-3K∞) exp(-(t-t')/τ_K)) 𝕁
#                                + (2μ∞ + (2μ₀-2μ∞) exp(-(t-t')/τ_μ)) 𝕂
function R_iso(t, tp)
    α = 3 * (1.0 + 4.0 * exp(-(t - tp) / 1.0))
    β = 2 * (0.5 + 1.5 * exp(-(t - tp) / 0.5))
    return TensISO{3}(α, β)
end
law_M = ViscoLaw(R_iso, :relaxation)
```

Pre-built constructors exist for common cases:

```julia
maxwell_iso(k, μ, η_k, η_μ)            # 3K e^{-t/η_k} 𝕁 + 2μ e^{-t/η_μ} 𝕂
kelvin_iso(k_inst, μ_inst, k_inf, μ_inf, η_k, η_μ)
heaviside_law(C)                       # elastic limit : R(t,t') = C·H(t-t')
```

## 2. Trapezoidal discretisation

Given a time grid `times = [t_1, …, t_n]`, the Stieltjes integral
becomes a `(B·n × B·n)` block-lower-triangular matrix:

```julia
times = collect(range(0.0, 5.0; length = 50))
M = trapezoidal_matrix(law_M, times)         # 300 × 300 (= 6 · 50)
```

[`volterra_inverse`](@ref) takes the relaxation discretisation to the
creep discretisation (and vice versa) via block forward substitution
on the diagonal blocks (size `B`). A `LowerTriangular` BLAS path is
selected internally for `B ≥ 2`.

## 3. Building an RVE and homogenising

```julia
rve = RVE(:M)
add_matrix!(rve, Ellipsoid(1.0), Dict(:C => law_M))
add_phase!(rve, :I, Ellipsoid(1.0, 1.0, 0.5),
           Dict(:C => heaviside_law(TensISO{3}(60.0, 20.0)));
           fraction = 0.2)

C_eff = homogenize_alv(rve, MoriTanaka(), :C; times = times)
```

[`homogenize_alv`](@ref) returns a `(6n × 6n)` block matrix. Schemes
supported : `Voigt`, `Reuss`, `Dilute`, `DiluteDual`, `MoriTanaka`,
`Maxwell`, `PonteCastanedaWillis`, `SelfConsistent`,
`AsymmetricSelfConsistent`, `DifferentialScheme`.

A symmetric companion exists for **conductivity / diffusion**
(2-tensor 3 × 3) — the dispatcher inspects the sample type returned
by `visco_eval(law, t, t)` and routes to the order-2 pipeline
automatically.

## 4. Cracks in ALV

```julia
add_phase!(rve, :C, PennyCrack(1.0), Dict(:C => law_M);
           density = 0.05, symmetrize = :iso)
```

`PennyCrack`, `EllipticCrack` and `RibbonCrack` are accepted (penny
limit: traction-free; interface-stiffness extension is on the
roadmap). The dispatcher pre-aggregates crack stiffness / compliance
contributions and routes them through the appropriate scheme branch.

## 5. Symmetry-class fast paths

When all phases share an iso / TI / ortho symmetry with compatible
axes, [`homogenize_alv`](@ref) automatically routes through a fast
path that solves the scheme algebra in the **structured** domain :

| Path  | Components       | Volterra inverse cost                |
|-------|------------------|--------------------------------------|
| ISO   | (α, β)           | 2 × scalar `n × n` forward solves    |
| TI    | (ℓ₁, …, ℓ₆)      | (2n × 2n) block + 2 scalar solves    |
| ORTHO | (o₁, …, o₁₂)     | (3n × 3n) block + 3 scalar solves    |

Detection is heuristic (`_is_iso_block` / `_is_ti_block` /
`_is_ortho_block`) — the user never asks for a fast path explicitly.
The output is still returned as a dense `(6n × 6n)` `Matrix{T}` for
backwards compatibility.

## 6. Structured ALV kernel types

For user code that wants to keep the compact storage **and** the
type information, the structured wrappers
[`ALVKernelISO`](@ref) / [`ALVKernelTI`](@ref) /
[`ALVKernelOrtho`](@ref) are `AbstractMatrix{T}` subtypes:

```julia
M = trapezoidal_matrix(law_M, times)
K_iso = ALVKernelISO(M)            # extracts (α, β), 18× cheaper storage

# Algebra closure stays in the structured class, no (6n×6n) materialisation
K_prod = K_iso * K_iso             # ALVKernelISO
K_inv  = volterra_inverse(K_iso)   # ALVKernelISO

# Auto-promotion iso ⊂ TI ⊂ ortho when mixing
K_TI = ALVKernelTI(K_iso)
K_O  = ALVKernelOrtho(K_iso)
K_iso + K_TI                       # ALVKernelTI
K_iso * K_O                        # ALVKernelOrtho

Matrix(K_iso)                      # back to dense (6n × 6n) on demand
```

These types are a **prototype**: they are fully usable for hand-rolled
ALV pipelines but `homogenize_alv` does not yet accept them as inputs
(use `Matrix(K)` to cross the boundary).

## 7. Reproducing ECHOES C++ benchmarks

`scripts/37_fluage_echoes_solid.jl` reproduces the multi-phase Maxwell
+ solidifying Maxwell + pore benchmark from the ECHOES manual.
`scripts/41_fluage_echoes_cracks.jl` covers all seven crack-aware ALV
schemes on a penny-crack RVE. The cross-check
`scripts/bench_echoes/benchmark_alv.jl` validates `homogenize_alv` vs
the C++ implementation via PyCall (rtol ≤ 1e-8 on the (1,1) Mandel
block, ≤ 1e-6 on the full matrix).
