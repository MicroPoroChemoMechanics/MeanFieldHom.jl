# Ageing linear viscoelasticity (ALV)

`MeanFieldHom.Viscoelasticity` extends the elastic homogenisation
framework to **ageing** linear viscoelastic media [@sanahuja2013;
@barthelemyIJSS2016; @barthelemyIJES2019]. The constitutive law is a
two-time kernel — relaxation `R(t,t')` or compliance / creep `J(t,t')`
— linking stress and strain through a Stieltjes integral

```math
\boldsymbol{\sigma}(t) = \int_{-\infty}^{t} \mathbb{R}(t,t')\,
                          \mathrm{d}\boldsymbol{\varepsilon}(t').
```

After discretisation on a time grid `t = (t_1, …, t_n)` via the
trapezoidal rule, every viscoelastic operator (kernel, Hill tensor,
localisation tensor, effective stiffness) becomes a **lower
block-triangular** matrix of size `(B·n) × (B·n)` with `B = 6` for
4-tensors (Mandel form) and `B = 3` for 2-tensors (conductivity /
diffusion). Causality is encoded in the lower-triangular structure;
products and Volterra inverses act on the time variable as ordinary
matrix products / forward-substitutions.

## Time-space decoupling for an isotropic matrix

For an **isotropic** ALV matrix the Hill polarisation tensor admits a
particularly clean factorisation [@barthelemyIJSS2016, App.] :

```math
\widetilde{\mathbb{P}}_{\mathcal{E}}(t,t')
   = \widetilde{(k+\tfrac{4}{3}\mu)}^{-\mathrm{vol}}_{(t,t')}\,
       \mathbb{U}^{\boldsymbol{A}}
   + \widetilde{\mu}^{-\mathrm{vol}}_{(t,t')}\,
       \bigl(\mathbb{V}^{\boldsymbol{A}} - \mathbb{U}^{\boldsymbol{A}}\bigr).
```

The geometric tensors `U^A`, `V^A` are **identical to the elastic
case** (already computed in
[`MeanFieldHom.Elasticity`](@ref)); only two scalar Volterra
inverses carry the time dependence. The full discrete Hill matrix
`P̃ ∈ ℝ^{6n×6n}` is then assembled by adding `n²` Mandel blocks of
size `6 × 6`.

## Symmetry classes and structured storage

ALV operators inherit the same symmetry classes as their elastic
counterparts. Closure under Volterra product / inverse defines the
following compact storage hierarchy:

| Class    | Stored components | (6n × 6n) entries | Storage cost   | Closure operation                       |
|----------|-------------------|-------------------|----------------|-----------------------------------------|
| ISO      | (α, β)            | 36 n²             | **2 n²** (18×) | scalar Volterra products / inverses     |
| TI       | (ℓ₁, …, ℓ₆)       | 36 n²             | **6 n²** (6×)  | (2n × 2n) block-Volterra + 2 scalars    |
| ORTHO    | (o₁, …, o₁₂)      | 36 n²             | **12 n²** (3×) | (3n × 3n) block-Volterra + 3 scalars    |
| Generic  | full (6n × 6n)    | 36 n²             | 36 n²          | (6n × 6n) block-LU                      |

Iso ⊂ TI ⊂ ortho ⊂ generic. The structured types
[`ALVKernelISO`](@ref), [`ALVKernelTI`](@ref) and
[`ALVKernelOrtho`](@ref) wrap these compact representations as
`AbstractMatrix{T}` so they flow through Julia generic matrix code
while preserving the storage savings and algebra closure.

## Cracks in the ALV pipeline

A flat crack (penny / elliptic / ribbon) carries no volume but
contributes a `ΔC̃_crack` term to the **numerator** of every scheme
[@barthelemyIJES2019, §4]. The ALV crack-opening-displacement (COD)
tensor `B̃` factorises like the iso Hill kernel : two scalar Volterra
matrices (one normal, one tangential compliance) times geometric
4-tensors. From `B̃` follow the compliance contribution `H̃` and the
stiffness contribution `Ñ = -C̃_M·H̃·C̃_M`.

Schemes integrate cracks as follows:

| Scheme                      | Crack treatment                                       |
|-----------------------------|-------------------------------------------------------|
| Voigt / Reuss               | ignored (zero-volume convention)                      |
| Dilute / DiluteDual         | additive `+ delta_C_cracks`                           |
| Mori-Tanaka / Maxwell / PCW | virtual phase with `A = 0`, `N = delta_C`, `f = 1`    |
| SC / ASC                    | re-evaluated against the running effective estimate   |

## Numerical example — minimal working ALV pipeline

```julia
using MeanFieldHom, TensND

# Iso Maxwell relaxation kernel
function R_iso(t, tp)
    α = 3 * (3.0 + 2.0 * exp(-(t - tp) / 1.0))
    β = 2 * (1.0 + 1.0 * exp(-(t - tp) / 0.5))
    return TensISO{3}(α, β)
end
law_M = ViscoLaw(R_iso, :relaxation)

# RVE = matrix + iso elastic spheres at f = 0.20
rve = RVE(:M)
add_matrix!(rve, Ellipsoid(1.0), Dict(:C => law_M))
add_phase!(rve, :I, Ellipsoid(1.0), Dict(:C => heaviside_law(TensISO{3}(60.0, 20.0)));
            fraction = 0.20)

# Time grid + Mori-Tanaka homogenisation
times = collect(range(0.0, 5.0; length = 50))
C_eff = homogenize_alv(rve, MoriTanaka(), :C; times = times)

# Read off the effective shear modulus history (column 1 = t' = 0)
α_eff, β_eff = iso_params_from_blocks(C_eff)
μ_eff_history = β_eff[:, 1] ./ 2
```

A worked end-to-end example (multi-phase Maxwell + solidifying gel +
pore) lives in `scripts/37_fluage_echoes_solid.jl`. Closed-form
validation against the Rabotnov / Mittag-Leffler benchmark of
[@barthelemyIJES2019, §5] is in `scripts/36_rabotnov_mittag_leffler.jl`.

## Symmetry-class fast paths — runnable comparison

```julia
# All phases iso → automatic iso fast path (2·n² storage instead of 36·n²)
using MeanFieldHom

α_eff, β_eff = iso_params_from_blocks(C_eff)         # extracted from the result
@assert MeanFieldHom.Viscoelasticity._is_iso_block(C_eff)

# To opt into the structured wrapper (compact storage + algebra closure):
K_iso = ALVKernelISO(C_eff)            # 2·n² entries instead of 36·n²
K_inv = volterra_inverse(K_iso)        # stays ALVKernelISO
Matrix(K_inv) ≈ volterra_inverse(C_eff; block_size = 6)   # cross-check
```

`scripts/42_alv_kernel_types.jl` walks through `ALVKernelISO` /
`ALVKernelTI` / `ALVKernelOrtho` with the algebra ladder iso ⊂ TI ⊂
ortho.

## References

[@sanahuja2013] [@barthelemyIJSS2016] [@barthelemyIJES2019]
