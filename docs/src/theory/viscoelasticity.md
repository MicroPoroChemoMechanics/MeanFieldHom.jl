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

| Scheme                           | Crack treatment                                             |
|----------------------------------|-------------------------------------------------------------|
| Voigt / Reuss                    | ignored (zero-volume convention)                            |
| Dilute / DiluteDual              | additive `+ ΔC̃_cracks`                                     |
| Mori-Tanaka / Maxwell / PCW      | virtual phase with `Ã = 0`, `Ñ = ΔC̃_crack`, `f = 1`        |
| SC / ASC                         | re-evaluated against the running effective estimate         |

## References

[@sanahuja2013] [@barthelemyIJSS2016] [@barthelemyIJES2019]
