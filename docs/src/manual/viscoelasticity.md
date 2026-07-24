# Viscoelastic homogenization — user manual

The ALV (ageing linear viscoelastic) pipeline reuses the [`RVE`](@ref)
machinery of the elastic side. Switching to viscoelasticity is
essentially a matter of replacing each phase property by a
[`ViscoLaw`](@ref) and passing a `times` grid to
[`homogenize_alv`](@ref).

This manual walks through eight use cases, each runnable as is.
A more elaborate version of every example exists under
`scripts/50_visco_law_basics.jl` … `scripts/59_alv_sensitivities.jl`.

## 1. Defining a constitutive law

A [`ViscoLaw`](@ref) wraps a two-argument kernel function `(t, t')` and
a mode flag (`:relaxation` for the relaxation kernel `R(t, t')`,
`:creep` for the compliance kernel `J(t, t')`). The kernel can return:

* a scalar (`Float64` / `Complex{Float64}`) — 1D ALV problems,
* a `TensISO{4, 3}` / `TensTI{4}` / `TensOrtho` / `Matrix{6×6}` — 3D
  4-tensor in Mandel form,
* a `TensISO{2, 3}` / `Matrix{3×3}` — 3D 2-tensor for conductivity /
  diffusion / permittivity (order-2 ALV).

### 1.1 Hand-rolled Maxwell isotropic relaxation

```julia
using MeanFieldHom, TensND

# Iso Maxwell relaxation:
#   R(t,t') = (3 K∞ + (3 K₀ - 3 K∞) exp(-(t-t')/τ_K)) 𝕁
#           + (2 μ∞ + (2 μ₀ - 2 μ∞) exp(-(t-t')/τ_μ)) 𝕂
const k₀ = 5.0;  const μ₀ = 2.0
const k∞ = 3.0;  const μ∞ = 1.0
const τ_K = 1.0; const τ_μ = 0.5

function R_iso(t, tp)
    α = 3 * (k∞ + (k₀ - k∞) * exp(-(t - tp) / τ_K))
    β = 2 * (μ∞ + (μ₀ - μ∞) * exp(-(t - tp) / τ_μ))
    return TensISO{3}(α, β)        # iso 4-tensor with parameters (3K, 2μ)
end
law_M = ViscoLaw(R_iso, :relaxation)
```

### 1.2 Pre-built constructors

```julia
# `maxwell_iso(K, μ, τ_K, τ_μ)` —  R = 3K·e^{-t/τ_K} 𝕁 + 2μ·e^{-t/τ_μ} 𝕂
law_max = maxwell_iso(5.0, 2.0, 1.0, 0.5)

# `kelvin_iso(K_∞, μ_∞, K₀, μ₀, τ_K, τ_μ)` — Kelvin (creep) iso
law_kel = kelvin_iso(3.0, 1.0, 5.0, 2.0, 1.0, 0.5)

# Elastic limit : R(t,t') = C · H(t-t')
law_el  = heaviside_law(TensISO{3}(15.0, 4.0))
```

### 1.3 Ageing kernels

The first argument is the current time `t`, the second is the loading
time `t'`. Ageing means the kernel depends on `t'`, not just on the
duration `t − t'` (basic linear viscoelasticity is the special case
where it depends only on `t − t'`):

```julia
# Sanahuja-style solidification : volume fraction of "active" gel grows
# with t' as `f_∞ · t'^α / (1 + t'^α)`.
const α_age = 4.0
const f_∞   = 0.3
@inline solidification(tp) = f_∞ * tp^α_age / (1 + tp^α_age)

function R_aging(t, tp)
    f = solidification(tp)
    α = 3 * (3.0 + (5.0 - 3.0) * f * exp(-(t - tp) / 1.0))
    β = 2 * (1.0 + (2.0 - 1.0) * f * exp(-(t - tp) / 0.5))
    return TensISO{3}(α, β)
end
law_aging = ViscoLaw(R_aging, :relaxation)
```

## 2. Trapezoidal discretization

Given a time grid `times = (t₁, …, tₙ)`, the Stieltjes integral becomes
a `(B·n × B·n)` lower-block-triangular matrix (`B = 6` for a 4-tensor
kernel, `B = 3` for a 2-tensor kernel, `B = 1` for scalar kernels):

```julia
times = collect(range(0.0, 5.0; length = 50))
M = trapezoidal_matrix(law_M, times)      # 300 × 300  (=  6·50)
```

[`volterra_inverse`](@ref) flips a relaxation matrix to the
corresponding creep matrix (and vice-versa) via block forward
substitution. A `LowerTriangular` BLAS path is selected internally for
`B ≥ 2`:

```julia
J = volterra_inverse(M; block_size = 6)            # creep (compliance) matrix
@assert isapprox(M * J,
                 [iszero(rem(i - j, 6)) && (i ÷ 6 == j ÷ 6) for i in 1:300, j in 1:300] |>
                  Matrix{Float64};
                 atol = 1e-10)                       # block-diagonal identity
```

## 3. Building an RVE and homogenising

```julia
# 50-step time grid; the matrix is the Maxwell iso law from §1.1
times = collect(range(0.0, 5.0; length = 50))

# Inclusions : aligned spheroids (oblate ratio 0.5), elastic
C_I = TensISO{3}(60.0, 20.0)

rve = RVE(:M)
add_matrix!(rve, Ellipsoid(1.0), Dict(:C => law_M))
add_phase!(rve, :I, Ellipsoid(1.0, 1.0, 0.5),
            Dict(:C => heaviside_law(C_I));
            fraction = 0.2)

# Homogenise
C_eff = homogenize_alv(rve, MoriTanaka(), :C; times = times)   # 300 × 300
```

[`homogenize_alv`](@ref) accepts: `Voigt`, `Reuss`, `Dilute`,
`DiluteDual`, `MoriTanaka`, `Maxwell`, `PonteCastanedaWillis`,
`SelfConsistent`, `AsymmetricSelfConsistent`, `DifferentialScheme`.

A symmetric companion exists for **conductivity / diffusion /
permittivity** (order-2 properties, 3 × 3 kernels) — the dispatcher
inspects the sample type returned by `visco_eval(law, t, t)` and
routes to the order-2 pipeline automatically (see §6).

## 4. Reading effective properties

The `(6n × 6n)` output is a Mandel block matrix. Use
[`iso_params_from_blocks`](@ref) to extract the iso `(α, β)` parameter
matrices:

```julia
α, β = iso_params_from_blocks(C_eff)         # n × n each
# 3 K_eff(t,t')   = α(t,t')
# 2 μ_eff(t,t')   = β(t,t')

# Effective shear modulus history : column index = t' = 0 (Heaviside step)
times_keep = times                            # (length n)
μ_eff_t = β[:, 1] ./ 2                        # μ(tᵢ, t₁) for i = 1..n
```

For a **uniaxial creep** test (unit longitudinal stress, all
components 1, 2, 3 stay the same in iso ; 4, 5, 6 are zero):

```julia
n = length(times)
J_eff = volterra_inverse(C_eff; block_size = 6)       # creep matrix
S = zeros(n * 6); for i in 1:n;  S[6 * (i - 1) + 1] = 1.0;  end
ε = J_eff * S
ε_xx_t = ε[1:6:end]                            # ε_xx(tᵢ)
```

For TI (axis = e₃), the 4-tensor Walpole parameters
`(ℓ₁, ℓ₂, ℓ₃, ℓ₄, ℓ₅, ℓ₆)` are extracted similarly:

```julia
ℓ = ti_params_from_blocks(C_eff)               # NTuple{6, Matrix}
```

## 5. Cracks in ALV

### 5.1 Traction-free penny crack

```julia
add_phase!(rve, :C, PennyCrack(1.0), Dict(:C => law_M);
            density = 0.05, symmetrize = :iso)
C_eff_cracks = homogenize_alv(rve, MoriTanaka(), :C; times = times)
```

`PennyCrack`, `EllipticCrack` and `RibbonCrack` are accepted. The
dispatcher pre-aggregates crack stiffness and compliance contributions
and routes them through the appropriate scheme branch.

### 5.2 Cracks with finite interface stiffness (Sevostianov)

For a flat crack carrying a **spring-like interface stiffness** with
time-dependent normal `Rn(t,t')` and tangential `Rt(t,t')` ageing
kernels, attach the interface laws as `:Rn` / `:Rt` properties on the
crack phase :

```julia
# Interface kernels — same Maxwell-iso ageing form as the matrix law
R_n_kernel(t, tp) = (1 + 0.1 * tp^0.4) *
                     (1.0e10 + (2.0e10 - 1.0e10) * exp(-(t - tp) / 2.0))
R_t_kernel(t, tp) = (1 + 0.1 * tp^0.2) *
                     (1.0e10 + (1.0e10 - 1.0e10) * exp(-(t - tp) / 3.0))
law_Rn = ViscoLaw(R_n_kernel, :relaxation)
law_Rt = ViscoLaw(R_t_kernel, :relaxation)

add_phase!(rve, :CRACK, PennyCrack(1.0),
            Dict(:C => law_M, :Rn => law_Rn, :Rt => law_Rt);
            density = 0.05, symmetrize = :iso)

C_eff = homogenize_alv(rve, MoriTanaka(), :C; times = times)
```

Behind the scenes the COD matrices `B̃_n`, `B̃_t` are post-corrected via
the Sevostianov spring identity

```math
\widetilde{\mathbf B}_{\text{eff}}
   = (b\,\mathbb K + \widetilde{\mathbf B}^{-1})^{-\text{vol}}
   = \widetilde{\mathbf B} \circ
     (\mathbb 1 + b\,\mathbb K \circ \widetilde{\mathbf B})^{-\text{vol}},
```

where `b = semi_minor(crack)` is the in-plane semi-axis. Limits :

| Interface | Behavior                                  |
|-----------|--------------------------------------------|
| `Rn = Rt = nothing`   | traction-free penny — `B̃` unchanged |
| `Rn, Rt → 0`          | recovers traction-free                |
| `Rn, Rt → ∞` (rigid)  | `B̃_eff → 0`, cracks behave as bonded |

### 5.3 Notes on Mori-Tanaka and Self-Consistent for cracks

Two well-established but **distinct** formulations co-exist in the
literature for crack-bearing RVEs :

* **Additive (Budiansky-O'Connell)** : MT adds the crack stiffness
  contribution `(4π/3)·ε·Ñ` to the numerator with a zero
  contribution to the denominator (cracks have no volume in the
  strain-concentration sum).  At convergence SC writes
  `J_eff = J_M + ε·H̃(C_eff)`.  This is what MFH currently
  implements.

* **Multiplicative (ECHOES)** : MT adds `(4π/3)·ε·H̃·C_0` to the
  *denominator* via the strain-strain concentration tensor
  (`strain_Strain = H̃·C_0`).  At convergence SC writes
  `C_eff = (B_E)·(A_E)^{-vol}` where the cracks contribute to `A_E`.

The two formulations differ at finite crack density (they coincide in
the dilute limit `ε → 0`).  At `d = 0.30, traction-free` for example,
MFH MT gives `ε_xx(t→∞) ≈ 0.481` while ECHOES MT gives `0.559`.
PCW happens to coincide between the two implementations (at least
numerically through the configurations exercised by
`scripts/60_alv_cracks_interface.jl`).

The MFH MT and SC additive forms are kept for internal consistency
with the (also-additive) MFH elastic MT.  Switching all four MFH
schemes (elastic, conduction, ALV, ALV-cracks) to the multiplicative
ECHOES form simultaneously is left to a follow-up PR — it requires a
coordinated refactor of `Schemes/contribution_helpers.jl`.

ECHOES C++ cross-check : `scripts/60_alv_cracks_interface.jl` runs the
same configuration through both implementations.  At low density the
two MT and SC variants agree numerically with PCW (`rtol ≤ 1e-3`);
at moderate density (`d ≥ 0.20`) the additive vs multiplicative
discrepancy becomes visible (a few % to ~14% depending on density).

A static (non-ageing) elastic + conductivity crack benchmark with
matrix-only interface stiffness is in
`scripts/15_cracks_iso_interface.jl`.

| Scheme                         | Crack treatment                                     |
|--------------------------------|-----------------------------------------------------|
| `Voigt`, `Reuss`               | ignored (zero-volume convention)                    |
| `Dilute`, `DiluteDual`         | additive `+ ΔC_cracks`                              |
| `MoriTanaka`, `Maxwell`, `PCW` | virtual phase with `A = 0`, `N = ΔC`, `f = 1`       |
| `SC`, `ASC`                    | re-evaluated against the running effective estimate |

A complete demo with **all seven** crack-aware ALV schemes lives in
`scripts/57_ageing_creep_cracks.jl`.

## 6. Order-2 ALV — conductivity / diffusion

Same API as the order-4 case, but the kernel returns a 2-tensor:

```julia
function K_iso_order2(t, tp)
    κ = 1.0 + 0.5 * exp(-(t - tp))
    return TensISO{2,3}(κ)
end
law_κ = ViscoLaw(K_iso_order2, :relaxation)

rve_κ = RVE(:M)
add_matrix!(rve_κ, Ellipsoid(1.0), Dict(:K => law_κ))
add_phase!(rve_κ, :I, Ellipsoid(1.0), Dict(:K => heaviside_law(TensISO{2,3}(5.0)));
            fraction = 0.3)

K_eff = homogenize_alv(rve_κ, MoriTanaka(), :K; times = times)   # 150 × 150 (= 3·n)
```

The dispatcher sees the 2-tensor sample and routes via the
order-2 pipeline ([`homogenize_alv_order2`](@ref) under the hood).
Result is a `(3n × 3n)` block matrix. See
`scripts/56_ageing_creep_order2.jl`.

## 7. Symmetry-class fast paths

When all phases share an iso / TI / ortho symmetry with compatible
axes, [`homogenize_alv`](@ref) automatically routes through a fast
path that solves the scheme algebra in the **structured** domain :

| Path  | Components       | Storage      | Volterra inverse cost                |
|-------|------------------|--------------|--------------------------------------|
| ISO   | (α, β)           | 2 n²         | 2 × scalar `n × n` forward solves    |
| TI    | (ℓ₁, …, ℓ₆)      | 6 n²         | (2n × 2n) block + 2 scalar solves    |
| ORTHO | (o₁, …, o₁₂)     | 12 n²        | (3n × 3n) block + 3 scalar solves    |

Detection is heuristic (`_is_iso_block` / `_is_ti_block` /
`_is_ortho_block`) — the user never asks for a fast path explicitly,
and the output is still a dense `(6n × 6n)` `Matrix{T}`.

For user code that wants to keep the compact storage and the type
information, the structured wrappers
[`ALVKernelISO`](@ref) / [`ALVKernelTI`](@ref) /
[`ALVKernelOrtho`](@ref) are `AbstractMatrix{T}` subtypes:

```julia
M = trapezoidal_matrix(law_M, times)
K_iso = ALVKernelISO(M)            # extracts (α, β), 18× cheaper storage

# Algebra closure stays in the structured class (no (6n × 6n)
# materialisation), with auto-promotion iso ⊂ TI ⊂ ortho
K_prod = K_iso * K_iso             # ALVKernelISO
K_inv  = volterra_inverse(K_iso)   # ALVKernelISO

K_TI = ALVKernelTI(K_iso)          # promote to TI form
K_O  = ALVKernelOrtho(K_iso)       # promote to ortho form
K_iso + K_TI                       # ALVKernelTI (auto-promote)
K_iso * K_O                        # ALVKernelOrtho

Matrix(K_iso)                      # back to dense (6n × 6n) on demand
```

These types are a **prototype**: they are fully usable for hand-rolled
ALV pipelines but `homogenize_alv` does not yet accept them as inputs
(use `Matrix(K)` to cross the boundary). See
`scripts/58_alv_kernel_types.jl` for a runnable demo.

## 8. Sensitivities (autodiff via ForwardDiff)

The pipeline supports `ForwardDiff.Dual` end-to-end so derivatives of
effective properties wrt RVE parameters are direct.

### 8.1 Sensitivity wrt volume fraction — recommended `set_param` lens

```julia
using ForwardDiff

# Build the RVE once with a Float64 placeholder fraction.
rve_base = RVE(:M)
add_matrix!(rve_base, Ellipsoid(1.0), Dict(:C => law_M))
add_phase!(rve_base, :I, Ellipsoid(1.0), Dict(:C => heaviside_law(TensISO{3}(60.0, 20.0)));
            fraction = 0.20)

# Differentiate by substituting a `Dual` value via `set_param`.
function eff_mu(f)
    rve_f = set_param(rve_base, AmountParameter(:I), f)
    R̃ = homogenize_alv(rve_f, MoriTanaka(), :C; times = times)
    _, β = iso_params_from_blocks(R̃)
    return β[end, end] / 2
end

dμ_df = ForwardDiff.derivative(eff_mu, 0.20)        # ≈ 1.66 (validated FD ≤ 1e-7)
```

### 8.2 Sensitivity wrt a material parameter — closure-captured

When the parameter lives **inside** the kernel function (e.g. a
modulus, relaxation time, ageing exponent), close it into the kernel
and differentiate normally. ForwardDiff lifts the parameter to `Dual`
through the closure:

```julia
function eff_mu_vs_μM(μ_M)
    function R(t, tp)
        TensISO{3}(15.0, 2 * μ_M * (0.5 + 1.5 * exp(-(t - tp) / 0.5)))
    end
    rve = RVE(:M)
    add_matrix!(rve, Ellipsoid(1.0), Dict(:C => ViscoLaw(R, :relaxation)))
    add_phase!(rve, :I, Ellipsoid(1.0), Dict(:C => heaviside_law(TensISO{3}(60.0, 20.0)));
                fraction = 0.20)
    R̃ = homogenize_alv(rve, MoriTanaka(), :C; times = times)
    _, β = iso_params_from_blocks(R̃)
    return β[end, end] / 2
end

dμ_dμM = ForwardDiff.derivative(eff_mu_vs_μM, 1.0)
```

### 8.3 Joint gradient over multiple parameters

```julia
function eff_mu_vs_p(p)
    f, k_M, μ_M = p
    rve = RVE(:M)
    add_matrix!(rve, Ellipsoid(1.0), Dict(:C => maxwell_iso(k_M, μ_M, 1.0, 0.5)))
    add_phase!(rve, :I, Ellipsoid(1.0), Dict(:C => heaviside_law(TensISO{3}(60.0, 20.0)));
                fraction = 0.20)
    rve_f = set_param(rve, AmountParameter(:I), f)
    R̃ = homogenize_alv(rve_f, MoriTanaka(), :C; times = times)
    _, β = iso_params_from_blocks(R̃)
    return β[end, end] / 2
end

∇ = ForwardDiff.gradient(eff_mu_vs_p, [0.20, 1.0, 1.0])    # 3-vector
```

The complete suite of sensitivity patterns lives in
`scripts/59_alv_sensitivities.jl`. Each derivative is validated against
a central finite difference at `rtol ≤ 1e-7`.

## 9. Validation against ECHOES C++

`scripts/53_ageing_creep_solid.jl` reproduces the multi-phase Maxwell +
solidifying Maxwell + pore benchmark from the ECHOES C++ manual.
`scripts/57_ageing_creep_cracks.jl` covers all seven crack-aware ALV
schemes on a penny-crack RVE.
`scripts/52_rabotnov_mittag_leffler.jl` validates the **Rabotnov /
Mittag-Leffler** closed-form benchmark of @barthelemyIJES2019 §5,
overlaying the analytical curves and reaching `rtol ≤ 1.3e-3` at
`n_times = 200` (trapezoidal-rule discretization accuracy).

The ECHOES Python module is callable from Julia via PyCall :

```julia
using PyCall
ml_dir = raw"<ECHOES_root>\tests\python\creep\mittag_leffler"
pushfirst!(PyVector(pyimport("sys")."path"), ml_dir)
ml_mod = pyimport("mittag_leffler")
I_Rabotnov(t, α, β) = Float64(ml_mod.I_Rabotnov(t, α, β)[])
# … then use `I_Rabotnov` inside a Julia `ViscoLaw` closure
```

Random-RVE cross-checks vs the C++ reference live in
`scripts/bench_echoes/benchmark.jl` (relative error `≤ 1e-8` on the
Mandel `(1, 1)` block, `≤ 1e-6` on the full matrix).
