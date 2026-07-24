# Layered spheroid — confocal harmonic series and imperfect interfaces

`MeanFieldHom.LayeredSpheroids` provides [`LayeredSpheroid`](@ref) — an
`N`-layer confocal spheroidal composite inclusion embedded in an
infinite isotropic matrix, **conduction only** (thermal / electric /
Darcy). It follows
[Barthélémy & Bignonnet (2020)](@cite barthelemyBignonnetIJES2020),
which extends the layered-sphere state-vector recurrence of
[herve1993](@cite) to spheroids: the non-spherical geometry forces the
imperfect (Kapitza / surface-conductive) interfaces to couple
*different* spherical-harmonic degrees, so the sphere's simple 2×2
transfer per mode becomes a truncated series with a `2𝒩 × 2𝒩`
transfer matrix per interface. Unlike [`LayeredSphere`](@ref), there is
no elastic counterpart: the harmonic decomposition is specific to the
scalar Laplace equation.

## Confocal spheroidal coordinates

Prolate spheroidal coordinates `(φ, p, q)`, revolution axis `n̂`, half
focal distance `c > 0`:

```math
x_1 = c\sqrt{1-p^2}\sqrt{q^2-1}\cos\varphi,\quad
x_2 = c\sqrt{1-p^2}\sqrt{q^2-1}\sin\varphi,\quad
x_3 = c\,p\,q,
```

with `-1 ≤ p ≤ 1`, `q ≥ 1`. Iso-`q` surfaces are confocal spheroids of
axis/disk semi-axes `ρ_a = c q`, `ρ_t = c\sqrt{q^2-1}` (aspect ratio
`ω = q/\sqrt{q^2-1} > 1`). **Oblate spheroids** (`ω < 1`) follow by the
formal substitution `c → -i c̄`, `q → i τ` (`c̄, τ` real): every prolate
formula carries over unchanged once evaluated with complex arithmetic —
`ρ_a = c̄ τ`, `ρ_t = c̄\sqrt{τ^2+1}`. `LayeredSpheroid` exploits this
directly: its confocal parameter `q` is stored as `T` (prolate) or
`Complex{T}` (oblate), and every downstream function (Legendre
recurrences, coupling integrals, transfer matrices) is written
generically over `Q <: Number` — no branch is ever taken on
prolate/oblate.

## Boundary value problem and harmonic series

`N` confocal layers of isotropic conductivity `k_ℓ`, interfaces
`ℐ_ℓ` at `q = q_ℓ` (`ℓ = 1,…,N`), embedded in a matrix of
conductivity `k_{N+1}` for `q > q_N`, under a remote uniform gradient
`H = H₁ê₁ + H₃ê₃`. The temperature in layer `ℓ` decomposes as a series
of spheroidal harmonics,

```math
T_\ell = c \sum_{m=0}^\infty \sum_{n=m}^\infty
P_n^m(p)\Big[a_{\ell,n}^m P_n^m(q) + b_{\ell,n}^m Q_n^m(q)\Big]\cos(m\varphi) + (\sin\text{-terms}),
```

`P_n^m`, `Q_n^m` associated Legendre functions of the first/second
kind. By symmetry the problem splits into an **axial** problem
(`H = H₃ê₃`, order `m = 0`, odd degrees only) and a **transverse** one
(`H = H₁ê₁`, order `m = 1`, odd degrees), solved by the exact same
machinery — an arbitrary in-plane `H` follows by rotation (axisymmetry
of the geometry).

## Interface conditions and the coupling matrices `I, J, K, L`

Three interface types couple the axial (or transverse) series:

- **Perfect** (`PerfectInterface`): temperature and normal flux
  continuous, diagonal in degree.
- **LC** (`KapitzaInterface(ρ)`, low-conducting): flux continuous, but
  the temperature jump `[T] = ρ·q_n` mixes ALL degrees through the
  integral `I_{ij}(q) = ∫_{-1}^1 P_i(x)P_j(x)/\sqrt{q^2-x^2}\,dx`
  (axial) or `J_{ij}` (transverse, `Pᵢ → Pᵢ¹`).
- **HC** (`SurfaceConductiveInterface(β)`, highly-conducting):
  temperature continuous, but the flux jump
  `[q_n] = -β\,\mathrm{div}_S(∇_S T)` mixes degrees through `J_{ij}`
  (axial) or `K_{ij} + L_{ij}/(q^2-1)` (transverse).

`KapitzaInterface`/`SurfaceConductiveInterface` reuse
[`LayeredSpheres`](@ref MeanFieldHom.LayeredSpheres)'s types and sign
convention (`ρ` a genuine thermal resistance, `β` a genuine surface
conductance — NOT the inverse convention some raw `echoes`
`interf_prop` values use for the low-conducting case).

Truncating the series at `𝒩` terms (odd degrees `1, 3, …, 2𝒩-1`), the
interface condition becomes a `2𝒩 × 2𝒩` linear map
`X_{ℓ+1} = R_ℓ X_ℓ` between the layers' coefficient vectors
`X_ℓ = [A_ℓ; B_ℓ]`,

```math
R_\ell = \mathcal J(k_{\ell+1}, q_\ell)^{-1}\Big(\mathcal J(k_\ell, q_\ell) + \delta\mathcal J^{\mathrm{LC/HC}}\Big),
\qquad
\mathcal J(k,q) = \begin{pmatrix}\mathcal J_P(q) & \mathcal J_Q(q) \\ k\,\mathcal J_{P'}(q) & k\,\mathcal J_{Q'}(q)\end{pmatrix},
```

`\mathcal J_R(q) = \mathrm{diag}(R_1(q), R_3(q), …, R_{2𝒩-1}(q))`, and
`δ\mathcal J^{LC}`/`δ\mathcal J^{HC}` full `𝒩×𝒩` blocks built from
`I`/`J` (LC) or `J`/`K+L/(q^2-1)` (HC), perturbing the temperature
(LC) or flux (HC) half of `\mathcal J`. Cumulating
`S_ℓ = R_ℓ⋯R_1`, imposing core regularity `B_1 = 0` and the unit
remote field `A_{N+1} = (±1, 0, …, 0)` (`+` axial, `−` transverse)
gives every layer's coefficients — this is exactly what
[`spheroid_state_sequence`](@ref MeanFieldHom.LayeredSpheroids.spheroid_state_sequence)
computes, and what [`local_temperature`](@ref
MeanFieldHom.LayeredSpheroids.local_temperature) /
[`local_gradient`](@ref MeanFieldHom.LayeredSpheroids.local_gradient) /
[`local_flux`](@ref MeanFieldHom.LayeredSpheroids.local_flux)
reconstruct pointwise from.

## Volume-averaged concentration tensors

The whole-particle averages depend only on the ratio
`(b/a) = b^{0/1}_{N+1,1}/a^{0/1}_{N+1,1}` at the OUTER boundary `q_N`
and two pairs of shape functions,

```math
\mathcal T_a(q) = \operatorname{arccoth}q - \tfrac1q,\qquad
\mathcal T_t(q) = \operatorname{arccoth}q - \tfrac{q}{q^2-1},\qquad
\mathcal U_a(q) = \mathcal T_t(q),\qquad
\mathcal U_t(q) = \operatorname{arccoth}q + \tfrac{2-q^2}{q(q^2-1)},
```

giving `⟨∇T⟩_Ω = A_Ω · H`, `⟨K∇T⟩_Ω = B_Ω · H` with (axial/transverse
diagonal in the spheroid's own frame, `TensND.TensTI{2,3}`)

```math
A_\Omega = \mathrm{diag}\big(1+\tfrac ba\big|_t \mathcal T_t(q_N),\; \ldots,\; 1+\tfrac ba\big|_a \mathcal T_a(q_N)\big),
\qquad
B_\Omega = k_{N+1}\,\mathrm{diag}\big(1+\tfrac ba\big|_t \mathcal U_t(q_N),\; \ldots\big).
```

`LayeredSpheroid` has **no Hill tensor** — like `LayeredSphere`, it
plugs into the mean-field schemes through `A_Ω`, `B_Ω` directly
(`gradient_gradient_loc`, `flux_gradient_loc`,
`conductivity_contribution` overrides,
`is_homogeneous_inclusion = false`). The size-independent contribution
satisfies the same invariant as the sphere: `N_K = B_Ω - K₀·A_Ω`.

**Equivalent particle.** [Barthélémy & Bignonnet (2020, §4)](@cite
barthelemyBignonnetIJES2020) define the equivalent homogeneous
(perfectly-bonded) particle by `k^{eq} = B_Ω · A_Ω^{-1}` — a
size-dependent quantity (unlike a homogeneous perfect-interface
spheroid's shape-only response), computed directly from
`gradient_gradient_loc`/`flux_gradient_loc` in `scripts/34_spheroid_equivalent_conductivity.jl`.
When every interface is perfect, only degree 1 survives and `k^{eq}`
follows the closed-form nested recursion of the paper's §3 (built from
the classical conduction depolarization factors already available via
[`tens_IA`](@ref)/[`hill_tensor`](@ref) for a single spheroid) — an
exact oracle used to validate the general series solution in
`test/LayeredSpheroids/test_conductivity.jl`.

## Numerical precision: quadrature vs. the original monomial series

The reference implementation
(`echoes_cpp/interface/python/py_inclusions/spheroid_nlayers.py`)
computes `I, J, K, L` by expanding `Pᵢ(x)Pⱼ(x)` etc. into monomials
(coefficients `γ, η, δ`, built by the recursions of the paper's
Appendix) and summing against `W_k(q) = ∫_{-1}^1 x^{2k}/\sqrt{q^2-x^2}\,dx`
(closed form via `₂F₁`). This is numerically **ill-conditioned**: the
monomial coefficients of a degree-`n` Legendre polynomial grow like
`10^{0.8n}`, while `I_{ii}(q) = O(1/n)`, so the summation
`I_{ij} = \sum_k γ_{2k} W_k(q)` must cancel `10^{0.8n}`-sized terms
down to an `O(1/n)` result. The paper's own rule of thumb — working
precision `\gtrsim 0.8(2𝒩-1)` decimal digits — is exactly this
cancellation bound, and is why the reference implementation resorts to
`mpmath` arbitrary precision once `𝒩 \gtrsim 10`.

`MeanFieldHom.jl`'s default backend
([`coupling_matrices`](@ref MeanFieldHom.LayeredSpheroids.coupling_matrices)`,
`method = :quadrature`) sidesteps the issue rather than reproducing it:
`I, J, K, L` are the paper's own closed-form INTEGRAL definitions
(e.g. `I_{ij}(q) = ∫_{-1}^1 P_i(x)P_j(x)/\sqrt{q^2-x^2}\,dx`), integrated
directly by Gauss quadrature (`QuadGK`) with `Pᵢ(x)`, `Pᵢ¹(x)` and
their derivatives evaluated at real `x ∈ [-1,1]` by the STABLE
three-term recurrence (never the monomial expansion). For prolate `q`
the integrand is smooth and bounded on `[-1,1]`; for oblate `q = iτ` it
is complex-analytic with no real singularity. Either way `Float64`
quadrature converges to essentially machine precision for any `𝒩` —
`scripts/33_spheroid_series_convergence.jl` demonstrates this directly
against the faithful BigFloat port of the original monomial series
(`method = :series`, kept as an independent validation oracle and
cross-checked to machine precision in
`test/LayeredSpheroids/test_coupling.jl`).

## Registering a new phase type — the pattern this follows

`LayeredSpheroid` is integrated with the schemes exactly like
`LayeredSphere` (see [layered_sphere.md](@ref) and
`src/LayeredSpheres/scheme_integration.jl` for the elastic/conduction
analog): `is_homogeneous_inclusion = false`, and overrides of
`gradient_gradient_loc`, `flux_gradient_loc`,
`conductivity_contribution` (all 3-argument, the declared phase
property ignored) route every dilute/Mori-Tanaka/self-consistent/
Maxwell/differential kernel through the confocal transfer-matrix
solution with no change to `src/Schemes/`.
