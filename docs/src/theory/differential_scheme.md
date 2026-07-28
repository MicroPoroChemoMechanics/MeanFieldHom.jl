# The differential scheme

The differential (or incremental) scheme builds the composite by
*repeated dilute incorporation*: an infinitesimal amount of each
inclusion phase is added to the current effective medium, which is
re-homogenized before the next increment. This page derives the
resulting ODE for an arbitrary number of phases and an arbitrary
incorporation trajectory ([Norris 1985](@cite norris1985)).

## Incorporation process

The RVE holds a matrix (index `0`) and `N` inclusion phases, of final
volume fractions ``f_i^\infty``. Each phase is grown along a fictitious
incorporation time,

```math
f_i : [0, T] \to [0, f_i^\infty], \qquad f_i(0) = 0, \quad f_i(T) = f_i^\infty,
```

non-decreasing, the matrix taking up the rest,

```math
f_0(t) = 1 - \sum_{i=1}^{N} f_i(t) .
```

Between ``t`` and ``t + \mathrm dt``, a volume fraction
``\mathrm d\varphi_i`` of the **current medium** is *replaced* by phase
`i`. The ``\mathrm d\varphi_i`` are the increments the scheme actually
applies; the ``\mathrm d f_i`` are what the user prescribes. They
differ, because removing ``\sum_j \mathrm d\varphi_j`` of the current
medium also removes the phase `i` already present in it:

```math
\mathrm d f_i = \mathrm d\varphi_i - f_i \sum_{j=1}^{N} \mathrm d\varphi_j ,
\qquad 1 \le i \le N .
```

In matrix form, with ``\mathbf U = (1, \dots, 1)^{\mathsf T}``,

```math
[\mathrm d f] = \big(\mathbb 1 - [f]\,\mathbf U^{\mathsf T}\big)\,[\mathrm d\varphi] .
```

The matrix to invert is a rank-one update of the identity, so the
Sherman-Morrison formula gives the inverse in closed form:

```math
[\mathrm d\varphi] = \Big(\mathbb 1 + \frac{[f]\,\mathbf U^{\mathsf T}}{1 - \mathbf U^{\mathsf T}[f]}\Big)[\mathrm d f] ,
```

that is, component-wise,

```math
\dot\varphi_i = \dot f_i + \frac{f_i}{1 - \sum_j f_j} \sum_{j=1}^{N} \dot f_j
             = \dot f_i - \frac{f_i}{f_0}\,\dot f_0 .
```

This is exactly what the RHS of the scheme evaluates at every step
(`_diff_ode_rhs!` in `src/Schemes/differential.jl`); nothing else about
the trajectory enters the volume balance.

The derivation above assumes every phase occupies a finite volume.
Crack families do not, and the balance has to be extended — see the
section on cracks below.

## The ODE

The incremental step replaces ``\mathrm d\varphi_i`` of the current
medium by phase `i`, each increment being a *dilute* problem posed in
that medium:

```math
\mathbb C^{hom}(t + \mathrm dt) = \mathbb C^{hom}(t)
   + \sum_{i=1}^{N} \mathrm d\varphi_i \, \mathbb N_i\big(\mathbb C^{hom}(t)\big) ,
```

where ``\mathbb N_i`` is the size-independent **stiffness contribution**
of phase `i` in the current medium,

```math
\mathbb N_i = \mathbb B_i - \mathbb C^{hom} : \mathbb A_i ,
\qquad
\langle \boldsymbol\sigma \rangle_i = \mathbb B_i : \mathbf E^\infty ,
\quad
\langle \boldsymbol\varepsilon \rangle_i = \mathbb A_i : \mathbf E^\infty .
```

Written this way, ``\mathbb N_i`` needs only the phase's concentration
tensors — **not** a Hill tensor, and not even a uniform phase property.
It is therefore available for every inclusion family the package
supports (ellipsoids, cylinders, flat cracks, multi-layer spheres and
spheroids), which is why the differential scheme accepts all of them
without special-casing. See [Localization](localization.md).

Hence the ODE integrated by [`DifferentialScheme`](@ref):

```math
\frac{\mathrm d \mathbb C^{hom}}{\mathrm d t}
  = \sum_{i=1}^{N} \dot\varphi_i \, \mathbb N_i(\mathbb C^{hom})
  = \sum_{i=1}^{N} \Big(\dot f_i + \frac{f_i \sum_j \dot f_j}{1 - \sum_j f_j}\Big)\,
    \mathbb N_i(\mathbb C^{hom}) ,
```

started from ``\mathbb C^{hom}(0) = \mathbb C_0``. The package
parametrizes the incorporation time as ``\tau = t/T \in [0, 1]``.

### Compliance form

The same process written on the compliance uses the **compliance
contribution** ``\mathbb H_i`` of each phase:

```math
\frac{\mathrm d \mathbb S^{hom}}{\mathrm d t}
  = \sum_{i=1}^{N} \dot\varphi_i \, \mathbb H_i(\mathbb S^{hom}) ,
\qquad
\mathbb H_i = - \mathbb S^{hom} : \mathbb N_i : \mathbb S^{hom} .
```

The two forms are *exactly* equivalent: the second is the image of the
first under ``\mathbb S = \mathbb C^{-1}``, term by term. They are not
numerically identical, though — the solver controls its error on
whichever variable is integrated. The compliance form is better
conditioned for a medium that softens towards percolation (porous or
cracked), the stiffness form for one that stiffens. Both are available
through the `formulation` keyword and can be compared on the same RVE:

```julia
homogenize(rve, DifferentialScheme(), :C)                             # stiffness
homogenize(rve, DifferentialScheme(; formulation = :compliance), :C)  # compliance
```

### Cracks: a departure from the volume-replacement picture

Flat cracks do not fit the derivation above, and the difference is
physical rather than technical. A penny-shaped crack of radius `a` and
aperture `c` has, for an aspect ratio ``X = c/a \to 0``,

```math
f_c = \frac{4\pi}{3}\,\varepsilon_c X \;\longrightarrow\; 0 ,
\qquad
\varepsilon_c = \frac{\mathcal N_c\, a^3}{V} = \mathcal O(1) ,
```

(``\mathcal N_c`` cracks of that family in the volume ``V``)

so its **volume fraction vanishes asymptotically while its density stays
finite** — and it is the density, not the fraction, that measures the
mechanical effect. Three consequences for the scheme:

1. **The volume balance must be rewritten.** Opening a crack removes no
   current medium, so cracks do not appear in the sum
   ``\sum_j \mathrm d\varphi_j`` and

   ```math
   f_0 = 1 - \sum_{i \in \mathcal S} f_i
   ```

   whatever the crack densities (``\mathcal S`` = solid phases). They are
   nevertheless *diluted* by it: the piece of current medium that gets
   replaced by solid material takes the cracks it contained with it, and
   the incoming material has none. Writing ``\mathrm d\varphi_c^\varepsilon``
   for the density of cracks actually created in the current medium,

   ```math
   \mathrm d f_i = \mathrm d\varphi_i - f_i \sum_{j \in \mathcal S} \mathrm d\varphi_j ,
   \qquad
   \mathrm d \varepsilon_c = \mathrm d\varphi_c^\varepsilon
                           - \varepsilon_c \sum_{j \in \mathcal S} \mathrm d\varphi_j .
   ```

   Both lines are the single rank-one relation
   ``[\mathrm d g] = (\mathbb 1 - [g]\,\mathbf U_{\mathcal S}^{\mathsf T})[\mathrm d\varphi]``
   on the stacked amounts ``[g] = (f_i ; \varepsilon_c)``, where
   ``\mathbf U_{\mathcal S}`` carries a `1` on solid entries and a `0` on
   crack entries — the *only* change to the manuscript's derivation. The
   same Sherman-Morrison inverse follows, with the density playing the
   role of the fraction:

   ```math
   \dot\varphi_i = \dot f_i + \frac{f_i}{f_0} \sum_{j \in \mathcal S} \dot f_j ,
   \qquad
   \dot\varphi_c^\varepsilon = \dot\varepsilon_c
      + \frac{\varepsilon_c}{f_0} \sum_{j \in \mathcal S} \dot f_j ,
   ```

   and the crack term of the ODE is

   ```math
   \frac{\mathrm d \mathbb C^{hom}}{\mathrm d \tau} \mathrel{+}=
      \sum_c \dot\varphi_c^\varepsilon \,
             \Delta\mathbb C^{crack}_c(\mathbb C^{hom}) .
   ```

   The crack correction vanishes whenever no solid phase grows at the
   same `τ` — a crack-only RVE, or a [`Sequential`](@ref) trajectory in
   which the cracks come after the solids — and is what makes
   ``\varepsilon_c(\tau)`` mean the density *actually reached* at `τ`,
   consistently with ``f_i(\tau)`` for solids.

2. **There is no saturation from the volume side.** A solid phase is
   bounded by ``\sum_i f_i < 1``; a crack density is not, so nothing
   geometric stops the integration — `f₀` is still 1 at any density.
   The effective stiffness decays towards zero asymptotically in
   ``\varepsilon_c`` without a finite threshold, unlike the
   self-consistent scheme, because each infinitesimal crack increment is
   introduced into an already-degraded medium.

3. **The compliance form is the natural one.** As ``X \to 0`` the
   crack's strain concentration tensor diverges, so the stiffness
   contribution ``\mathbb N_c`` is a limit of a divergent quantity,
   while the compliance contribution ``\mathbb H_c`` — built from the
   crack-opening-displacement tensor, see [COD tensors](cod_tensors.md)
   — stays finite and is what the package actually computes:

   ```math
   \Delta \mathbb S^{crack}_c = \frac{4\pi}{3}\,\varepsilon_c\,\mathbb H_c ,
   \qquad
   \Delta \mathbb C^{crack}_c = -\,\mathbb C^{hom} : \Delta\mathbb S^{crack}_c : \mathbb C^{hom} ,
   ```

   with the Budiansky-O'Connell prefactor (`4π/3` for 3D penny and
   elliptic cracks, `π` for 2D ribbons). The stiffness form of the ODE
   therefore evaluates the crack term by passing through `ℍ_c` anyway;
   `formulation = :compliance` merely stops undoing that conversion,
   which is why it behaves better on heavily cracked media.

Apart from that, cracks follow the chosen trajectory exactly like solid
phases, through ``\varepsilon_c(\tau)``: `Proportional`, `Sequential`
and the explicit paths all apply to a `CrackDensity` amount, the target
being the final density instead of the final volume fraction.

## Trajectories

With a single inclusion phase the trajectory is immaterial: any
monotone ``f_1(\tau)`` traverses the same states, and only the final
fraction matters. With two or more phases, the *order of incorporation*
is a modeling choice — the effective property at `τ = 1` depends on it.
Three cases, following the same order as the derivation above.

### 1. Homothetic growth

All phases grow proportionally, ``f_i(t) = \lambda(t)\, f_i^\infty``
with ``\lambda(0) = 0``, ``\lambda(T) = 1`` (the default
[`Proportional`](@ref), where ``\lambda(\tau) = \tau``). Then
``f_0 = 1 - \lambda(1 - f_0^\infty)`` and

```math
\dot\varphi_i = \frac{\dot\lambda\, f_i^\infty}{1 - \lambda\,(1 - f_0^\infty)} ,
\qquad
\varphi_i = -\ln\!\big(1 - \lambda(1 - f_0^\infty)\big)\,
             \frac{f_i^\infty}{1 - f_0^\infty} .
```

This closed form is a useful check on the volume balance alone,
independently of any tensor algebra; the single-phase case reduces to
the classical ``\varphi = -\ln(1 - f)``.

### 2. Successive phases

Each phase is grown to completion before the next one starts —
[`Sequential`](@ref), given the incorporation order. Phase `i` owns a
contiguous slice of `τ` and is frozen outside it.

### 3. Arbitrary trajectory

The user prescribes ``f_i(\tau)`` for every inclusion phase and
``f_0 = 1 - \sum_i f_i`` follows — [`Path`](@ref) for callables
(differentiated by `ForwardDiff`) or [`CustomPath`](@ref) for tabulated
values (piecewise-linear):

```julia
DifferentialScheme(; trajectory = Path(:I1 => τ -> τ^2, :I2 => τ -> 2τ - τ^2))
```

## Conduction and viscoelasticity

The derivation never used the tensor order: replacing
``(\mathbb C, \mathbb N)`` by ``(\mathbf K, \mathbf N_K)`` gives the
conduction / diffusion scheme, and the compliance form becomes the
resistivity form. Both orders are implemented.

In ageing linear viscoelasticity the same ODE holds on the discrete
Volterra block matrices, the products being Volterra products
(`differential_alv` for the relaxation tensor, `differential_alv_order2`
for conduction) — see [Viscoelasticity](viscoelasticity.md). One
restriction is specific to ALV: the ALV Hill kernel is built for an
isotropic reference, and the reference of the differential scheme is its
*running* medium, so every phase must keep that medium isotropic
(spherical inclusions, multi-layer spheres, or any shape with
`symmetrize = :iso`). The ODE raises an explicit error otherwise.

## Numerical resolution

The ODE is integrated by the SciML stack
(`OrdinaryDiffEq.solve`), with an adaptive 5th-order Runge-Kutta
(`Tsit5`) by default. The algorithm and its parameters are user-facing:

```julia
DifferentialScheme(; alg = Vern9(), abstol = 1e-12, reltol = 1e-10)
DifferentialScheme(; maxiters = 10^7)     # any other kwarg reaches `solve`
```

The ODE state is the canonical components of the running estimate, in
the smallest symmetry class that estimate can stay in — two numbers for
an isotropic medium, five to eight for a transversely isotropic one, the
full Mandel matrix otherwise. That class is not the matrix's: a phase
whose contribution is less symmetric drags the running estimate out of
it at the first step (a crack, or simply an aligned spheroid, whose
``\mathbb A_\mathrm{dil}`` is transversely isotropic even between two
isotropic materials). The scheme therefore probes each phase's
contribution once before integrating and sizes the state accordingly.

Two configurations remain out of reach, for the same backend reason:
aligned *triaxial* ellipsoids, and a solid phase combined with an
*aligned* crack family. In both the running estimate is fully
anisotropic from the first step, so the solid Hill tensors have to come
from the numerical residue algorithm — which is degenerate exactly at
the isotropic starting point. `Dilute` and `MoriTanaka` fail the same
way given such a reference, so this is a property of that backend rather
than of the scheme; the scheme reports it with an explicit error.
Giving the phase an orientation average (`symmetrize = :iso`) keeps the
running medium isotropic and resolves both cases.

`nsteps` does **not** set the integration step — it only sets the
density of saved points along `τ`, which [`differential_path`](@ref)
returns:

```julia
τ, Cs = differential_path(rve, DifferentialScheme(; nsteps = 200), :C)
```

Implicit algorithms must be given a non-AD Jacobian
(`Rosenbrock23(autodiff = AutoFiniteDiff())`): the RHS calls the
Hill-tensor backends, which are not differentiable with respect to the
ODE state.
