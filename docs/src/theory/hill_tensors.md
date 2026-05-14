# Hill polarisation tensors

This page reproduces the Hill-tensor framework of the
[Echoes manual](https://github.com/jeanfrancoisbarthelemy/echoes) —
chapter *Eshelby and Hill polarization tensors* and appendix
*Hill polarization tensors* — and points to the corresponding
MeanFieldHom implementation. All expressions, conventions and
references are aligned verbatim on the Echoes Quarto book; MFH
adds a cylinder extension and a 2-D isotropic analytical path that are
flagged as such.

## The Eshelby inclusion problem

Let ``\mathbb R^3`` be filled with a homogeneous linear elastic medium
of stiffness ``\mathbb C``. An ellipsoid ``\mathcal E_{\mathbf A}``
centered at the origin is described by a second-order invertible
tensor ``\mathbf A`` such that ``\mathbf A^{\!T}\!\cdot\mathbf A`` is
symmetric positive-definite, with eigenvalues
``\rho_1=a\ge\rho_2=b\ge\rho_3=c`` (semi-axes) and orthonormal
eigenvectors ``(\hat{\mathbf e}_i^{\mathbf A})``:

```math
\mathbf x\in\mathcal E_{\mathbf A}
\;\Longleftrightarrow\;
\mathbf x\cdot(\mathbf A^{\!T}\!\cdot\mathbf A)^{-1}\!\cdot\mathbf x\le 1,
\qquad
\mathbf A^{\!T}\!\cdot\mathbf A
=\sum_{i=1}^{3}\rho_i^{2}\,\hat{\mathbf e}_i^{\mathbf A}\otimes\hat{\mathbf e}_i^{\mathbf A}.
```

A uniform polarisation stress ``\boldsymbol\tau`` prescribed inside
``\mathcal E_{\mathbf A}`` (and zero outside) drives the boundary
value problem. Following [eshelby1957](@cite), the strain solution is
uniform inside the ellipsoid and reads

```math
\forall\mathbf x\in\mathcal E_{\mathbf A},\qquad
\boldsymbol\varepsilon(\mathbf x) = -\mathbb P:\boldsymbol\tau,
```

where ``\mathbb P=\mathbb P(\mathbf A,\mathbb C)`` is the **Hill
polarisation tensor**. Introducing the equivalent eigenstrain
``\boldsymbol\varepsilon^\star = -\mathbb C^{-1}:\boldsymbol\tau``
gives

```math
\boldsymbol\varepsilon(\mathbf x)=\mathbb S:\boldsymbol\varepsilon^\star,
\qquad
\mathbb S = \mathbb P:\mathbb C,
```

the classical **Eshelby tensor** form. A dual statement involves the
**second Hill tensor** ``\mathbb Q``:

```math
\boldsymbol\sigma(\mathbf x) = -\mathbb Q:\boldsymbol\varepsilon^\star,
\qquad
\mathbb Q = \mathbb C - \mathbb C:\mathbb P:\mathbb C.
```

## Newton-potential integrals

Three geometric integrals — depending only on ``\mathbf A`` — factor
every analytical Hill formula:

```math
\mathbf I^{\mathbf A}
= \frac{\det\mathbf A}{4\pi}\!\int_{\|\hat{\boldsymbol\xi}\|=1}
\frac{\hat{\boldsymbol\xi}\otimes\hat{\boldsymbol\xi}}
     {\|\mathbf A\cdot\hat{\boldsymbol\xi}\|^{3}}\,\mathrm dS_{\xi}
```

```math
\mathbb U^{\mathbf A}
= \frac{\det\mathbf A}{4\pi}\!\int_{\|\hat{\boldsymbol\xi}\|=1}
\frac{\hat{\boldsymbol\xi}\otimes\hat{\boldsymbol\xi}\otimes
      \hat{\boldsymbol\xi}\otimes\hat{\boldsymbol\xi}}
     {\|\mathbf A\cdot\hat{\boldsymbol\xi}\|^{3}}\,\mathrm dS_{\xi}
```

```math
\mathbb V^{\mathbf A}
= \frac{\det\mathbf A}{4\pi}\!\int_{\|\hat{\boldsymbol\xi}\|=1}
\frac{\hat{\boldsymbol\xi}\stackrel{s}{\otimes}\mathbf 1
      \stackrel{s}{\otimes}\hat{\boldsymbol\xi}}
     {\|\mathbf A\cdot\hat{\boldsymbol\xi}\|^{3}}\,\mathrm dS_{\xi}
\;=\;
\tfrac{1}{2}\bigl(\mathbf 1\,\underline{\boxtimes}\,\mathbf I^{\mathbf A}
                 +\mathbf I^{\mathbf A}\,\underline{\boxtimes}\,\mathbf 1\bigr).
```

They are exposed through [`tens_IA`](@ref), [`tens_UA`](@ref) and
[`tens_VA`](@ref). The intrinsic change-of-variable linking the
``\hat{\boldsymbol\xi}`` and
``\hat{\boldsymbol\zeta}=\mathbf A\cdot\hat{\boldsymbol\xi}/\|\cdot\|``
parameterisations is detailed in [barthelemyIJSS2016](@cite).

### Diagonal coefficients ``I_i^{\mathbf A}``, ``I_{ij}^{\mathbf A}``

By symmetry, ``\mathbf A`` and ``\mathbf I^{\mathbf A}`` share their
eigenvectors:

```math
\mathbf I^{\mathbf A}
= \sum_{i=1}^{3} I_i^{\mathbf A}\,
  \hat{\mathbf e}_i^{\mathbf A}\otimes\hat{\mathbf e}_i^{\mathbf A}.
```

The coefficients ``I_i^{\mathbf A}`` and the secondary coefficients
``I_{ij}^{\mathbf A}`` — identified with Newton-potential integrals
[kellogg1929](@cite), [eshelby1957](@cite), [parnell2016](@cite) and
re-written in the generic form of [barthelemyIJSS2016](@cite),
[barthelemyIJES2020_hilltrans](@cite) — admit closed forms in each
symmetry class (triaxial, prolate, oblate, sphere, cylinder). They are
tabulated in the Echoes appendix; for the triaxial case they involve
the complete elliptic integrals ``\mathcal F(\theta,\kappa)`` and
``\mathcal E(\theta,\kappa)`` of first and second kind
[abramowitz1972](@cite), with

```math
\theta = \arcsin\sqrt{1-\tfrac{c^{2}}{a^{2}}},
\qquad
\kappa = \sqrt{\tfrac{a^{2}-b^{2}}{a^{2}-c^{2}}}.
```

The following identities are always satisfied [eshelby1957](@cite):

```math
\sum_i I_i^{\mathbf A} = 1,
\qquad
3\,I_{ii}^{\mathbf A} + \sum_{j\ne i} I_{ij}^{\mathbf A}
= \frac{1}{\rho_i^{2}},
\qquad
3\,\rho_i^{2}\,I_{ii}^{\mathbf A} + \sum_{j\ne i}\rho_j^{2}\,I_{ij}^{\mathbf A}
= 3\,I_i^{\mathbf A}.
```

In MFH, the closed forms are evaluated by
[`MeanFieldHom.Core.newton_potential_3d`](@ref) and
[`MeanFieldHom.Core.newton_potential_2d`](@ref), with the elliptic
integrals [`ell_K`](@ref) / [`ell_E`](@ref) provided by the
`MeanFieldHom.Elliptic` submodule (Carlson symmetric forms
[carlson1995](@cite)).

### Kelvin–Mandel components of ``\mathbb U^{\mathbf A}``, ``\mathbb V^{\mathbf A}``

In the principal frame the generic expressions read

```math
U^{\mathbf A}_{iiii} = \tfrac{3}{2}\bigl(I_i^{\mathbf A}-\rho_i^{2}\,I_{ii}^{\mathbf A}\bigr),
```

```math
U^{\mathbf A}_{iijj} = U^{\mathbf A}_{ijij} = U^{\mathbf A}_{ijji}
= \tfrac{1}{2}\bigl(I_j^{\mathbf A}-\rho_i^{2}\,I_{ij}^{\mathbf A}\bigr)
= \tfrac{1}{2}\bigl(I_i^{\mathbf A}-\rho_j^{2}\,I_{ij}^{\mathbf A}\bigr)
\quad (i\ne j),
```

```math
V^{\mathbf A}_{iiii} = I_i^{\mathbf A},
\qquad
V^{\mathbf A}_{ijij} = V^{\mathbf A}_{ijji}
= \tfrac{1}{4}\bigl(I_i^{\mathbf A}+I_j^{\mathbf A}\bigr)
\quad (i\ne j).
```

Spherical limit ``\mathbf A = \mathbf 1``:

```math
\mathbb U^{\mathbf 1} = \tfrac{1}{3}\mathbb J + \tfrac{2}{15}\mathbb K,
\qquad
\mathbb V^{\mathbf 1} = \tfrac{1}{3}\mathbb I.
```

## Hill tensor in elasticity

### General expression (Willis 1977 / Mura 1987)

For an arbitrary matrix stiffness ``\mathbb C``, the elastic Hill
polarisation tensor is [willis1977](@cite), [mura1987](@cite)

```math
\mathbb P(\mathbf A,\mathbb C)
= \frac{\det\mathbf A}{4\pi}\!\int_{\|\hat{\boldsymbol\xi}\|=1}
\frac{\hat{\boldsymbol\xi}\stackrel{s}{\otimes}
      \bigl(\hat{\boldsymbol\xi}\cdot\mathbb C\cdot\hat{\boldsymbol\xi}\bigr)^{-1}
      \stackrel{s}{\otimes}\hat{\boldsymbol\xi}}
     {\|\mathbf A\cdot\hat{\boldsymbol\xi}\|^{3}}\,\mathrm dS_{\xi}.
```

Inverting the acoustic tensor
``\hat{\boldsymbol\xi}\cdot\mathbb C\cdot\hat{\boldsymbol\xi}`` is the
source of all computational work.

### Isotropic matrix

With a bulk modulus ``k=E/(3(1-2\nu))``, shear modulus
``\mu=E/(2(1+\nu))`` and Lamé first parameter ``\lambda=k-2\mu/3``, the
isotropic stiffness reads

```math
\mathbb C = 3k\,\mathbb J + 2\mu\,\mathbb K = 3\lambda\,\mathbb I + 2\mu\,\mathbb K.
```

Substituting in the general expression gives [willis1977](@cite)

```math
\mathbb P\bigl(\mathbf A,\,3\lambda\,\mathbb I + 2\mu\,\mathbb K\bigr)
= \frac{1}{\lambda+2\mu}\,\mathbb U^{\mathbf A}
+ \frac{1}{\mu}\,\bigl(\mathbb V^{\mathbf A}-\mathbb U^{\mathbf A}\bigr).
```

For a **sphere** ``\mathbf A = \mathbf 1`` this collapses to the
classical Eshelby closed form:

```math
\mathbb P\bigl(\mathbf 1,\,3k\,\mathbb J + 2\mu\,\mathbb K\bigr)
= \frac{1}{3k+4\mu}\,
\left(\mathbb J + \frac{3(k+2\mu)}{5\mu}\,\mathbb K\right).
```

Implementation:
[`src/Elasticity/hill_3d_iso.jl`](https://codeberg.org/MicroPoroChemoMechanics/MeanFieldHom.jl/src/branch/main/src/Elasticity/hill_3d_iso.jl),
triggered by `method = :auto` when `C₀::TensISO`.

### Anisotropic matrix

When ``\mathbb C`` is arbitrarily anisotropic, no closed form of the
above integral is available in general and one must resort to
numerical cubature [ghahremani1977](@cite), [gavazzi1990](@cite),
[masson2008](@cite). MFH implements two algorithm traits that mirror
the Echoes `NUMINT` / `RESIDUES` options:

- `DECUHR` — the 2-D surface integral is handled by the adaptive
  cubature of [espelid1994](@cite). ForwardDiff-safe.
- `Residue` — the 2-D cubature is reduced to a 1-D quadrature by the
  Cauchy residue theorem applied to the inner ``\varphi`` loop, as
  derived in [masson2008](@cite). Float64 only (the polynomial root
  finder used for the inner sum is not ForwardDiff-compatible).

#### Transversely isotropic matrix coaxial with a spheroid (analytical)

For the special case of a TI matrix whose symmetry axis is parallel to
the spheroid axis, MFH provides a fully analytical closed-form path
based on [barthelemyIJES2020_hilltrans](@cite). The Hill tensor
admits the Walpole-basis decomposition
``\mathbb P = P_1 W_1 + P_2 W_2 + P_3 (W_3+W_4) + P_5 W_5 + P_6 W_6``
whose six coefficients depend on the aspect ratio
``\omega = (\text{axial})/(\text{transverse})`` and the five
independent TI elastic constants
``(C_{1111}, C_{1122}, C_{1133}, C_{3333}, C_{2323})`` through six
elementary integrals (closed-form combinations of `acosh` and complex
square roots).

**Selection** — the dispatcher
[`MeanFieldHom.Core._resolve_algo`](@ref) routes a `TensTI{4}` matrix
combined with a coaxial `Ellipsoid{3, Spherical|Prolate|Oblate}` to the
`Analytical` algorithm trait by default. Coaxiality is detected via
the helper `_ti_coaxial(C₀, ell)`. Non-coaxial spheroids and triaxial
ellipsoids fall back to `Residue` (the default for general
anisotropy).

| Algorithm     | Symbol            | Use case                                            | Speed        | ForwardDiff |
|---------------|-------------------|-----------------------------------------------------|--------------|-------------|
| `Analytical`  | `:auto` (default) | TI matrix coaxial with spheroid                     | O(1)         | Yes         |
| `Residue`     | `:residues`       | Anisotropic, default fallback                       | ~ µs         | No          |
| `DECUHR`      | `:decuhr`         | Anisotropic, ForwardDiff-friendly numerical         | ~ ms         | Yes         |
| `NestedQuadGK`| `:nestedquadgk`   | Historical nested-1D-QuadGK, kept for benchmarking | ~ ms         | Yes         |

For other symmetry classes (orthotropic matrix, non-coaxial TI,
generic anisotropic), analytical paths exist in the literature
([withers1989](@cite), [pouya2000](@cite), [pouya2006](@cite),
[suvorov2002](@cite)) but are not yet implemented in MFH.

### Cylinder (MFH extension)

!!! note "Extension over Echoes"
    MFH exposes an explicit `Cylinder` inclusion type corresponding to
    the limit ``a\to\infty`` of a prolate spheroid with transverse
    semi-axes ``b\ge c>0`` and axis ``\hat{\mathbf e}_1``. The
    expressions below are obtained by passing to that limit in the
    generic triaxial formulas [mura1987](@cite), §11.22.

The Newton-potential coefficients become

```math
I_1^{\text{cyl}} = 0,\qquad
I_2^{\text{cyl}} = \frac{c}{b+c},\qquad
I_3^{\text{cyl}} = \frac{b}{b+c},
```

```math
I_{22}^{\text{cyl}} = \frac{c(2b+c)}{3\,b^{2}(b+c)^{2}},\quad
I_{33}^{\text{cyl}} = \frac{b(b+2c)}{3\,c^{2}(b+c)^{2}},\quad
I_{23}^{\text{cyl}} = \frac{1}{(b+c)^{2}},
```

```math
I_{11}^{\text{cyl}} = I_{12}^{\text{cyl}} = I_{13}^{\text{cyl}} = 0
\quad\text{(with }a^{2} I_{12}^{\text{cyl}}\to I_2^{\text{cyl}},\;
a^{2} I_{13}^{\text{cyl}}\to I_3^{\text{cyl}}\text{ as }a\to\infty).
```

Substituting in the Kelvin–Mandel formulas gives a block-diagonal
``\mathbb U^{\text{cyl}}`` and ``\mathbb V^{\text{cyl}}`` whose first
row and column vanish (``U^{\text{cyl}}_{1ikl} = V^{\text{cyl}}_{1ikl} = 0``).
In the isotropic case

```math
\mathbb P^{\text{cyl}}
= \frac{1}{\lambda+2\mu}\,\mathbb U^{\text{cyl}}
+ \frac{1}{\mu}\,\bigl(\mathbb V^{\text{cyl}}-\mathbb U^{\text{cyl}}\bigr),
\qquad
P^{\text{cyl}}_{1jkl}\equiv 0,
```

— no polarisation is transmitted along the cylinder axis.

For the circular cylinder ``b=c`` the non-zero ``\mathbb V`` components
reduce to

```math
V^{\text{cyl},\circ}_{2222} = V^{\text{cyl},\circ}_{3333}
= V^{\text{cyl},\circ}_{2323} = \tfrac{1}{2},
\qquad
V^{\text{cyl},\circ}_{1313} = V^{\text{cyl},\circ}_{1212} = \tfrac{1}{4}.
```

Implementation:
[`hill_3d_cylinder_iso.jl`](https://codeberg.org/MicroPoroChemoMechanics/MeanFieldHom.jl/src/branch/main/src/Elasticity/hill_3d_cylinder_iso.jl).
For arbitrarily anisotropic matrices the Masson polynomial degenerates
at the cylinder limit (one root at infinity); MFH therefore routes
`Cylinder` + `AbstractTens{4,3}` through a dedicated 1-D quadrature
(`CylinderQuadrature` trait,
[`hill_3d_cylinder_aniso.jl`](https://codeberg.org/MicroPoroChemoMechanics/MeanFieldHom.jl/src/branch/main/src/Elasticity/hill_3d_cylinder_aniso.jl)).

### 2-D plane strain (MFH extension)

!!! note "Extension over Echoes"
    MFH also handles the plane-strain limit ``\hat{\boldsymbol\xi}\in S^{1}``
    with a ``1/(2\pi)`` prefactor in place of ``1/(4\pi)``. The
    isotropic version is analytical; the anisotropic one uses the
    Masson residue reduction on the inner line integral.

## Hill tensor in conductivity

The Eshelby framework extends verbatim to transport problems (heat
conduction, mass diffusion, electric conduction, ...): a uniform
polarisation flux inside the ellipsoid produces a uniform temperature
gradient, and a **2nd-order** Hill polarisation tensor
``\mathbf P(\mathbf A,\mathbf K)`` plays the role of ``\mathbb P``
[willis1977](@cite).

### General expression

For a conductivity tensor ``\mathbf K`` [willis1977](@cite):

```math
\mathbf P(\mathbf A,\mathbf K)
= \frac{\det\mathbf A}{4\pi}\!\int_{\|\hat{\boldsymbol\xi}\|=1}
\frac{\hat{\boldsymbol\xi}\otimes\hat{\boldsymbol\xi}}
     {(\hat{\boldsymbol\xi}\cdot\mathbf K\cdot\hat{\boldsymbol\xi})\,
      \|\mathbf A\cdot\hat{\boldsymbol\xi}\|^{3}}\,\mathrm dS_{\xi}.
```

The 2nd-order Eshelby tensor is ``\mathbf s = \mathbf P\cdot\mathbf K``.

### Isotropic matrix

If ``\mathbf K = K\,\mathbf 1``:

```math
\mathbf P(\mathbf A, K\,\mathbf 1) = \frac{\mathbf I^{\mathbf A}}{K}.
```

For a sphere, ``\mathbf I^{\mathbf 1} = \tfrac{1}{3}\mathbf 1`` gives
``\mathbf P(\mathbf 1, K\,\mathbf 1) = \tfrac{1}{3K}\mathbf 1``,
``\mathbf s(\mathbf 1, K\,\mathbf 1) = \tfrac{1}{3}\mathbf 1``
(Eshelby sphere, independent of ``K``).

### Anisotropic matrix — square-root transformation

Unlike the 4th-order elastic case, the 2nd-order Hill tensor admits a
**closed-form** expression for any matrix anisotropy. Following
[giraudMOM2019](@cite) (equivalent derivation by Green's function in
[barthelemyTIPM2009](@cite)), the square-root ``\mathbf K^{1/2}`` of
``\mathbf K = \sum_i K_i\,\hat{\mathbf e}_i^{\mathbf K}\otimes\hat{\mathbf e}_i^{\mathbf K}``
is introduced:

```math
\mathbf K^{1/2}
= \sum_i\sqrt{K_i}\,\hat{\mathbf e}_i^{\mathbf K}\otimes\hat{\mathbf e}_i^{\mathbf K},
\qquad
\mathbf K^{-1/2} = (\mathbf K^{1/2})^{-1}.
```

The change of variable
``\hat{\boldsymbol\zeta}\to\mathbf K^{1/2}\cdot\mathbf A^{-1}\cdot\hat{\boldsymbol\zeta}/\|\cdot\|``
reduces the general expression to the Newton-potential integral of a
**fictitious ellipsoid** of shape tensor
``\mathbf A\cdot\mathbf K^{-1/2}``:

```math
\mathbf P(\mathbf A,\mathbf K)
= \mathbf K^{-1/2}\cdot\mathbf I^{\mathbf A\cdot\mathbf K^{-1/2}}\cdot\mathbf K^{-1/2}.
```

The semi-axes and principal directions of the fictitious ellipsoid are
obtained by diagonalising
``\mathbf K^{-1/2}\cdot\mathbf A^{\!T}\!\cdot\mathbf A\cdot\mathbf K^{-1/2}``.
Implementation:
[`src/Conductivity/hill_order2_3d.jl`](https://codeberg.org/MicroPoroChemoMechanics/MeanFieldHom.jl/src/branch/main/src/Conductivity/hill_order2_3d.jl).

## Eshelby tensor

The **Eshelby tensor** is the contraction of the Hill polarisation
tensor with the reference-medium stiffness or conductivity:

```math
\mathbb S(\mathbf A,\mathbb C_{0}) = \mathbb P(\mathbf A,\mathbb C_{0}) : \mathbb C_{0}
\qquad\text{(order 4, elasticity)},
```

```math
\mathbf s(\mathbf A,\mathbf K_{0}) = \mathbf P(\mathbf A,\mathbf K_{0})\cdot\mathbf K_{0}
\qquad\text{(order 2, conductivity / diffusion)}.
```

Public entry point [`eshelby_tensor`](@ref MeanFieldHom.Core.eshelby_tensor):

```julia
using MeanFieldHom, TensND

C₀ = TensISO{3}(3k, 2μ)                    # elastic reference
S  = eshelby_tensor(Ellipsoid(1.0), C₀)    # 4th-order

K₀ = TensISO{3}(1.0)                        # conductivity reference
s  = eshelby_tensor(Ellipsoid(1.0), K₀)    # 2nd-order
```

## Dispatch and implementation notes

| `(inclusion, C₀)`                    | `:auto` selects        | alternative(s)             | ForwardDiff |
| :----------------------------------- | :--------------------- | :------------------------- | :---------: |
| `Ellipsoid{3}, TensISO`              | `Analytical`           | —                          |     ✓       |
| `Ellipsoid{3}, TensTI` (aligned)     | `Analytical` (MFH)     | `:residues`, `:decuhr`      |     ✓       |
| `Ellipsoid{3}, AbstractTens{4,3}`    | `Residue` (Float64)    | `:decuhr`                  |  ✓ (decuhr) |
| `Cylinder, TensISO`                  | `Analytical`           | —                          |     ✓       |
| `Cylinder, AbstractTens{4,3}`        | `CylinderQuadrature`   | (residue degenerates)      |     ✓       |
| `Ellipsoid{2}, TensISO`              | `Analytical`           | —                          |     ✓       |
| `Ellipsoid{2}, AbstractTens{4,2}`    | `Analytical` (residue) | —                          |  (Float64)  |
| `Ellipsoid{3}, AbstractTens{2,3}`    | `Analytical` (K⁻¹ᐟ²)   | —                          |     ✓       |
| `Cylinder, AbstractTens{2,3}`        | `Analytical`           | —                          |     ✓       |

Entry point: [`hill_tensor`](@ref). Shape tensor retrieval:
[`shape_tensor`](@ref). Builders of geometric auxiliaries:
[`tens_IA`](@ref), [`tens_UA`](@ref), [`tens_VA`](@ref).
