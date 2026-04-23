# Crack opening displacement and compliance tensors

This page follows the Echoes manual chapter *Crack compliance tensors*
(and the chapter on stress/displacement intensity factors) and
points to the MeanFieldHom implementation. All expressions,
conventions and references are aligned on the Echoes Quarto book.

The computational pipeline is deliberately narrow:

```
   B   ──►   H = (3/4)·n̂ ⊗ˢ B ⊗ˢ n̂   ──►   K̂, N̂   (SIF / DIF)
 (COD)       (crack compliance)          (fracture quantities)
```

All quantities reduce to the **COD tensor** ``\mathbf B`` and the
resolved stress ``\boldsymbol\Sigma\cdot\hat{\mathbf n}`` on the
crack plane — no linear transformation on ``\mathbf A``, ``\mathbb C``
is introduced.

## Crack geometry

A flat elliptic crack is the limit of a flat spheroidal inclusion as
its smallest aspect ratio tends to zero. Keeping the Echoes
parameterisation of the ellipsoid by the shape tensor ``\mathbf A``,

```math
\mathbf A
= \hat{\boldsymbol\ell}\otimes\hat{\boldsymbol\ell}
+ \eta\,\hat{\mathbf m}\otimes\hat{\mathbf m}
+ \omega\,\hat{\mathbf n}\otimes\hat{\mathbf n},
\qquad
\eta=\frac{b}{a},\quad\omega=\frac{c}{a}\to 0,
```

with in-plane orthonormal frame
``(\hat{\boldsymbol\ell},\hat{\mathbf m})`` along the major/minor
semi-axes ``a\ge b`` and unit normal ``\hat{\mathbf n}``. MFH
supports:

- **Elliptic cracks** of aspect ratio ``\eta\in(0,1]``
  ([`EllipticCrack`](@ref)) — including the circular **penny**
  ``\eta=1`` ([`PennyCrack`](@ref));
- **Ribbon cracks** (tunnel cracks) infinite along
  ``\hat{\boldsymbol\ell}`` of half-width ``b`` along
  ``\hat{\mathbf m}`` ([`RibbonCrack`](@ref)).

## Crack compliance ``\mathbb H`` and COD tensor ``\mathbf B``

The **COD tensor** ``\mathbf B`` is the size-independent symmetric
2-tensor defined from the averaged displacement jump on the crack
surface:

```math
\frac{1}{S}\int_{S}[\![\mathbf u]\!]\,\mathrm dS
\;=\; b\,\mathbf B\cdot(\boldsymbol\Sigma\cdot\hat{\mathbf n}),
```

where ``b`` is the semi-minor in-plane semi-axis (``b\ge c\to 0``).

The **crack compliance contribution tensor** is defined from the
second Hill tensor ``\mathbb Q = \mathbb C - \mathbb C:\mathbb P:\mathbb C``
through the limit
[kachanov1992](@cite), [kachanov1993](@cite), [sevostianov2002](@cite),
[barthelemyIJES2021](@cite):

```math
\mathbb H \;=\; \lim_{c/b\to 0}\, \frac{c}{b}\,\mathbb Q^{-1}.
```

The limit is finite (the divergent components of ``\mathbb Q^{-1}``
scale as ``b/c`` so that ``(c/b)\mathbb Q^{-1}`` stays bounded). The
prefactor ``c/b`` – rather than the ``c/a`` used in earlier works –
ensures a **uniform definition across 2D and 3D geometries**
(``b`` is the only always-finite semi-axis for both elliptic cracks
with ``c\to 0`` and ribbon cracks with ``a\to\infty``).

The consistency ``\int\boldsymbol\varepsilon\,\mathrm dV =
\int[\![\mathbf u]\!]\stackrel{s}{\otimes}\hat{\mathbf n}\,\mathrm dS``
combined with the COD definition yields a **geometric factorisation of
``\mathbb H``** through ``\mathbf B``:

```math
\mathbb H \;=\; \frac{c\,S}{V}\,
\hat{\mathbf n}\stackrel{s}{\otimes}\mathbf B\stackrel{s}{\otimes}\hat{\mathbf n},
```

with ``S`` the crack surface and ``V`` the embedding ellipsoidal
volume.  Evaluating ``cS/V`` for each geometry:

- **Elliptic 3D crack** (``S=\pi ab``, ``V=\tfrac{4}{3}\pi abc``):
  ``cS/V = 3/4``, so

  ```math
  \mathbb H^{\mathcal E}
  \;=\; \tfrac{3}{4}\,\hat{\mathbf n}\stackrel{s}{\otimes}\mathbf B
                                 \stackrel{s}{\otimes}\hat{\mathbf n}.
  ```

- **Ribbon 2D crack** (``S=4ab``, ``V=2\pi abc``, ``a\to\infty``):
  ``cS/V = 2/\pi``, so

  ```math
  \mathbb H^{\mathcal R}
  \;=\; \tfrac{2}{\pi}\,\hat{\mathbf n}\stackrel{s}{\otimes}\mathbf B
                                 \stackrel{s}{\otimes}\hat{\mathbf n}.
  ```

Consistency of ``\mathbb H`` between the ribbon limit ``\eta\to 0`` of
the elliptic case and the intrinsic 2D derivation fixes the relation

```math
\mathbf B^{\mathcal R} \;=\; \tfrac{3\pi}{8}\,\lim_{\eta\to 0}\mathbf B^{\mathcal E}.
```

``\mathbf B`` is size-independent (it only depends on ``\eta`` and on
the orientation of the crack plane). It is the quantity MFH computes
directly via [`cod_tensor`](@ref) (alias [`B_tensor`](@ref));
the bridge ``\mathbf B\leftrightarrow\mathbb H`` is handled by
[`compliance_from_cod`](@ref) / [`cod_from_compliance`](@ref)
([`src/Cracks/cod_H_bridge.jl`](https://github.com/MicMacTools/MeanFieldHom.jl/tree/main/src/Cracks/cod_H_bridge.jl)),
which dispatches on the crack type (elliptic or ribbon) to apply the
correct geometric factor.

## Analytical COD tensor — isotropic matrix

For an isotropic matrix ``\mathbb C = 3k\,\mathbb J + 2\mu\,\mathbb K``
with Young's modulus ``E`` and Poisson ratio ``\nu``, the Echoes
closed-form expression of ``\mathbf B`` in the crack-local frame
``(\hat{\boldsymbol\ell},\hat{\mathbf m},\hat{\mathbf n})`` is

```math
\mathbf B
= B_{nn}\,\hat{\mathbf n}\otimes\hat{\mathbf n}
+ B_{mm}\,\hat{\mathbf m}\otimes\hat{\mathbf m}
+ B_{\ell\ell}\,\hat{\boldsymbol\ell}\otimes\hat{\boldsymbol\ell},
```

```math
\begin{aligned}
B_{nn}       &= \frac{8\,\eta\,(1-\nu^{2})}{3E}\,\frac{1}{\mathcal E_\eta},\\[6pt]
B_{mm}       &= \frac{8\,\eta\,(1-\nu^{2})}{3E}\,
                \frac{1-\eta^{2}}
                     {\bigl(1-(1-\nu)\eta^{2}\bigr)\mathcal E_\eta
                      - \nu\,\eta^{2}\,\mathcal K_\eta},\\[6pt]
B_{\ell\ell} &= \frac{8\,\eta\,(1-\nu^{2})}{3E}\,
                \frac{1-\eta^{2}}
                     {\bigl(1-\nu-\eta^{2}\bigr)\mathcal E_\eta
                      + \nu\,\eta^{2}\,\mathcal K_\eta},
\end{aligned}
```

where ``\mathcal K_\eta=\mathcal K(\sqrt{1-\eta^{2}})`` and
``\mathcal E_\eta=\mathcal E(\sqrt{1-\eta^{2}})`` are the complete
elliptic integrals of first and second kind [abramowitz1972](@cite).
In MFH these are provided by [`ell_K`](@ref) and [`ell_E`](@ref).

For the **circular penny** ``\eta=1``, ``\mathcal K_1 = \mathcal E_1 = \pi/2``
and the formulas collapse to

```math
B_{nn} = \frac{16\,(1-\nu^{2})}{3\pi E},
\qquad
B_{mm} = B_{\ell\ell} = \frac{B_{nn}}{1-\nu/2}.
```

Implementation:
[`src/Cracks/cod_analytical.jl`](https://github.com/MicMacTools/MeanFieldHom.jl/tree/main/src/Cracks/cod_analytical.jl),
selected by `method = :auto` when `C₀::TensISO{4,3}`.

## Transversely-isotropic matrix

When the matrix is transversely isotropic and the TI axis is
**aligned with the crack normal** ``\hat{\mathbf n}``, Echoes also
handles ``\mathbf B`` analytically. The closed-form expressions
involve the engineering parameters ``(E,\nu_{1},\nu_{2},H,\Gamma)``
defined on the compliance ``\mathbb S=\mathbb C^{-1}``
[hoenig1978](@cite), [kanaun2009](@cite), [barthelemyIJES2021](@cite);
they reduce to the isotropic case for
``\nu_{1}=\nu_{2}=\nu``, ``H=\Gamma=1``.

MFH exposes the TI closed form in
[`src/Cracks/cod_analytical.jl`](https://github.com/MicMacTools/MeanFieldHom.jl/tree/main/src/Cracks/cod_analytical.jl);
the detailed algebraic expressions (auxiliary coefficients
``R_{ijkl}``, ``\sigma_\gamma``) are documented inline in that file to
keep the present theory page close to the Echoes manual.

!!! note "Misaligned TI / elliptic orthotropy — not yet documented"
    The more general cases of a TI matrix whose axis is **not** aligned
    with the crack normal, or of an elliptic-orthotropic matrix, are
    not covered here. They will be handled in a future extension.

## Anisotropic matrix — numerical paths

For an arbitrarily anisotropic matrix there is no closed form of
``\mathbf B`` in general. Following [barthelemyIJSS2009](@cite), the
limit ``\omega\to 0`` is resolved by extracting the
**first-order term** of the Taylor expansion of the Hill polarisation
tensor in ``\omega``; that term admits an integral representation on
the unit circle of the crack plane which MFH evaluates via two
algorithm traits:

- `DECUHR` — adaptive 2-D cubature of [espelid1994](@cite);
  ForwardDiff-safe. Entry point
  [`src/Cracks/green_decuhr.jl`](https://github.com/MicMacTools/MeanFieldHom.jl/tree/main/src/Cracks/green_decuhr.jl).
- `Residue` — Cauchy-residue reduction to a 1-D quadrature, as in
  [masson2008](@cite) adapted to the crack kernel; Float64 only
  (PolynomialRoots). Entry point
  [`src/Cracks/green_residue.jl`](https://github.com/MicMacTools/MeanFieldHom.jl/tree/main/src/Cracks/green_residue.jl).

`method = :auto` picks `Residue` on anisotropic Float64 inputs and
falls back to `DECUHR` for symbolic or `ForwardDiff.Dual` scalars.

## Compliance contribution to the effective stiffness

The size-independent contribution tensor returned by
[`compliance_contribution`](@ref)`(crack, C₀)` is ``\mathbb H`` itself
— not the full dilute correction ``\Delta\mathbb S``.  The dilute
correction is recovered by applying the Budiansky density
``\varepsilon`` through [`delta_compliance`](@ref):

- **Elliptic 3D crack**: Budiansky density
  ``\varepsilon^{3\mathrm d} = N\,a\,b^{2}`` (number density × major ×
  minor² semi-axes), and

  ```math
  \Delta\mathbb S \;=\; \tfrac{4\pi}{3}\,\varepsilon^{3\mathrm d}\,\mathbb H^{\mathcal E}.
  ```

- **Ribbon 2D crack**: Budiansky density
  ``\varepsilon^{2\mathrm d} = N\,b^{2}`` (number per unit area × minor²
  semi-axis), and

  ```math
  \Delta\mathbb S \;=\; \pi\,\varepsilon^{2\mathrm d}\,\mathbb H^{\mathcal R}.
  ```

Implementation: [`compliance_contribution`](@ref) returns ``\mathbb H``
directly (convention shared with Echoes); [`delta_compliance`](@ref)
assembles ``\Delta\mathbb S`` from ``\mathbb H`` and ``\varepsilon``
with dispatch on the crack shape
([`src/Cracks/compliance.jl`](https://github.com/MicMacTools/MeanFieldHom.jl/tree/main/src/Cracks/compliance.jl)).

## Conductivity crack resistivity

The transport analog of the elastic crack compliance is the
**crack resistivity contribution tensor** ``\mathbf R``
[kachanov2018](@cite), assembled from a scalar thermal COD ``b`` and
an effective rank-1 direction ``\hat{\mathbf w}``; see
[Thermal cracks — COD scalar and resistivity contribution](thermal_cracks.md)
for the full derivation.  The geometric factors ``3/4`` (elliptic) and
``2/\pi`` (ribbon) coincide with those of the elasticity case.

## Stress and displacement intensity factors

At a point ``\mathbf x^{\star}_{0}`` of the crack front, with outer
in-plane normal ``\hat{\boldsymbol\nu}`` and tangent
``\hat{\boldsymbol\tau}=\hat{\mathbf n}\wedge\hat{\boldsymbol\nu}``,
the asymptotic expansions of the displacement jump and traction read
[irwin1957](@cite), [kassir1968](@cite), [willis1968](@cite):

```math
[\![\mathbf u]\!](\mathbf x^{\star}_{0}+r\hat{\boldsymbol\nu})
\underset{r\to 0^{-}}{\sim}
8\,\sqrt{\tfrac{-r}{2\pi}}\,\hat{\mathbf N},
\qquad
\mathbf t(\mathbf x^{\star}_{0}+r\hat{\boldsymbol\nu})
\underset{r\to 0^{+}}{\sim}
\frac{\hat{\mathbf K}}{\sqrt{2\pi r}}.
```

``\hat{\mathbf N}`` is the **displacement intensity factor** (DIF)
vector, ``\hat{\mathbf K}`` the **stress intensity factor** (SIF)
vector. Their normalisation is chosen so that the local energy release
rate [barnett1972](@cite), [rice1989](@cite) reads
``G = \hat{\mathbf K}\cdot\hat{\mathbf N}``.

A central identity of the general anisotropic theory
[kanaun1981](@cite), [kunin1983](@cite), [kanaun2009](@cite) is that
the SIF and DIF are purely local and are exchanged by the COD tensor
``\mathbf B^{\mathcal R}`` of the **ribbon crack tangent** to the real
crack at the observation point:

```math
\hat{\mathbf K}
= \pi\,\bigl(\mathbf B^{\mathcal R}(\hat{\boldsymbol\nu},\hat{\mathbf n})\bigr)^{-1}
\cdot\hat{\mathbf N}.
```

This relation holds regardless of the matrix anisotropy and of the
far-field loading.

### Elliptic crack

The crack-plane parameterisation

```math
\hat{\mathbf y}^{\star}_{0}
= \cos\theta_y\,\hat{\boldsymbol\ell} + \sin\theta_y\,\hat{\mathbf m}
```

defines the tip outer normal

```math
\hat{\boldsymbol\nu}
= \frac{\mathbf S^{\dagger}\cdot\hat{\mathbf y}^{\star}_{0}}
       {\|\mathbf S^{\dagger}\cdot\hat{\mathbf y}^{\star}_{0}\|},
```

where ``\mathbf S`` is the 2-D semi-axis tensor and
``\mathbf S^{\dagger}`` its pseudo-inverse. The SIF and DIF vectors at
the front are

```math
\hat{\mathbf N}^{\mathcal E}
= \tfrac{3}{8}\sqrt{\pi b}\,
  \sqrt{b\,\|\mathbf S^{\dagger}\cdot\hat{\mathbf y}^{\star}_{0}\|}\,
  \mathbf B^{\mathcal E}(\hat{\mathbf m},\hat{\mathbf n},\eta)
  \cdot\boldsymbol\Sigma\cdot\hat{\mathbf n},
```

```math
\hat{\mathbf K}^{\mathcal E}
= \tfrac{3}{8}\pi^{3/2}\sqrt{b}\,
  \sqrt{b\,\|\mathbf S^{\dagger}\cdot\hat{\mathbf y}^{\star}_{0}\|}\,
  \bigl(\mathbf B^{\mathcal R}(\hat{\boldsymbol\nu},\hat{\mathbf n})\bigr)^{-1}
  \cdot\mathbf B^{\mathcal E}(\hat{\mathbf m},\hat{\mathbf n},\eta)
  \cdot\boldsymbol\Sigma\cdot\hat{\mathbf n},
```

with the dimensionless prefactor

```math
b\,\|\mathbf S^{\dagger}\cdot\hat{\mathbf y}^{\star}_{0}\|
= \sqrt{\eta^{2}\cos^{2}\theta_y + \sin^{2}\theta_y}.
```

### Ribbon crack

For a ribbon crack with
``\hat{\boldsymbol\nu}=\pm\hat{\mathbf m}`` the prefactor is unity:

```math
\hat{\mathbf N}^{\mathcal R}
= \sqrt{\tfrac{b}{\pi}}\,
  \mathbf B^{\mathcal R}(\hat{\mathbf m},\hat{\mathbf n})
  \cdot\boldsymbol\Sigma\cdot\hat{\mathbf n},
\qquad
\hat{\mathbf K}^{\mathcal R}
= \sqrt{\pi b}\,\boldsymbol\Sigma\cdot\hat{\mathbf n}.
```

The SIF of an infinite ribbon crack is **independent of the matrix
stiffness**.

### Mode decomposition

The classical ``(K_{I},K_{II},K_{III})`` decomposition on
``(\hat{\mathbf n},\hat{\boldsymbol\nu},\hat{\boldsymbol\tau})`` reads

```math
K_{I}   = |\hat{\mathbf K}\cdot\hat{\mathbf n}|,\quad
K_{II}  = |\hat{\mathbf K}\cdot\hat{\boldsymbol\nu}|,\quad
K_{III} = |\hat{\mathbf K}\cdot\hat{\boldsymbol\tau}|.
```

Evaluation is handled by [`sif`](@ref) and [`dif`](@ref); see
[`src/Cracks/sif.jl`](https://github.com/MicMacTools/MeanFieldHom.jl/tree/main/src/Cracks/sif.jl).

## Dispatch and implementation notes

| `(crack, C₀)`                                   | `:auto` selects        | alternative(s)    | ForwardDiff |
| :---------------------------------------------- | :--------------------- | :---------------- | :---------: |
| `EllipticCrack / RibbonCrack, TensISO`          | `Analytical`           | —                 |     ✓       |
| `EllipticCrack / RibbonCrack, TensTI` (aligned) | `Analytical` (MFH)     | —                 |     ✓       |
| `EllipticCrack / RibbonCrack, AbstractTens{4,3}`| `Residue` (Float64)    | `:decuhr`         |  ✓ (decuhr) |

Entry points: [`cod_tensor`](@ref) / [`B_tensor`](@ref) for ``\mathbf B``,
[`compliance_contribution`](@ref) for ``\mathbb H``,
[`delta_compliance`](@ref) for ``\Delta\mathbb S``,
[`sif`](@ref) / [`dif`](@ref) for the front quantities.
