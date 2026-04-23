# Thermal cracks — COD scalar and resistivity contribution

This page gives the analogue of the **elasticity crack quantities** for
the 2nd-order (conductivity / diffusion / Darcy) problem.  It is the
direct transposition of the theory in
[Crack opening displacement and compliance tensors](cod_tensors.md)
to the scalar-potential problem, where the driving field is a
vector (heat flux or gradient) rather than a symmetric 2-tensor.

## Elasticity ↔ Conductivity — correspondence table

| Elasticity (4-tensor problem)                                    | Conductivity (2-tensor problem)                                         |
| ---------------------------------------------------------------- | ----------------------------------------------------------------------- |
| Displacement ``\mathbf u`` — 1-tensor                            | Temperature ``T`` — scalar                                              |
| Stress ``\boldsymbol\sigma = \mathbb C:\boldsymbol\varepsilon``  | Heat flux ``\mathbf q = -\mathbf K_0\cdot\nabla T``                     |
| Stiffness ``\mathbb C`` — 4-tensor, 21 independent components    | Conductivity ``\mathbf K_0`` — 2-tensor, 6 independent components       |
| Hill tensor ``\mathbb P`` — 4-tensor                             | Hill tensor ``\mathbf P`` — 2-tensor                                    |
| COD tensor ``\mathbf B`` — **2-tensor** (6 components)           | COD scalar ``b`` — **scalar** (1 component)                             |
| Kachanov factorisation ``\mathbb H = k\,\hat{\mathbf n}\stackrel{s}{\otimes}\mathbf B\stackrel{s}{\otimes}\hat{\mathbf n}``, ``k=3/4`` (elliptic), ``k=2/\pi`` (ribbon) | Rank-1 factorisation ``\mathbf R = k\,b\,\hat{\mathbf w}\otimes\hat{\mathbf w}``, same ``k``, ``\hat{\mathbf w}\parallel\mathbf K_0^{-1/2}\hat{\mathbf n}`` |
| Dilute contribution ``\Delta\mathbb S = (4\pi/3)\varepsilon^{3\mathrm d}\mathbb H`` (elliptic), ``= \pi\varepsilon^{2\mathrm d}\mathbb H`` (ribbon) | Dilute contribution ``\Delta\mathbf R = (4\pi/3)\varepsilon^{3\mathrm d}\mathbf R`` (elliptic), ``= \pi\varepsilon^{2\mathrm d}\mathbf R`` (ribbon) |
| Sextic acoustic polynomial (Masson 2008)                         | Quadratic acoustic form ``\xi\cdot\mathbf K_0\cdot\xi`` → **analytical** |
| Stress intensity factors ``K_I, K_{II}, K_{III}``                | Heat-flux intensity factor ``K_T`` — scalar (mode I analogue only)      |
| Displacement intensity factor ``\hat{\mathbf N}``                | Temperature intensity factor — scalar ``[T]_\text{avg}``                |

Why a scalar ``b`` and not a tensor? In the 2nd-order problem, the
driving field ``\nabla T`` is a vector and the jump ``[T]`` across the
crack is a scalar; only the **normal component** of the heat flux
``\mathbf q\cdot\hat{\mathbf n}`` produces a non-trivial jump.  A
single scalar ``b`` captures the full crack flexibility — as opposed to
the 6-component ``\mathbf B`` of the 4-tensor problem, which has to
resolve sliding and shear modes.

Why a rank-1 direction ``\hat{\mathbf w} \ne \hat{\mathbf n}`` in
general?  Because the null space of the block ``\mathbf K_0 - \mathbf
K_0\mathbf P(0)\mathbf K_0`` — which emerges in the ``\omega\to 0``
limit of the Hill tensor — is aligned with
``\mathbf K_0^{-1/2}\hat{\mathbf n}``, **not** with ``\hat{\mathbf n}``.
The two coincide whenever ``\hat{\mathbf n}`` is an eigenvector of
``\mathbf K_0`` (isotropic matrix, TI aligned, orthotropic with
``\hat{\mathbf n}`` along a symmetry axis).

## Crack geometry

Same geometric families as in the elasticity chapter:

- **Elliptic cracks** of aspect ratio ``\eta = b/a \in (0, 1]``
  ([`EllipticCrack`](@ref)), including the circular penny ``\eta = 1``
  ([`PennyCrack`](@ref)).
- **Ribbon cracks** (tunnel cracks) of half-width ``b``
  ([`RibbonCrack`](@ref)), infinite along ``\hat{\boldsymbol\ell}``.

## Hill tensor Taylor expansion and block-matrix limit

Mirror of the elasticity derivation of
[Barthélémy (2009)](@cite barthelemyIJSS2009).  For a flat ellipsoidal
inclusion of aspect ratio ``\omega\to 0``, the 2nd-order Hill tensor
admits the expansion

```math
\mathbf P(\omega) = \mathbf P(0) + \omega\,\boldsymbol\Pi + o(\omega),
\qquad
\mathbf P(0) = \bigl(\mathbf K_0^{-1/2}\hat{\mathbf n}\bigr)
               \otimes
               \bigl(\mathbf K_0^{-1/2}\hat{\mathbf n}\bigr)
```

which is rank-1.  The acoustic block is then

```math
\boldsymbol\Lambda(\omega) = \mathbf K_0 - \mathbf K_0\mathbf P(\omega)\mathbf K_0
                           = \boldsymbol\Lambda(0) + \omega\,\boldsymbol\Lambda_1 + o(\omega),
```

with ``\boldsymbol\Lambda(0) = \mathbf K_0 - (\mathbf K_0^{1/2}\hat{\mathbf n})
(\mathbf K_0^{1/2}\hat{\mathbf n})^T``.  Its null space is spanned by
``\mathbf v = \mathbf K_0^{-1/2}\hat{\mathbf n}`` (one-dimensional in
the 2-tensor case — to be contrasted with the 3-dimensional null space
in the elasticity problem).  The limit

```math
\mathbf R
= \lim_{\omega\to 0}\omega\,\boldsymbol\Lambda(\omega)^{-1}
= \frac{1}{Y_{22}}\,\hat{\mathbf w}\otimes\hat{\mathbf w},
\qquad
\hat{\mathbf w} = \frac{\mathbf K_0^{-1/2}\hat{\mathbf n}}
                        {\sqrt{\hat{\mathbf n}\cdot\mathbf K_0^{-1}\hat{\mathbf n}}}
```

is rank-1, with ``Y_{22} = \mathbf v\cdot\boldsymbol\Lambda_1\cdot\mathbf v /
\|\mathbf v\|^2``.

## Closed-form COD scalar ``b``

### Isotropic matrix

For ``\mathbf K_0 = k_0\,\mathbf 1`` the square-root transform is
trivial and ``\hat{\mathbf w} = \hat{\mathbf n}``:

```math
\boxed{\;
b_{\text{ell}}^{\text{iso}}
= \frac{\eta}{\pi\,k_0\,\mathcal E_\eta}
\;},
\qquad
\mathcal E_\eta = \mathcal E\!\bigl(\sqrt{1-\eta^{2}}\bigr),
```

with ``\mathcal E`` the complete elliptic integral of the second kind.
Penny limit ``\eta = 1``: ``b = 2/(\pi^{2} k_0)``.  Ribbon:
``b = 2/(\pi k_0)`` (direct 2-D computation, not the ``\eta\to 0``
limit of the elliptic formula).

### Anisotropic matrix (K⁻¹ᐟ² transform)

Applying the square-root change of variable
``\tilde{\mathbf x} = \mathbf K_0^{-1/2}\mathbf x`` reduces the problem
to the isotropic one with a **transformed crack shape**.  Let
``\tilde{\mathbf A} = \mathbf A\cdot\mathbf K_0^{-1/2}`` with
``\mathbf A = \mathbf R_c\,\mathrm{diag}(a,b,0)\,\mathbf R_c^{T}``, and
denote its singular values ``\sigma_1 \ge \sigma_2 \ge \sigma_3 = 0``.
The transformed in-plane aspect ratio is
``\eta_t = \sigma_2/\sigma_1 \in (0,1]`` and the transformed ellipse
first-kind integral parameter is ``k_t^{2} = 1 - \eta_t^{2}``.  Then

```math
\boxed{\;
b_{\text{ell}}^{\text{aniso}}
= \frac{\sigma_2\,(\hat{\mathbf n}\cdot\mathbf K_0^{-1}\hat{\mathbf n})
         \,\sqrt{\hat{\mathbf n}\cdot\mathbf K_0\hat{\mathbf n}}}
       {\pi\,a_\text{max}\,\mathcal E_{\eta_t}}
\;},
\qquad
a_\text{max} = \max(a,b).
```

It reduces to the isotropic formula above when
``\mathbf K_0 = k_0\,\mathbf 1`` (then ``\sigma_2 = \min(a,b)/\sqrt{k_0}``
and ``\hat{\mathbf n}\cdot\mathbf K_0\hat{\mathbf n} = k_0``,
``\hat{\mathbf n}\cdot\mathbf K_0^{-1}\hat{\mathbf n} = 1/k_0``).
For a TI matrix aligned with ``\hat{\mathbf n}``
(``\mathbf K_0 = \mathrm{diag}(k_t, k_t, k_n)`` in the crack frame, penny):

```math
b_{\text{penny}}^{\text{aligned TI}}
= \frac{2}{\pi^{2}\,\sqrt{k_t k_n}}
\qquad\text{(geometric mean of }k_t\text{ and }k_n\text{)}.
```

### Ribbon crack — 2D formula

Only the ``(\hat{\mathbf m}, \hat{\mathbf n})`` transverse block of
``\mathbf K_0`` enters the formula:

```math
\boxed{\;
b_{\text{ribbon}}^{\text{aniso}}
= \frac{2}{\pi\,\sqrt{\det\bigl(\mathbf K_0|_{(\hat{\mathbf m},\hat{\mathbf n})}\bigr)}}
\;}
```

which reduces to ``b = 2/(\pi k_0)`` for an isotropic matrix.

## Resistivity contribution ``\mathbf R`` and dilute correction ``\Delta\mathbf R``

The size-independent **crack resistivity contribution tensor** is
assembled from the scalar ``b`` and the effective direction
``\hat{\mathbf w}``:

```math
\mathbf R^{\mathcal E} = \tfrac{3}{4}\,b\,\hat{\mathbf w}\otimes\hat{\mathbf w}
\qquad\text{(elliptic)},
\qquad
\mathbf R^{\mathcal R} = \tfrac{2}{\pi}\,b\,\hat{\mathbf w}\otimes\hat{\mathbf w}
\qquad\text{(ribbon)}.
```

The geometric prefactors ``3/4`` and ``2/\pi`` are the same as in the
elasticity case (they come from ``cS/V`` evaluated on the ellipsoidal
and ribbon geometries — see
[Crack compliance and COD tensor](cod_tensors.md#Crack-compliance-H-and-COD-tensor-B)).
In iso and aligned-TI cases ``\hat{\mathbf w} = \hat{\mathbf n}``.

[`compliance_contribution`](@ref)`(crack, K₀)` returns ``\mathbf R``
directly.  The dilute resistivity correction to the effective
resistivity of the cracked conductor is obtained via
[`delta_resistivity`](@ref)`(crack, R, ε)`:

```math
\Delta\mathbf R = \tfrac{4\pi}{3}\,\varepsilon^{3\mathrm d}\,\mathbf R
\qquad\text{(elliptic, }\varepsilon^{3\mathrm d} = Nab^{2}\text{)},
\qquad
\Delta\mathbf R = \pi\,\varepsilon^{2\mathrm d}\,\mathbf R
\qquad\text{(ribbon, }\varepsilon^{2\mathrm d} = Nb^{2}\text{)}.
```

These reduce to the Sevostianov–Kachanov expressions
(see [Sevostianov & Kachanov (2002)](@cite sevostianov2002),
 [Kachanov (2018)](@cite kachanov2018)).

## Intensity factors

Thermal analogues of the elastic stress / displacement intensity
factors:

- **Heat-flux intensity factor** ``K_T`` (mode I analogue, scalar):
  the singular crack-tip field scales as ``\sim K_T/\sqrt{r}``.  For a
  ribbon crack of half-width ``b`` and remote flux ``\mathbf q^\infty``,
  ``K_T = \sqrt{\pi b}\,(\hat{\mathbf n}\cdot\mathbf q^{\infty})``.
  For an elliptic crack the formula involves the tangent-ribbon COD
  ratio exactly as in the elasticity case (``b^\mathcal E/b^\mathcal R``
  replaces ``\mathbf B^\mathcal E(\mathbf B^\mathcal R)^{-1}``).
- **Temperature intensity factor** (scalar):
  ``[T]_\text{avg} = b\,(\hat{\mathbf n}\cdot\mathbf q^{\infty})``.

See [`sif`](@ref) and [`dif`](@ref) for the full signatures (dispatched
on ``\mathbf K_0::\texttt{AbstractTens\{2,3\}}``).
