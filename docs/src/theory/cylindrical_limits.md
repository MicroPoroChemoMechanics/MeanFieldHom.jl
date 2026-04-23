# Cylindrical limits of the Hill tensor

`MeanFieldHom.jl` treats an infinite cylinder as a first-class
inclusion type (`Cylinder{S, T, B}`), obtained by passing `a → ∞` in
the ellipsoid family. This page documents how the Newton-potential
coefficients and the Hill polarisation tensor specialise in that
limit, matching the cylinder column of the Echoes appendix
([Mura 1987](@cite mura1987), §11.22).

## Geometry and local frame

Let the cylinder axis be ``\hat{\mathbf e}_1`` in the inclusion's
local basis and let ``(b, c)`` with ``b \ge c > 0`` be the two
transverse semi-axes (along ``\hat{\mathbf e}_2`` and
``\hat{\mathbf e}_3``). The family contains two shape traits:

- `CircularCylindrical` when ``b = c`` — transversely isotropic
  response, returned as a `TensWalpole` with axis ``\hat{\mathbf e}_1``;
- `EllipticCylindrical` when ``b > c`` — orthotropic response,
  returned as a `TensOrtho`.

## Newton-potential coefficients in the cylinder limit

Passing to the limit ``a\to\infty`` in the triaxial ``I_i^{\mathbf A}``
/ ``I_{ij}^{\mathbf A}`` formulas yields the Echoes-appendix cylinder
column (normalised convention ``\sum_i I_i^{\mathbf A} = 1``):

```math
I_1^{\text{cyl}} = 0,\qquad
I_2^{\text{cyl}} = \frac{c}{b+c},\qquad
I_3^{\text{cyl}} = \frac{b}{b+c},
```

```math
I_{22}^{\text{cyl}} = \frac{c(2b+c)}{3\,b^{2}(b+c)^{2}},\qquad
I_{33}^{\text{cyl}} = \frac{b(b+2c)}{3\,c^{2}(b+c)^{2}},\qquad
I_{23}^{\text{cyl}} = \frac{1}{(b+c)^{2}},
```

```math
I_{11}^{\text{cyl}} = I_{12}^{\text{cyl}} = I_{13}^{\text{cyl}} = 0,
\qquad
a^{2}\,I_{12}^{\text{cyl}} \xrightarrow[a\to\infty]{} I_2^{\text{cyl}},
\qquad
a^{2}\,I_{13}^{\text{cyl}} \xrightarrow[a\to\infty]{} I_3^{\text{cyl}}.
```

The products ``a^{2}\,I_{12}`` and ``a^{2}\,I_{13}`` remain finite in
the limit and are the key to the consistency identity
``3\rho_{i}^{2}I_{ii}^{\mathbf A}+\sum_{j\ne i}\rho_{j}^{2}I_{ij}^{\mathbf A}
=3I_i^{\mathbf A}`` at the cylinder endpoint.

The circular sub-case ``b = c`` is evaluated on a dedicated branch to
bypass any ``(b^{2}-c^{2})^{-1}``-style intermediate — the expressions
above are regular but the triaxial formula they descend from is not,
hence the dispatch split via the `CylindricalShape` trait.

!!! note "Storage convention"
    Internally `newton_potential_3d_cylinder` returns the **raw**
    (un-normalised) kernel, i.e. the values above multiplied by
    ``4\pi``; the normalising division is applied at the
    [`tens_IA`](@ref) call site.  The expressions on this page match
    the Echoes manual normalisation (``\sum_i I_i^{\text{cyl}}=1``).

## Consequences for the Hill tensor

Two independent calculation paths converge on the same Hill tensor.

**Isotropic matrix.** Substituting the cylinder coefficients above in
the Kelvin–Mandel forms of ``\mathbb U^{\mathbf A}`` and
``\mathbb V^{\mathbf A}`` yields the Eshelby / Hill closed form of
[Mura 1987](@cite mura1987), §11.22. The first row and column vanish,
``P^{\text{cyl}}_{1jkl}\equiv 0``, reflecting the fact that no
polarisation is transmitted along the cylinder axis. See
[hill_tensors.md](hill_tensors.md) for the explicit matrices.

**Anisotropic matrix.** For a general 3D anisotropic stiffness
``\mathbb C``, the Willis integral reduces to a 1-D quadrature over
the transverse plane:

```math
\mathbb P^{\text{cyl}}
= \frac{b\,c}{2\pi}\!\int_{0}^{2\pi}
\frac{\hat{\boldsymbol\xi}\stackrel{s}{\otimes}
      \bigl(\hat{\boldsymbol\xi}\cdot\mathbb C\cdot\hat{\boldsymbol\xi}\bigr)^{-1}
      \stackrel{s}{\otimes}\hat{\boldsymbol\xi}}
     {\|\mathbf A\cdot\hat{\boldsymbol\xi}\|^{3}}\,\mathrm d\varphi,
\qquad
\hat{\boldsymbol\xi} = \bigl(0,\,\cos\varphi,\,\sin\varphi\bigr).
```

The axial component of ``\hat{\boldsymbol\xi}`` is identically zero,
so the double surface cubature of the general 3-D `DECUHR` path
collapses to a single-variable `QuadGK` integral that remains
ForwardDiff-compatible.

The in-plane components of this tensor coincide exactly with the
solution of the 2D plane-strain problem produced by the 2-D
anisotropic Hill routine — the cylinder geometry is the 3-D
realisation of the 2-D ellipse.

## Residue algorithm unavailable

The residue representation of [Masson 2008](@cite masson2008) relies
on the six complex roots of the acoustic polynomial along
``\hat{\boldsymbol\xi}_3``. For an infinite cylinder one root escapes
to infinity and the polynomial degenerates — the algorithm is
therefore **not applicable**. Calling
`hill_tensor(Cylinder(…), C₀; method=:residue)` silently falls back
to the dedicated 1-D quadrature via `CylinderQuadrature`, preserving
user ergonomics.

## Cross-references

- `Cylinder`, `CircularCylindrical`, `EllipticCylindrical` in the
  [Elasticity API](../api/elasticity.md).
- `newton_potential_3d_cylinder` in the
  [Core API](../api/core.md).
- Manual page: [cylindrical inclusions](../manual/cylindrical_inclusions.md).
