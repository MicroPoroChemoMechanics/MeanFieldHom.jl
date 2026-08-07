# [Notation and conventions](@id th-notation)

This page fixes the notation used throughout the documentation. It is
deliberately short: every symbol that is *not* listed here is redefined on the
page where it appears.

## Tensor order is carried by the typeface

The convention is that of the
[Echoes manual](https://jfbarthelemy.github.io/echoes/), so that
formulas can be compared side by side with it.

| Object | Typeset as | Example |
| :----- | :--------- | :------ |
| scalar | italic | ``k``, ``\mu``, ``\nu``, ``\omega``, ``\eta``, ``\chi`` |
| vector (order 1) | underlined | ``\underline{u}``, ``\underline{n}``, ``\underline{\xi}`` |
| tensor of order 2 | bold | ``\boldsymbol{A}``, ``\boldsymbol{B}``, ``\boldsymbol{K}``, ``\boldsymbol{\sigma}``, ``\boldsymbol{\varepsilon}`` |
| tensor of order 4 | blackboard bold | ``\mathbb{C}``, ``\mathbb{P}``, ``\mathbb{Q}``, ``\mathbb{H}``, ``\mathbb{S}`` |
| set, geometry, shape function | calligraphic | ``\mathcal{E}_{\boldsymbol{A}}``, ``\mathcal{G}_i``, ``\mathcal{T}_a``, ``\mathcal{K}_\eta`` |

The single most useful thing to remember: **an underline is a vector, bold is
order 2, blackboard bold is order 4**. This is what makes an expression such as

```math
\mathbb{H} = \tfrac{3}{4}\,
\underline{n}\stackrel{s}{\otimes}\boldsymbol{B}\stackrel{s}{\otimes}\underline{n}
```

readable at a glance: an order-4 tensor ``\mathbb{H}`` is built from an order-2
tensor ``\boldsymbol{B}`` and a unit vector ``\underline{n}``.

## Operators

| Operator | Meaning |
| :------- | :------ |
| ``\underline{u}\cdot\underline{v}``, ``\boldsymbol{A}\cdot\underline{u}`` | simple contraction (one index) |
| ``\boldsymbol{A}:\boldsymbol{B}``, ``\mathbb{C}:\boldsymbol{\varepsilon}`` | double contraction (two indices) |
| ``\underline{u}\otimes\underline{v}`` | tensor product |
| ``\underline{u}\stackrel{s}{\otimes}\boldsymbol{A}\stackrel{s}{\otimes}\underline{u}`` | symmetrized tensor product (order 4, both symmetries) |
| ``\boldsymbol{A}\stackrel{s}{\boxtimes}\boldsymbol{B}`` | symmetrized box product |
| ``\boldsymbol{A}^{\!T}`` | transpose |
| ``[\![\,\underline{u}\,]\!]`` | jump across an interface, ``\underline{u}^+-\underline{u}^-`` |
| ``\langle\,\cdot\,\rangle_\Omega`` | volume average over ``\Omega`` |
| ``\mathrm{d}S_\xi`` | surface element in the ``\underline{\xi}`` parametrization |

## Isotropic and transversely isotropic bases

Any isotropic order-4 tensor is a combination of the spherical and deviatoric
projectors, with ``\mathbb{I}`` the order-4 identity:

```math
\mathbb{J} = \tfrac{1}{3}\,\boldsymbol{1}\otimes\boldsymbol{1},
\qquad
\mathbb{K} = \mathbb{I} - \mathbb{J},
\qquad
\mathbb{C}^{\mathrm{iso}} = 3k\,\mathbb{J} + 2\mu\,\mathbb{K}.
```

A transversely isotropic order-4 tensor needs five coefficients. The
documentation uses the **Walpole basis** ``(\mathbb{W}_1,\dots,\mathbb{W}_6)``,
in which a *major-symmetric* tensor has ``\mathbb{W}_3`` and ``\mathbb{W}_4``
sharing a single coefficient — hence five independent numbers written
``(P_1, P_2, P_3, P_5, P_6)``, with no ``P_4``:

```math
\mathbb{P} = P_1\,\mathbb{W}_1 + P_2\,\mathbb{W}_2
           + P_3\,(\mathbb{W}_3+\mathbb{W}_4)
           + P_5\,\mathbb{W}_5 + P_6\,\mathbb{W}_6 .
```

This is the storage used by `TensND.TensTI{4}` and by
[`hill_tensor`](@ref) on a transversely isotropic matrix
([walpole1981](@cite), [barthelemyIJES2020_hilltrans](@cite)).

## Ellipsoid geometry

An ellipsoid is described by an invertible order-2 **shape tensor**
``\boldsymbol{A}`` such that ``\boldsymbol{A}^{\!T}\!\cdot\boldsymbol{A}`` is
symmetric positive definite, with eigenvalues the squared semi-axes:

```math
\underline{x}\in\mathcal{E}_{\boldsymbol{A}}
\iff
\underline{x}\cdot(\boldsymbol{A}^{\!T}\!\cdot\boldsymbol{A})^{-1}
\cdot\underline{x}\le 1,
\qquad
\boldsymbol{A}^{\!T}\!\cdot\boldsymbol{A}
= \sum_{i=1}^{3}\rho_i^{2}\,
  \underline{e}^{\boldsymbol{A}}_i\otimes\underline{e}^{\boldsymbol{A}}_i .
```

Semi-axes are ordered ``a = \rho_1 \ge b = \rho_2 \ge c = \rho_3``, with
``(\underline{e}^{\boldsymbol{A}}_i)`` the orthonormal eigenvectors. Two aspect
ratios recur:

- ``\eta = b/a`` — in-plane aspect ratio (of a crack, in particular);
- ``\omega`` — the *flatness* ``c/a``, or the axial-to-transverse ratio of a
  spheroid, depending on the page. **It is always redefined locally**, because
  the two usages differ.

The figure below is drawn from an actual [`Ellipsoid`](@ref) instance, so the
guides are the stored ``\rho_i`` in the stored frame — the same ones every
closed form on the following pages is written in. Semi-axes
``(a, b, c) = (3, 1.5, 0.8)``, hence ``\eta = 0.5`` and ``\omega \approx 0.27``.

```@setup notation
using MeanFieldHomogenization
include(joinpath(pkgdir(MeanFieldHomogenization), "scripts", "common", "docviz.jl"))
```

```@example notation
plotly_scene(shape_traces(Ellipsoid(3.0, 1.5, 0.8)); uid = "notation-ellipsoid",
    height = 430, title = "Semi-axes a ≥ b ≥ c and the principal frame")
```

## Storage: Kelvin–Mandel

Order-2 and order-4 tensors are stored in the **Kelvin–Mandel** convention
(orthonormal 6-dimensional basis, off-diagonal components carrying
``\sqrt{2}`` and ``2`` factors), not in the engineering Voigt convention. This
makes a double contraction an ordinary matrix product and a tensor inverse an
ordinary matrix inverse. The practical consequences when reading printed
components are spelled out in
[A storage convention worth knowing](../tutorials/first_estimate.md).

``\mathrm{Mat}(\mathbb{C})`` denotes the ``6\times 6`` Kelvin–Mandel matrix of an
order-4 tensor, and
``\mathrm{Mat}\bigl(\mathbb{C}, \underline{e}^{\boldsymbol{A}}_i\bigr)`` the same
matrix expressed in the ellipsoid's principal frame.

## Two rules this documentation follows

**Every symbol is defined where it is used.** A page never relies on a symbol
introduced only on another page, even at the cost of repeating a definition.

!!! warning "No formula without a traceable source"
    Every expression in this documentation is either (i) accompanied by a
    citation to published work, or (ii) derived explicitly on the page from
    expressions that are. Where a convention differs between references — the
    crack opening displacement tensor ``\boldsymbol{B}`` is the notable case,
    see [Crack opening displacement](cod_tensors.md) — the competing
    conventions are named and the one implemented by `MeanFieldHomogenization` is stated.
