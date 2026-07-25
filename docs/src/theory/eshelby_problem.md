# The Eshelby inclusion problem

Everything in this section rests on one result. This page states it, and defines
the three tensors it produces — ``\mathbb{P}``, ``\mathbb{S}``, ``\mathbb{Q}``.
Their closed forms are the subject of the next page,
[Hill polarization tensors](hill_tensors.md).

## The problem

Fill ``\mathbb{R}^3`` with a homogeneous linear elastic medium of stiffness
``\mathbb{C}``, and single out an ellipsoid ``\mathcal{E}_{\boldsymbol{A}}``
centered at the origin (shape tensor ``\boldsymbol{A}``, semi-axes
``a\ge b\ge c``; see [Notation](notation.md#Ellipsoid-geometry)). Prescribe a
uniform **polarization stress** ``\boldsymbol{\tau}`` inside the ellipsoid and
zero outside, with no remote loading.

[eshelby1957](@cite) showed that the resulting strain field is **uniform inside
the ellipsoid**. This is the whole reason mean-field homogenization works, and
it is specific to the ellipsoid: no other shape has it.

```math
\forall\,\underline{x}\in\mathcal{E}_{\boldsymbol{A}},\qquad
\boldsymbol{\varepsilon}(\underline{x}) = -\,\mathbb{P}:\boldsymbol{\tau}.
```

The order-4 tensor ``\mathbb{P} = \mathbb{P}(\boldsymbol{A},\mathbb{C})`` is the
**Hill polarization tensor** [hill1963](@cite), [willis1977](@cite). It depends
on nothing but the ellipsoid's shape and orientation, through
``\boldsymbol{A}``, and the reference medium, through ``\mathbb{C}``. In
particular it is **independent of the ellipsoid's size** — only the ratios of
the semi-axes matter.

## Three tensors, one contraction apart

The same solution is written three ways in the literature. Knowing which is
which avoids most confusion when comparing formulas across papers.

**Eshelby tensor ``\mathbb{S}``.** Introduce the equivalent *eigenstrain*
(stress-free strain) ``\boldsymbol{\varepsilon}^{\star} =
-\mathbb{C}^{-1}:\boldsymbol{\tau}``. Then

```math
\boldsymbol{\varepsilon}(\underline{x}) = \mathbb{S}:\boldsymbol{\varepsilon}^{\star},
\qquad
\boxed{\;\mathbb{S} = \mathbb{P}:\mathbb{C}\;}
```

which is Eshelby's original form. ``\mathbb{S}`` is dimensionless;
``\mathbb{P}`` has the dimension of a compliance.

**Second Hill tensor ``\mathbb{Q}``.** Asking for the *stress* inside the
inclusion rather than the strain gives the dual statement

```math
\boldsymbol{\sigma}(\underline{x}) = -\,\mathbb{Q}:\boldsymbol{\varepsilon}^{\star},
\qquad
\boxed{\;\mathbb{Q} = \mathbb{C} - \mathbb{C}:\mathbb{P}:\mathbb{C}\;}
```

``\mathbb{Q}`` is what degenerates in a controlled way when the inclusion
becomes flat, which is why the crack theory is built on it rather than on
``\mathbb{P}`` — see [Crack opening displacement](cod_tensors.md).

In `MeanFieldHom` these are [`hill_tensor`](@ref), [`eshelby_tensor`](@ref) and
`hill_dual`.

## The transport counterpart

Replace elasticity by a scalar diffusion problem — heat conduction, mass
diffusion, electric conduction, Darcy flow. The unknown is a scalar potential
``T``, the flux is ``\underline{q} = -\boldsymbol{K}\cdot\nabla T``, and the
reference property is an order-2 conductivity ``\boldsymbol{K}``. Eshelby's
uniformity result holds unchanged: a uniform polarization flux inside the
ellipsoid produces a **uniform gradient** inside it, and an order-2 Hill tensor
``\boldsymbol{P}(\boldsymbol{A},\boldsymbol{K})`` plays the role of
``\mathbb{P}`` [willis1977](@cite):

```math
\nabla T(\underline{x}) = -\,\boldsymbol{P}\cdot\underline{\tau}_q,
\qquad
\boldsymbol{s} = \boldsymbol{P}\cdot\boldsymbol{K},
\qquad
\boldsymbol{Q} = \boldsymbol{K} - \boldsymbol{K}\cdot\boldsymbol{P}\cdot\boldsymbol{K}.
```

The two problems are handled by the same functions in `MeanFieldHom`, which
dispatch on the order of the property tensor passed in: an order-4
``\mathbb{C}`` selects the elastic path, an order-2 ``\boldsymbol{K}`` the
transport one.

```julia
using MeanFieldHom, TensND

sphere = Ellipsoid(1.0)

C₀ = TensISO{3}(3 * 70.0, 2 * 30.0)   # order 4: elasticity
P₄ = hill_tensor(sphere, C₀)          # → order-4 tensor
S₄ = eshelby_tensor(sphere, C₀)       # = P₄ : C₀

K₀ = TensISO{3}(2.5)                  # order 2: conduction
P₂ = hill_tensor(sphere, K₀)          # → order-2 tensor
s₂ = eshelby_tensor(sphere, K₀)       # = P₂ ⋅ K₀
```

## Why this matters for a real material

A real heterogeneous material is not one ellipsoid in an infinite medium. The
step from this idealized problem to an estimate of effective properties is the
subject of the next two pages, and it has two parts:

1. each inclusion is treated as if it were alone in an infinite *reference*
   medium — this is what makes ``\mathbb{P}`` usable, and it is exactly the
   approximation that distinguishes one mean-field scheme from another
   ([Homogenization schemes](homogenization.md));
2. the choice of that reference medium is the scheme: the matrix itself
   (Mori–Tanaka), the effective medium being sought (self-consistent), or
   something in between ([Localization](localization.md)).

!!! note "The uniformity result is exact"
    Nothing is approximated on *this* page: for a single ellipsoidal inclusion
    in an infinite homogeneous medium, ``\mathbb{P}`` is exact. All the
    approximation in mean-field homogenization enters later, when many
    interacting inclusions are replaced by many non-interacting ones.
