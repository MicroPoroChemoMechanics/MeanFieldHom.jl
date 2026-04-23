# Theoretical overview

`MeanFieldHom` implements the Eshelby / Hill polarisation machinery of
the [Echoes manual](https://github.com/jeanfrancoisbarthelemy/echoes)
in pure Julia, with matching conventions, expressions and
bibliography.

## Hill polarisation tensors

For an ellipsoid ``\mathcal E_{\mathbf A}`` embedded in a reference
medium of stiffness ``\mathbb C``, the 4th-order Hill polarisation
tensor is [willis1977](@cite), [mura1987](@cite):

```math
\mathbb P(\mathbf A,\mathbb C)
= \frac{\det\mathbf A}{4\pi}\!\int_{\|\hat{\boldsymbol\xi}\|=1}
\frac{\hat{\boldsymbol\xi}\stackrel{s}{\otimes}
      \bigl(\hat{\boldsymbol\xi}\cdot\mathbb C\cdot\hat{\boldsymbol\xi}\bigr)^{-1}
      \stackrel{s}{\otimes}\hat{\boldsymbol\xi}}
     {\|\mathbf A\cdot\hat{\boldsymbol\xi}\|^{3}}\,\mathrm dS_{\xi}.
```

For a 2nd-order conductivity ``\mathbf K`` the counterpart is obtained
by replacing ``(\hat{\boldsymbol\xi}\cdot\mathbb C\cdot\hat{\boldsymbol\xi})^{-1}``
with ``\hat{\boldsymbol\xi}\otimes\hat{\boldsymbol\xi}/
(\hat{\boldsymbol\xi}\cdot\mathbf K\cdot\hat{\boldsymbol\xi})``.

The Eshelby tensor follows from a single contraction:
``\mathbb S = \mathbb P:\mathbb C`` (order 4) or
``\mathbf s = \mathbf P\cdot\mathbf K`` (order 2). See the theory page
[Hill polarisation tensors](hill_tensors.md) for the closed forms and
the dispatch table.

## Flat cracks

A flat crack is the ``\omega\to 0`` limit of a flat spheroid; the
**crack compliance tensor** and the size-independent
**crack-opening-displacement (COD) tensor** ``\mathbf B`` are linked
by [kachanov1992](@cite), [barthelemyIJES2021](@cite):

```math
\mathbb H = \lim_{\omega\to 0}\,\omega\,\mathbb Q^{-1}
         = \tfrac{3}{4}\,\hat{\mathbf n}\stackrel{s}{\otimes}\mathbf B
                               \stackrel{s}{\otimes}\hat{\mathbf n},
```

with ``\mathbb Q = \mathbb C - \mathbb C:\mathbb P:\mathbb C`` the
second Hill tensor. Stress and displacement intensity factors at the
front are purely local and are exchanged by the COD tensor of the
tangent ribbon crack [kanaun1981](@cite), [kunin1983](@cite):
``\hat{\mathbf K} = \pi\,(\mathbf B^{\mathcal R})^{-1}\cdot\hat{\mathbf N}``.
Details in [COD tensors](cod_tensors.md).

## Elliptic integrals

Triaxial Newton-potential coefficients involve complete elliptic
integrals of first and second kind [abramowitz1972](@cite),
evaluated through Carlson symmetric forms [carlson1995](@cite). See
the dedicated theory page [Elliptic integrals](elliptic_integrals.md).
