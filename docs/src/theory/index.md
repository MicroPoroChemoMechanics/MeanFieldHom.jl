# [Theory — reading path](@id th-index)

`MeanFieldHom` computes effective properties of heterogeneous materials by
mean-field homogenization. This section states the theory it implements, in the
order in which it is built. Every page is self-contained on notation
([Notation and conventions](notation.md)) and every formula is either cited or
derived on the page.

## The chain, in one paragraph

An **inclusion** embedded in an infinite reference medium responds to a remote
load in a way entirely captured by one object, the **Hill polarization tensor**
``\mathbb{P}`` — this is [Eshelby's result](eshelby_problem.md), and
``\mathbb{P}`` depends only on the inclusion *shape* and the reference *moduli*
([Hill polarization tensors](hill_tensors.md)). From ``\mathbb{P}`` follows the
**localization tensor**, which says how much of the remote load each phase
actually sees, and hence each phase's **contribution** to the effective
stiffness ([Localization](localization.md)). Assembling those contributions
under an assumption about how phases interact gives a **homogenization scheme**
— dilute, Mori–Tanaka, self-consistent, differential, and the bounds
([Homogenization schemes](homogenization.md)).

Everything else in this section is a specialization of that chain:

| Page | Specialization |
| :--- | :------------- |
| [The finite Eshelby cell](corrected_cell.md) | when no closed form exists: the inclusion is solved on a *finite* cell and the truncation bias removed by its own dipole far field |
| [Differential scheme](differential_scheme.md) | incorporation as an ODE in a fictitious time, and what a crack — which has no volume to replace — does to it |
| [Crack opening displacement](cod_tensors.md) | the flat-inclusion limit: a crack has no volume, so it is described by ``\boldsymbol{B}`` and ``\mathbb{H}`` instead of a volume fraction |
| [Thermal cracks](thermal_cracks.md) | the same limit for scalar transport |
| [Layered sphere](layered_sphere.md) | a *composite* inclusion: no Hill tensor exists, the response is assembled by a radial recurrence |
| [Layered spheroid](layered_spheroid.md) | the same for confocal spheroids, where imperfect interfaces couple harmonic degrees |
| [Ageing linear viscoelasticity](viscoelasticity.md) | moduli become Volterra operators; the algebra of the chain is unchanged |
| [Elliptic integrals](elliptic_integrals.md) | mathematical appendix: the special functions the closed forms need |

## Where the two physics meet

The Eshelby framework applies verbatim to elasticity (order-4 tensors) and to
scalar transport — heat conduction, diffusion, electric conduction, Darcy flow
(order-2 tensors). The documentation treats them in parallel rather than in
separate sections, because the algebra is identical:

| | Elasticity | Transport |
| :--- | :--- | :--- |
| property | stiffness ``\mathbb{C}`` | conductivity ``\boldsymbol{K}`` |
| Hill tensor | ``\mathbb{P}(\boldsymbol{A},\mathbb{C})`` (order 4) | ``\boldsymbol{P}(\boldsymbol{A},\boldsymbol{K})`` (order 2) |
| Eshelby tensor | ``\mathbb{S} = \mathbb{P}:\mathbb{C}`` | ``\boldsymbol{s} = \boldsymbol{P}\cdot\boldsymbol{K}`` |
| crack descriptor | COD tensor ``\boldsymbol{B}``, compliance ``\mathbb{H}`` | COD scalar ``b``, resistivity ``\boldsymbol{R}`` |

One asymmetry is worth knowing in advance: for an **arbitrarily anisotropic**
matrix the order-2 Hill tensor has a *closed form* (via a square-root
transformation), whereas the order-4 one requires numerical cubature. This is
why the conductivity paths in `MeanFieldHom` are analytical far more often than
the elastic ones.

## Relation to the Echoes manual

`MeanFieldHom` is a pure-Julia reimplementation of the Eshelby/Hill machinery of
the [Echoes manual](https://github.com/jeanfrancoisbarthelemy/echoes), and the
Hill-tensor pages follow its appendix closely — same expressions, same
conventions, same bibliography. Two families of difference are flagged
explicitly wherever they occur:

- **extensions**: infinite cylinders and 2-D plane strain as first-class
  inclusion types, an analytical transversely isotropic path, automatic
  differentiation and symbolic number types throughout;
- **conventions that genuinely differ**: the crack opening displacement tensor
  ``\boldsymbol{B}`` and the compliance ``\mathbb{H}`` are the notable case.
  `MeanFieldHom` computes ``\boldsymbol{B}`` first and derives ``\mathbb{H}``
  from it, where Echoes computes ``\mathbb{H}`` directly and never forms
  ``\boldsymbol{B}``. The competing normalizations are named and compared in
  [Crack opening displacement](cod_tensors.md).
