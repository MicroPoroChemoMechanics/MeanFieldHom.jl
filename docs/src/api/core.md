# [API — Core](@id api-core)

```@docs
MeanFieldHom.Core
MeanFieldHom.AbstractInclusion
MeanFieldHom.AbstractEllipsoidalInclusion
MeanFieldHom.AbstractCrack
MeanFieldHom.AbstractLayeredInclusion
MeanFieldHom.AbstractAlgorithm
MeanFieldHom.Analytical
MeanFieldHom.Residue
MeanFieldHom.DECUHR
MeanFieldHom.CylinderQuadrature
MeanFieldHom.MaterialSymmetry
MeanFieldHom.material_symmetry
MeanFieldHom.Core.extract_iso_moduli
MeanFieldHom.Core.extract_ti_moduli
MeanFieldHom.Core.newton_potential_3d
MeanFieldHom.Core.newton_potential_2d
MeanFieldHom.Core.newton_potential_3d_cylinder
MeanFieldHom.Core.dimension
MeanFieldHom.Core.element_type
MeanFieldHom.Core.inclusion_basis
MeanFieldHom.Core.shape_trait
MeanFieldHom.Core.shape_tensor
```

## Exact rotation-group averages

Exact averages of a tensor over a rotation group — the *exact* counterpart of the
best-fit projections in [API — Schemes](schemes.md). See
[Symmetrization showcase](../tutorials/generated/symmetrization.md) for
the comparison between the two.

```@docs
MeanFieldHom.Core.isotropify
MeanFieldHom.Core.transverse_isotropify
```
