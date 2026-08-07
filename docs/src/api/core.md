# [API — Core](@id api-core)

```@docs
MeanFieldHomogenization.Core
MeanFieldHomogenization.AbstractInclusion
MeanFieldHomogenization.AbstractEllipsoidalInclusion
MeanFieldHomogenization.AbstractCrack
MeanFieldHomogenization.AbstractLayeredInclusion
MeanFieldHomogenization.AbstractAlgorithm
MeanFieldHomogenization.Analytical
MeanFieldHomogenization.Residue
MeanFieldHomogenization.DECUHR
MeanFieldHomogenization.CylinderQuadrature
MeanFieldHomogenization.MaterialSymmetry
MeanFieldHomogenization.material_symmetry
MeanFieldHomogenization.Core.extract_iso_moduli
MeanFieldHomogenization.Core.extract_ti_moduli
MeanFieldHomogenization.Core.newton_potential_3d
MeanFieldHomogenization.Core.newton_potential_2d
MeanFieldHomogenization.Core.newton_potential_3d_cylinder
MeanFieldHomogenization.Core.dimension
MeanFieldHomogenization.Core.element_type
MeanFieldHomogenization.Core.inclusion_basis
MeanFieldHomogenization.Core.shape_trait
MeanFieldHomogenization.Core.shape_tensor
```

## Exact rotation-group averages

Exact averages of a tensor over a rotation group — the *exact* counterpart of the
best-fit projections in [API — Schemes](schemes.md). See
[Symmetrization showcase](../tutorials/generated/symmetrization.md) for
the comparison between the two.

```@docs
MeanFieldHomogenization.Core.isotropify
MeanFieldHomogenization.Core.transverse_isotropify
```
