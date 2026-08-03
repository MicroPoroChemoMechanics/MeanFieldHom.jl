# [API — Laminate](@id api-laminate)

## The cell

```@docs
MeanFieldHom.Laminates
MeanFieldHom.Laminates.Laminate
MeanFieldHom.Laminates.Layer
MeanFieldHom.Laminates.add_layer!
MeanFieldHom.Laminates.layer_names
MeanFieldHom.Laminates.layer_property
MeanFieldHom.Laminates.layer_property_raw
MeanFieldHom.Laminates.layer_thickness
MeanFieldHom.Laminates.laminate_period
MeanFieldHom.Laminates.laminate_basis
MeanFieldHom.Laminates.laminate_normal
MeanFieldHom.Laminates.validate_laminate
```

## Anisotropic interfaces

A plane, unlike a sphere, imposes no symmetry on the interface, so a laminate
also accepts tensor-valued interface properties.

```@docs
MeanFieldHom.Laminates.AnisotropicSpringInterface
MeanFieldHom.Laminates.AnisotropicMembraneInterface
MeanFieldHom.Laminates.AnisotropicSurfaceConductiveInterface
```

## Fields and per-layer tensors

```@docs
MeanFieldHom.Laminates.laminate_hill
MeanFieldHom.Laminates.layer_strain_localization
MeanFieldHom.Laminates.layer_stress_localization
MeanFieldHom.Laminates.layer_gradient_localization
MeanFieldHom.Laminates.layer_flux_localization
MeanFieldHom.Laminates.interface_jump
```

## Parameter lenses

```@docs
MeanFieldHom.Laminates.ThicknessParameter
MeanFieldHom.Laminates.thickness
MeanFieldHom.Laminates.InterfaceParameter
MeanFieldHom.Laminates.interface_param
```

## Ageing viscoelasticity

```@docs
MeanFieldHom.Viscoelasticity.laminate_alv
```

## The block algebra

The kernel behind the cell, in `Core`, so that the ageing-viscoelastic
laminate reuses it with the Volterra inversion substituted.

```@docs
MeanFieldHom.Core.KM_IP
MeanFieldHom.Core.KM_OP
MeanFieldHom.Core.plane_pinv
MeanFieldHom.Core.plane_pinv2
MeanFieldHom.Core.flat_hill
MeanFieldHom.Core.acoustic_tensor
MeanFieldHom.Core.compliance_op_block
MeanFieldHom.Core.laminate_stiffness
MeanFieldHom.Core.laminate_conductivity
MeanFieldHom.Core.laminate_strain_localization
MeanFieldHom.Core.laminate_stress_localization
MeanFieldHom.Core.laminate_gradient_localization
MeanFieldHom.Core.laminate_flux_localization
```

## The cell contract and declarative nesting

```@docs
MeanFieldHom.Core.AbstractHomogenizationCell
MeanFieldHom.Core.validate_cell
MeanFieldHom.Core.cell_member_names
MeanFieldHom.Core.cell_container_property
MeanFieldHom.Core.cell_set_property
MeanFieldHom.Core.Homogenized
MeanFieldHom.Core.NestedParameter
MeanFieldHom.Core.nested
MeanFieldHom.Core.resolve_property
MeanFieldHom.Core.has_nested_property
MeanFieldHom.Core.MAX_NESTING
```
