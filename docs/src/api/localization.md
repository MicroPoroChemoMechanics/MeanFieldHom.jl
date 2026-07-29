# API — Localization & contribution

## Localization tensors

```@docs
MeanFieldHom.strain_strain_loc
MeanFieldHom.stress_strain_loc
MeanFieldHom.strain_stress_loc
MeanFieldHom.stress_stress_loc
MeanFieldHom.gradient_gradient_loc
MeanFieldHom.flux_gradient_loc
MeanFieldHom.gradient_flux_loc
MeanFieldHom.flux_flux_loc
```

## Contribution tensors

```@docs
MeanFieldHom.stiffness_contribution
MeanFieldHom.compliance_contribution
MeanFieldHom.conductivity_contribution
MeanFieldHom.resistivity_contribution
```

## The amount × contribution seam

```@docs
MeanFieldHom.delta_stiffness
MeanFieldHom.delta_compliance
MeanFieldHom.delta_conductivity
MeanFieldHom.delta_resistivity
```

For a flat inclusion the four three-argument seams share one geometric
prefactor, [`crack_density_factor`](@ref MeanFieldHom.Cracks.crack_density_factor).

## Inclusion traits and bundled seams

```@docs
MeanFieldHom.is_homogeneous_inclusion
MeanFieldHom.loc_and_stiffness
MeanFieldHom.loc_and_stress_average
MeanFieldHom.compliance_and_stiffness_contribution
```

## Custom (user-defined) inclusions

See [Custom inclusions](@ref man-custom-inclusions) for the tutorial and
[Adding a new inclusion](@ref) for the full contract.

```@docs
MeanFieldHom.CustomInclusions
MeanFieldHom.Core.AbstractCustomInclusion
MeanFieldHom.CustomInclusion
MeanFieldHom.CustomShape
MeanFieldHom.check_inclusion_interface
```

## Finite-element inclusions

Requires `Ferrite`, `FerriteGmsh` and `Gmsh`. Two morphologies, one method —
see [Finite-element inclusions](@ref man-fe-inclusions) for the elliptical
crack and [A recycled-concrete aggregate](@ref app-recycled-aggregate) for the
sphere with an off-centre core.

```@docs
MeanFieldHom.FiniteElements
MeanFieldHom.FECache
MeanFieldHom.fe_assembly_count
MeanFieldHom.fe_reset!
```

### Elliptical crack (3-D)

```@docs
MeanFieldHom.FEEllipticCrack
MeanFieldHom.FEMeshOptions
MeanFieldHom.fe_mesh_report
MeanFieldHom.fe_cod_breakdown
```

### Sphere with an off-centre core (axisymmetric Fourier)

```@docs
MeanFieldHom.FEExcenteredSphere
MeanFieldHom.FEAxiMeshOptions
MeanFieldHom.fe_axi_localization
MeanFieldHom.fe_axi_breakdown
MeanFieldHom.fe_axi_mesh_report
MeanFieldHom.FiniteElements.core_radius
MeanFieldHom.FiniteElements.core_offset
MeanFieldHom.FiniteElements.tensor_order
MeanFieldHom.FiniteElements.ExcenteredSphereShape
```

### Green function of the corrected boundary condition

```@docs
MeanFieldHom.Core.green_gradient_iso
MeanFieldHom.Core.dipole_displacement_iso
```
