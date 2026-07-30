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

Requires a finite-element backend: `Ferrite`, `FerriteGmsh` and `Gmsh`, or —
for the axisymmetric morphology — `Gridap` and `GridapGmsh`. Two morphologies,
one method: see [Finite-element inclusions](@ref man-fe-inclusions) for the
elliptical crack and [A recycled-concrete aggregate](@ref app-recycled-aggregate)
for the sphere with an off-centre core.

```@docs
MeanFieldHom.FiniteElements
MeanFieldHom.FECache
MeanFieldHom.fe_assembly_count
MeanFieldHom.fe_reset!
```

### Choosing a backend

```@docs
MeanFieldHom.FEBackend
MeanFieldHom.AutoBackend
MeanFieldHom.FerriteBackend
MeanFieldHom.GridapBackend
```

### Writing a backend

A backend is these nine methods and nothing else. The Fourier operators, the
boundary data, the fixed point of the corrected boundary condition and the
memoization are shared, and the driver closes the strain operator and the
azimuthal projection over the mode before handing them over — so an
implementation never sees a Fourier mode or a physics, only "this many scalar
fields, this operator, this projection".

`ext/MeanFieldHomGridapExt/` is the shorter of the two implementations, at
about 280 lines, and is the one to read first.

```@docs
MeanFieldHom.FiniteElements._build_gmsh_axi_model
MeanFieldHom.FiniteElements.fe_axi_grid
MeanFieldHom.FiniteElements.fe_axi_grid_counts
MeanFieldHom.FiniteElements.fe_axi_region_volume
MeanFieldHom.FiniteElements.fe_axi_mode
MeanFieldHom.FiniteElements.fe_axi_dof_split
MeanFieldHom.FiniteElements.fe_axi_set_dirichlet!
MeanFieldHom.FiniteElements.fe_axi_stiffness
MeanFieldHom.FiniteElements.fe_axi_average
MeanFieldHom.FiniteElements._resolve_backend
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
