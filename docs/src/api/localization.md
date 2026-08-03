# [API — Localization & contribution](@id api-localization)

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
MeanFieldHom.Core.is_homogeneous_inclusion
MeanFieldHom.loc_and_stiffness
MeanFieldHom.loc_and_stress_average
MeanFieldHom.compliance_and_stiffness_contribution
```

## Custom (user-defined) inclusions

See [Custom inclusions](@ref man-custom-inclusions) for the tutorial and
[Adding a new inclusion](@ref dev-adding-inclusion) for the full contract.

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
for the sphere with an off-center core.

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

A backend is sixteen methods and nothing else — nine for the axisymmetric
solve, seven for the crack. The Fourier operators, the boundary data, the fixed
point of the corrected boundary condition and the memoization are shared, and
the driver closes the strain operator and the azimuthal projection over the
mode before handing them over, so an implementation never sees a Fourier mode
or a physics: only "this many scalar fields, this operator, this projection".

`ext/MeanFieldHomGridapExt/` is the shorter of the two implementations and the
one to read first.

### The axisymmetric solve

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

### The crack

```@docs
MeanFieldHom.FiniteElements._build_gmsh_crack_model
MeanFieldHom.FiniteElements._weld_msh_crack_front
MeanFieldHom.FiniteElements.fe_crack_grid
MeanFieldHom.FiniteElements.fe_crack_counts
MeanFieldHom.FiniteElements.fe_crack_space
MeanFieldHom.FiniteElements.fe_crack_dof_split
MeanFieldHom.FiniteElements.fe_crack_set_dirichlet!
MeanFieldHom.FiniteElements.fe_crack_stiffness
MeanFieldHom.FiniteElements.fe_crack_mean_jump
```

### Elliptical crack (3-D)

```@docs
MeanFieldHom.FEEllipticCrack
MeanFieldHom.FEMeshOptions
MeanFieldHom.fe_mesh_report
MeanFieldHom.fe_cod_breakdown
```

### Sphere with an off-center core (axisymmetric Fourier)

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

## Neural-surrogate inclusions

The fourth route into the contract: the response comes out of a trained network.
See [Neural-surrogate inclusions](@ref man-neural-inclusions) for the tutorial.
Evaluating needs nothing beyond the package; *training* needs
`import Lux, Optimisers, Zygote`.

```@docs
MeanFieldHom.NeuralInclusions
MeanFieldHom.NeuralHillInclusion
MeanFieldHom.NeuralLocalizationInclusion
MeanFieldHom.NeuralInclusions.NeuralShape
```

### The surrogate

```@docs
MeanFieldHom.NeuralSurrogate
MeanFieldHom.Provenance
MeanFieldHom.worst_error
MeanFieldHom.NeuralInclusions.check_domain
MeanFieldHom.NeuralInclusions.predict_components
```

### What the network predicts

The symmetry class, the major symmetry, the homogeneity in the reference moduli
and the frame are *enforced* by these types rather than fitted — see
[What is exact, and what is fitted](@ref man-neural-inclusions).

```@docs
MeanFieldHom.NeuralInclusions.AbstractHillClass
MeanFieldHom.HillISO
MeanFieldHom.HillTI
MeanFieldHom.HillOrtho
MeanFieldHom.HillISO2
MeanFieldHom.HillTI2
MeanFieldHom.NeuralInclusions.AbstractOutputSpec
MeanFieldHom.DimensionlessHill
MeanFieldHom.AffineHill
MeanFieldHom.NeuralInclusions.ncomponents
MeanFieldHom.NeuralInclusions.tensor_order
MeanFieldHom.NeuralInclusions.nterms
MeanFieldHom.NeuralInclusions.noutputs
MeanFieldHom.NeuralInclusions.needs_nu
MeanFieldHom.NeuralInclusions.build
MeanFieldHom.NeuralInclusions.components
MeanFieldHom.NeuralInclusions.decode
MeanFieldHom.NeuralInclusions.material_coeffs
MeanFieldHom.NeuralInclusions.dimensionless_scale
MeanFieldHom.NeuralInclusions.hill_class
MeanFieldHom.NeuralInclusions.output_spec
MeanFieldHom.NeuralInclusions.apply_transform
MeanFieldHom.NeuralInclusions.invert_transform
MeanFieldHom.NeuralInclusions._feature
MeanFieldHom.NeuralInclusions.raw_features
MeanFieldHom.NeuralInclusions._class_frame
MeanFieldHom.NeuralInclusions._canonical_axes
MeanFieldHom.NeuralInclusions._spheroid_axis_index
```

### Sampling and labeling

```@docs
MeanFieldHom.SampleBox
MeanFieldHom.Dataset
MeanFieldHom.generate_dataset
MeanFieldHom.NeuralInclusions.sample_box
MeanFieldHom.NeuralInclusions.grid_box
MeanFieldHom.NeuralInclusions.halton
MeanFieldHom.NeuralInclusions.feature_index
MeanFieldHom.fit_scaling
MeanFieldHom.validate_surrogate
MeanFieldHom.report_surrogate
MeanFieldHom.component_labels
```

### Training

`train_surrogate` is the seam of the `MeanFieldHomLuxExt` extension: the method
below is the fallback that raises when the extension is not loaded.

```@docs
MeanFieldHom.TrainingOptions
MeanFieldHom.train_surrogate
MeanFieldHom.assemble_surrogate
MeanFieldHom.NeuralInclusions.network_widths
```

### The network, and its serialization

```@docs
MeanFieldHom.NeuralInclusions.MLP
MeanFieldHom.NeuralInclusions.NNDense
MeanFieldHom.NeuralInclusions.glorot_mlp
MeanFieldHom.NeuralInclusions.softplus
MeanFieldHom.NeuralInclusions.activation
MeanFieldHom.NeuralInclusions.activation_name
MeanFieldHom.NeuralInclusions.layer_widths
MeanFieldHom.NeuralInclusions.layer_activations
MeanFieldHom.NeuralInclusions.nparams
MeanFieldHom.save_surrogate
MeanFieldHom.load_surrogate
MeanFieldHom.model_path
MeanFieldHom.shipped_models
MeanFieldHom.NeuralInclusions.SURROGATE_FORMAT
MeanFieldHom.NeuralInclusions.MODEL_DIR
```
