# [Roadmap](@id dev-roadmap)

## Shipped

- Mean-field schemes: Voigt/Reuss bounds, dilute, Mori–Tanaka, Maxwell,
  Ponte-Castañeda–Willis, self-consistent (Anderson + Newton),
  asymmetric self-consistent, differential.
- Representative volume element (RVE) assembly and effective-property
  pipelines mirroring the reference C++ RVE assembly.
- Concentric multi-layer sphere (`LayeredSphere` via
  [`AbstractLayeredInclusion`](@ref)): Hervé-Zaoui bulk / shear /
  conductivity recurrences, five interface types (perfect, spring,
  membrane, Kapitza, surface-conductive), volume-average and pointwise
  localization fields.
- Ageing linear viscoelasticity (ALV): time-domain Volterra pipeline for
  every scheme, structured ISO/TI/ortho fast paths, ALV cracks and the
  ALV layered sphere (bulk **and** shear recurrences).
- Exact rotation-group symmetrization (ISO / TI) of concentration tensors,
  preserving non-major-symmetric content (`TensTI{4,T,8}`), for arbitrary
  multi-axis orientation distributions inside every scheme kernel.
- User-defined inclusions and algorithms: a levelled, documented contract
  ([Adding a new inclusion](@ref dev-adding-inclusion)), the neutral
  [`AbstractCustomInclusion`](@ref) branch, the callback-driven
  [`CustomInclusion`](@ref MeanFieldHom.CustomInclusion), the
  [`check_inclusion_interface`](@ref MeanFieldHom.check_inclusion_interface)
  conformance checker, and `shape_trait`-based inheritance of the crack
  algebra (a user crack needs only `cod_tensor`).
- Real-space Kelvin Green gradient and dipole far field for an isotropic
  matrix ([`green_gradient_iso`](@ref MeanFieldHom.Core.green_gradient_iso),
  [`dipole_displacement_iso`](@ref MeanFieldHom.Core.dipole_displacement_iso)) —
  the boundary correction that makes a finite numerical Eshelby cell behave
  like an infinite medium.
- ForwardDiff sensitivities across all elastic and ALV schemes (fractions,
  moduli, and inclusion geometry).
- NonlinearSolve.jl backend for the self-consistent fixed point
  (`MeanFieldHomNonlinearSolveExt`): any SciML algorithm
  (`NewtonRaphson`, `TrustRegion`, …) can solve `SelfConsistent` /
  `AsymmetricSelfConsistent`, through an implicit-function-theorem lift
  that keeps `derivative`/`gradient`/`jacobian` exact and free of nested
  `ForwardDiff.Dual`s regardless of algorithm.

## Open

- Extended-COD crack model in conduction: resistive cracks (linear-spring
  analog) **and** conductive cracks (elastic-membrane analog), via a
  tensorial conduction COD.
- Multi-layer extensions: coated cylinders, anisotropic per-layer moduli,
  excentered spheres.
- `PairwiseDistribution` (Willis 1982) envelope for the PCW scheme.
- Native Anderson acceleration with memory > 1, replacing the current
  `AndersonDefault` (currently Picard with relaxation, memory = 1).
- Optional structured `TensTI{4,T,8}` fast path for the ALV TI schemes.
- Viscoelastic constitutive laws in the Laplace–Carson domain.
- Finite-element inclusions, behind the `FEBackend` contract
  (`MeanFieldHomFerriteExt`, `MeanFieldHomGridapExt`), both with the
  first-order corrected boundary condition of
  [adessinaIJES2017](@cite) and an
  isotropic reference medium: the **elliptical crack** in 3-D tetrahedra
  (3 + 3 crack declination) and the **sphere with an off-centre core** in
  axisymmetric Fourier elements (the general polarization fixed point).
  Open extensions — anisotropic reference medium (Pan-Chou or Barnett-Willis
  Green gradient); more than one inclusion, or a non-spherical envelope, in the
  axisymmetric cell; automatic differentiation through the solve.
