"""
    MeanFieldHom

Julia package for mean-field homogenisation of heterogeneous materials.

`MeanFieldHom` unifies the computation of Hill polarisation tensors for
ellipsoidal inhomogeneities, crack opening displacement (COD) tensors, stress
and displacement intensity factors, and — in the near future — full
homogenisation schemes, representative volume elements (RVEs), and
viscoelastic constitutive laws, sharing a common abstraction for
inclusions, algorithms, and material symmetry classes.

# Sub-modules

- `MeanFieldHom.Core`         — abstractions (`AbstractInclusion`,
  `AbstractAlgorithm`, `MaterialSymmetry`), shared numerics
  (Green / Newton kernels, Masson-style residue, DECUHR integrand), modulus
  extractors, and central dispatch.
- `MeanFieldHom.Elasticity`   — Hill polarisation for ellipsoidal inclusions
  (2D / 3D, isotropic and anisotropic matrix).
- `MeanFieldHom.Cracks`       — COD tensors, compliance contributions, SIF
  and DIF for elliptic and ribbon cracks.
- `MeanFieldHom.Conductivity` — 2nd-order Hill tensor for conductivity /
  diffusion problems.
- `MeanFieldHom.Schemes`      — placeholder for future mean-field
  homogenisation schemes (dilute, Mori–Tanaka, self-consistent, PCW, …).

# Shared generic interface

```julia
hill_tensor(ell::AbstractEllipsoidalInclusion, C₀; method=:auto, ...)
cod_tensor(crack::AbstractCrack, C₀; method=:auto, ...)
compliance_contribution(crack, C₀, ε; method=:auto, ...)
sif(crack, C₀, Σ; method=:auto, ...)
dif(crack, C₀, Σ; method=:auto, ...)
```

All high-level entry points share the same algorithmic traits
(`Analytical`, `Residue`, `DECUHR`) and the same material-symmetry dispatch
rules. See the developer documentation (`docs/src/developer/`) for guidance
on extending the package with new inclusions, algorithms or schemes.
"""
module MeanFieldHom

include("Core/Core.jl")
include("Elasticity/Elasticity.jl")
include("Cracks/Cracks.jl")
include("Conductivity/Conductivity.jl")
include("Schemes/Schemes.jl")

using .Core
using .Elasticity
using .Cracks
using .Conductivity
using .Schemes

# ── Abstractions ─────────────────────────────────────────────────────────────
export AbstractInclusion, AbstractEllipsoidalInclusion, AbstractCrack
export AbstractLayeredInclusion
export AbstractAlgorithm, Analytical, Residue, DECUHR, Auto
export MaterialSymmetry, IsotropicSym, TransverselyIsotropicSym,
    OrthotropicSym, GeneralAnisotropicSym
export material_symmetry, dimension, inclusion_basis, shape_trait

# ── Elasticity ───────────────────────────────────────────────────────────────
export Ellipsoid
export EllipsoidShape, Spherical, Prolate, Oblate, Triaxial, Circular, Elliptic
export tens_IA, tens_UA, tens_VA
export hill_tensor

# ── Cracks ───────────────────────────────────────────────────────────────────
export CrackShape, Penny, EllipticShape, Ribbon
export EllipticCrack, RibbonCrack, PennyCrack
export crack_basis, aspect_ratio, semi_major, semi_minor, crack_normal
export cod_tensor, B_tensor
export cod_from_compliance, compliance_from_cod, cod_from_deltaS
export compliance_contribution
export sif, dif

# ── Backwards-compat aliases ─────────────────────────────────────────────────
const HillAlgorithm = AbstractAlgorithm
const CrackAlgorithm = AbstractAlgorithm
export HillAlgorithm, CrackAlgorithm

end # module
