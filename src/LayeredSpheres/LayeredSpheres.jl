"""
    MeanFieldHom.LayeredSpheres

Isotropic `n`-layer spherical composite inclusion (core + concentric
shells) embedded in an infinite matrix, with perfect or imperfect
interfaces (spring / surface-elastic / Kapitza / surface-conductive).
Public entry points: [`LayeredSphere`](@ref), the bulk / shear
localisation and layer-average utilities.
"""
module LayeredSpheres

using LinearAlgebra
using TensND

import ..Core
using ..Core
const MFH_Core = Core

import ..Elasticity   # for single-layer shear delegation to `Ellipsoid`

# Core-level generics we extend (strain_strain_loc etc. are declared in
# Core/abstractions.jl so every sub-module can attach methods uniformly).
import ..Core: strain_strain_loc, stress_strain_loc, strain_stress_loc,
    stress_stress_loc, gradient_gradient_loc, flux_gradient_loc,
    gradient_flux_loc, flux_flux_loc,
    stiffness_contribution, conductivity_contribution,
    resistivity_contribution, is_homogeneous_inclusion

include("interfaces.jl")
include("geometry.jl")
include("interface_transfer.jl") # interface jump matrices (bulk + shear)
include("bulk_recurrence.jl")    # bulk state-vector transfer + localisation
include("shear_recurrence.jl")   # multi-layer shear stub
include("conductivity.jl")       # Y₁-harmonic conductivity state-vector
include("averages.jl")
include("scheme_integration.jl") # concentration tensors → mean-field schemes

# ── Localisation / contribution overrides for LayeredSphere ─────────────────

"""
    strain_strain_loc(sphere::LayeredSphere, C₀::TensISO{4,3}; layer::Int) -> Tens{4,3}

Per-layer strain-strain localisation tensor in an ISO `LayeredSphere`.
Returns the isotropic 4-tensor `A_k = α_k J + β_k K` for the requested
layer.  `layer` must be in `1..N`.
"""
function strain_strain_loc(
        sphere::LayeredSphere{T, N},
        C₀::TensND.TensISO{4, 3};
        layer::Int,
        kw...,
    ) where {T, N}
    1 ≤ layer ≤ N || throw(BoundsError(sphere, layer))
    κ₀, μ₀ = _iso_bulk_shear(C₀)
    α_k = _bulk_localization(sphere, κ₀, μ₀)[layer]
    β_k = _shear_localization(sphere, C₀)[layer]
    return TensISO{3}(α_k, β_k)
end

# =============================================================================
#  Contribution tensors (iso composite sphere in iso matrix)
# =============================================================================

"""
    stiffness_contribution(sphere, C₀) -> Tens{4,3}

Size-independent stiffness contribution tensor of the composite sphere
relative to the matrix `C₀`.  The dilute-scheme effective stiffness
is `C_eff = C₀ + f · N_sphere` where `f` is the volume fraction of the
composite sphere.  For ISO materials this reduces to two scalar
contributions (bulk + shear):

```
N_bulk  = Σ_k f_k (κ_k - κ₀) α_k      → contributes to the `J` part
N_shear = Σ_k f_k (μ_k - μ₀) β_k      → contributes to the `K` part
```
"""
function Core.stiffness_contribution(
        sphere::LayeredSphere{T, N},
        C₀::TensND.TensISO{4, 3};
        kw...,
    ) where {T, N}
    κ₀, μ₀ = _iso_bulk_shear(C₀)
    κμ = _bulk_layer_moduli(sphere)
    α = _bulk_localization(sphere, κ₀, μ₀)
    β = _shear_localization(sphere, C₀)
    f = ntuple(k -> layer_volume_fraction(sphere, k), Val(N))
    N_bulk = sum(f[k] * (κμ[k][1] - κ₀) * α[k] for k in 1:N)
    N_shear = sum(f[k] * (κμ[k][2] - μ₀) * β[k] for k in 1:N)
    return TensISO{3}(3 * N_bulk, 2 * N_shear)
end

# =============================================================================
#  Conductivity LayeredSphere overrides
# =============================================================================

"""
    gradient_gradient_loc(sphere::LayeredSphere, K₀; layer)

Per-layer gradient-gradient localisation tensor for an isotropic
`LayeredSphere` embedded in an isotropic matrix of conductivity `K₀`.
Returns the scalar `α_k` packed as `TensISO{3}(α_k)` (isotropic
2-tensor), satisfying `<∇T>_layer = α_k · ∇T∞`.
"""
function gradient_gradient_loc(
        sphere::LayeredSphere{T, N},
        K₀::TensND.TensISO{2, 3};
        layer::Int,
        kw...,
    ) where {T, N}
    1 ≤ layer ≤ N || throw(BoundsError(sphere, layer))
    k₀ = _iso_scalar(K₀)
    α_k = _cond_localization(sphere, k₀)[layer]
    return TensISO{3}(α_k)
end

"""
    conductivity_contribution(sphere::LayeredSphere, K₀) -> Tens{2,3}

Size-independent conductivity contribution tensor of the composite
sphere:  `N_K = Σ_k f_k (k_k - k_0) α_k`.
"""
function Core.conductivity_contribution(
        sphere::LayeredSphere{T, N},
        K₀::TensND.TensISO{2, 3};
        kw...,
    ) where {T, N}
    k₀ = _iso_scalar(K₀)
    α = _cond_localization(sphere, k₀)
    k_layers = _cond_layer_moduli(sphere)
    f = ntuple(k -> layer_volume_fraction(sphere, k), Val(N))
    N_K = sum(f[k] * (k_layers[k] - k₀) * α[k] for k in 1:N)
    return TensISO{3}(N_K)
end

# ── Exports ─────────────────────────────────────────────────────────────────
export LayeredSphere, AbstractInterface, PerfectInterface
export SpringInterface, MembraneInterface
export KapitzaInterface, SurfaceConductiveInterface
export layer_count, layer_radius, layer_modulus, layer_interface,
    layer_volume_fraction, outer_radius
export layer_strain_average, sphere_strain_average, cumulative_strain_average
export layer_stiffness_average, layer_compliance_average,
    layer_conductivity_average, layer_resistivity_average

end # module
