# =============================================================================
#  scheme_integration.jl — plug a `LayeredSphere` into the mean-field schemes.
#
#  A composite sphere has NO Hill tensor: it is not an ellipsoidal
#  inhomogeneity with a uniform eigenstrain.  What it does have — and what the
#  schemes actually need — is a **concentration (localization) tensor**, which
#  the layered recurrences already provide per layer:
#
#      <ε>_k = α_k · ε∞_sph + β_k · ε∞_dev
#
#  The generic `strain_strain_loc(::AbstractInclusion, …)` in `localization.jl`
#  builds `A` from `hill_tensor`, so without the specializations below a
#  `LayeredSphere` phase would fall into it and fail.  The methods here
#  short-circuit that path, exactly as the conductivity side already does.
#
#  Two quantities are needed by the schemes:
#
#    * the whole-inclusion concentration tensor, a volume average over layers
#          A_Ω = (Σ_k f_k α_k) 𝕁 + (Σ_k f_k β_k) 𝕂 ;
#    * the stiffness contribution, which must be assembled **layer by layer**
#          N = Σ_k f_k (C_k − C₀) : A_k ,
#      and is *not* `(C₁ − C₀) : A_Ω` — the inclusion is heterogeneous, so no
#      single `C₁` represents it.
# =============================================================================

"""
    is_homogeneous_inclusion(::LayeredSphere) -> false

A composite sphere has no single representative property: its average stress
must be summed over layers. See [`stress_strain_loc`](@ref).
"""
Core.is_homogeneous_inclusion(::LayeredSphere) = false

"""
    _layer_iso_pairs(sphere) -> NTuple{N, Tuple}

Per-layer `(α, β)` pairs of the isotropic stiffnesses, i.e. `C_k = α_k 𝕁 + β_k 𝕂`.
"""
@inline function _layer_iso_pairs(sphere::LayeredSphere{T, N}) where {T, N}
    return ntuple(k -> TensND.get_data(layer_modulus(sphere, k)), Val(N))
end

@inline function _layer_fractions(sphere::LayeredSphere{T, N}) where {T, N}
    return ntuple(k -> layer_volume_fraction(sphere, k), Val(N))
end

"""
    _layer_localizations(sphere, C₀) -> (α, β, f)

Per-layer bulk (`α_k`) and deviatoric (`β_k`) localization scalars together
with the layer volume fractions.
"""
function _layer_localizations(
        sphere::LayeredSphere{T, N},
        C₀::TensND.TensISO{4, 3},
    ) where {T, N}
    κ₀, μ₀ = _iso_bulk_shear(C₀)
    α = _bulk_localization(sphere, κ₀, μ₀)
    β = _shear_localization(sphere, C₀)
    return α, β, _layer_fractions(sphere)
end

"""
    strain_strain_loc(sphere::LayeredSphere, C₁, C₀; kw...) -> TensISO{4,3}

Whole-inclusion **strain concentration tensor** of a composite sphere embedded
in the isotropic reference `C₀`:

```
A_Ω = (Σ_k f_k α_k) 𝕁 + (Σ_k f_k β_k) 𝕂 ,   <ε>_Ω = A_Ω : ε∞ .
```

`C₁` is accepted for signature compatibility with the generic
`strain_strain_loc(::AbstractInclusion, C₁, C₀)` used by the scheme dispatch,
but is **ignored**: the moduli of a composite sphere live in its layers
(`layer_modulus`), not in a single phase tensor.

For a single layer this reduces exactly to the Eshelby result for a sphere.
"""
function strain_strain_loc(
        sphere::LayeredSphere{T, N},
        ::TensND.AbstractTens{4, 3},
        C₀::TensND.TensISO{4, 3};
        kw...,
    ) where {T, N}
    α, β, f = _layer_localizations(sphere, C₀)
    return TensISO{3}(sum(f[k] * α[k] for k in 1:N), sum(f[k] * β[k] for k in 1:N))
end

"""
    stiffness_contribution(sphere::LayeredSphere, C₁, C₀; kw...) -> TensISO{4,3}

Size-independent **stiffness contribution tensor** of a composite sphere,

```
N_C = Σ_k f_k (C_k − C₀) : A_k .
```

Assembled layer by layer: a composite sphere is heterogeneous, so the usual
`(C₁ − C₀) : A` of a homogeneous inhomogeneity does not apply. `C₁` is ignored,
as in [`strain_strain_loc`](@ref).
"""
function stiffness_contribution(
        sphere::LayeredSphere{T, N},
        ::TensND.AbstractTens{4, 3},
        C₀::TensND.TensISO{4, 3};
        kw...,
    ) where {T, N}
    α, β, f = _layer_localizations(sphere, C₀)
    C_k = _layer_iso_pairs(sphere)
    α₀, β₀ = TensND.get_data(C₀)
    a = sum(f[k] * (C_k[k][1] - α₀) * α[k] for k in 1:N)
    b = sum(f[k] * (C_k[k][2] - β₀) * β[k] for k in 1:N)
    return TensISO{3}(a, b)
end

"""
    stress_strain_loc(sphere::LayeredSphere, C₁, C₀; kw...) -> TensISO{4,3}

Whole-inclusion **average stress** per unit remote strain,

```
⟨C:ε⟩_Ω = (Σ_k f_k C_k : A_k) : ε∞ ,
```

assembled layer by layer. This is what the self-consistent and Mori-Tanaka
kernels need; it is *not* `C₁ : A_Ω`. `C₁` is ignored (see
[`strain_strain_loc`](@ref)).

Consistency: `stiffness_contribution = stress_strain_loc - C₀ : strain_strain_loc`.
"""
function stress_strain_loc(
        sphere::LayeredSphere{T, N},
        ::TensND.AbstractTens{4, 3},
        C₀::TensND.TensISO{4, 3};
        kw...,
    ) where {T, N}
    α, β, f = _layer_localizations(sphere, C₀)
    C_k = _layer_iso_pairs(sphere)
    a = sum(f[k] * C_k[k][1] * α[k] for k in 1:N)
    b = sum(f[k] * C_k[k][2] * β[k] for k in 1:N)
    return TensISO{3}(a, b)
end

"""
    flux_gradient_loc(sphere::LayeredSphere, K₁, K₀; kw...) -> TensISO{2,3}

Conductivity counterpart of [`stress_strain_loc`](@ref):
`⟨k∇T⟩_Ω = (Σ_k f_k k_k α_k) · ∇T∞`. `K₁` is ignored.
"""
function flux_gradient_loc(
        sphere::LayeredSphere{T, N},
        ::TensND.AbstractTens{2, 3},
        K₀::TensND.TensISO{2, 3};
        kw...,
    ) where {T, N}
    α = _cond_localization(sphere, _iso_scalar(K₀))
    k_layers = _cond_layer_moduli(sphere)
    f = _layer_fractions(sphere)
    return TensISO{3}(sum(f[k] * k_layers[k] * α[k] for k in 1:N))
end

"""
    gradient_gradient_loc(sphere::LayeredSphere, K₁, K₀; kw...) -> TensISO{2,3}

Whole-inclusion **gradient concentration tensor**
`α_Ω = Σ_k f_k α_k`, the conductivity counterpart of
[`strain_strain_loc`](@ref). `K₁` is ignored (see there).

The per-layer form is available as
`gradient_gradient_loc(sphere, K₀; layer = k)`.
"""
function gradient_gradient_loc(
        sphere::LayeredSphere{T, N},
        ::TensND.AbstractTens{2, 3},
        K₀::TensND.TensISO{2, 3};
        kw...,
    ) where {T, N}
    α = _cond_localization(sphere, _iso_scalar(K₀))
    f = _layer_fractions(sphere)
    return TensISO{3}(sum(f[k] * α[k] for k in 1:N))
end

"""
    conductivity_contribution(sphere::LayeredSphere, K₁, K₀; kw...) -> TensISO{2,3}

Three-argument form matching the scheme dispatch; `K₁` is ignored and the
computation is delegated to the two-argument method.
"""
function Core.conductivity_contribution(
        sphere::LayeredSphere{T, N},
        ::TensND.AbstractTens{2, 3},
        K₀::TensND.TensISO{2, 3};
        kw...,
    ) where {T, N}
    return Core.conductivity_contribution(sphere, K₀; kw...)
end

"""
    layer_stiffness_average(sphere) -> TensISO{4,3}

Voigt (volume) average of the layer stiffnesses, `Σ_k f_k C_k`. This is what
the Voigt bound needs for a composite sphere: the declared phase property does
not represent it.
"""
function layer_stiffness_average(sphere::LayeredSphere{T, N}) where {T, N}
    C_k = _layer_iso_pairs(sphere)
    f = _layer_fractions(sphere)
    return TensISO{3}(
        sum(f[k] * C_k[k][1] for k in 1:N),
        sum(f[k] * C_k[k][2] for k in 1:N),
    )
end

"""
    layer_compliance_average(sphere) -> TensISO{4,3}

Reuss (volume) average of the layer compliances, `Σ_k f_k C_k⁻¹`.
"""
function layer_compliance_average(sphere::LayeredSphere{T, N}) where {T, N}
    C_k = _layer_iso_pairs(sphere)
    f = _layer_fractions(sphere)
    return TensISO{3}(
        sum(f[k] / C_k[k][1] for k in 1:N),
        sum(f[k] / C_k[k][2] for k in 1:N),
    )
end

"""
    layer_conductivity_average(sphere) -> TensISO{2,3}

Voigt average of the layer conductivities, `Σ_k f_k k_k`.
"""
function layer_conductivity_average(sphere::LayeredSphere{T, N}) where {T, N}
    k_layers = _cond_layer_moduli(sphere)
    f = _layer_fractions(sphere)
    return TensISO{3}(sum(f[k] * k_layers[k] for k in 1:N))
end

"""
    layer_resistivity_average(sphere) -> TensISO{2,3}

Reuss average of the layer resistivities, `Σ_k f_k / k_k`.
"""
function layer_resistivity_average(sphere::LayeredSphere{T, N}) where {T, N}
    k_layers = _cond_layer_moduli(sphere)
    f = _layer_fractions(sphere)
    return TensISO{3}(sum(f[k] / k_layers[k] for k in 1:N))
end
