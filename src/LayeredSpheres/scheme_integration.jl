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
    _membrane_surface_stress(sphere, C₀) -> (a_surf, b_surf)

Contribution of Gurtin–Murdoch surface stress on the dual
([`MembraneInterface`](@ref)) interfaces to the volume-averaged stress of
the composite sphere, per unit remote strain, split into bulk (`𝕁`) and
shear (`𝕂`) scalars.  From the average-stress theorem with a coherent
surface, `⟨σ⟩_Ω = Σ_k f_k C_k:A_k + (1/V) Σ_Γ ∮_Γ σˢ dS`; the surface
integrals are (with `κs = λs + μs`, `r` the interface radius, `R` the
outer radius):

```
bulk :  4 κs · u_r(r) · r / R³
shear:  (−6κs U + 18κs W + 36μs W) · r / (5 R³)         (× 3/2, 𝕂-projection)
```

`u_r(r)` is the bulk radial amplitude (normalised by the far-field `A∞`);
`U(r), W(r)` are the deviatoric displacement amplitudes at the interface
(already normalised to a unit remote deviatoric far field).
"""
function _membrane_surface_stress(
        sphere::LayeredSphere{T, N}, C₀::TensND.TensISO{4, 3}
    ) where {T, N}
    κ₀, μ₀ = _iso_bulk_shear(C₀)
    radii = sphere.radii
    R³ = radii[N]^3

    # Any membrane interface present?  (cheap short-circuit)
    has_membrane = any(k -> layer_interface(sphere, k) isa MembraneInterface, 1:N)
    has_membrane || return (zero(T), zero(T))

    # Bulk amplitudes u_r(r_k), normalised by the far-field A∞.
    inside_b, s_b = _bulk_state_seq(sphere, κ₀, μ₀)
    A_inf, _ = _bulk_extract_AB(radii[N], κ₀, μ₀, s_b[1], s_b[2])
    # Deviatoric state amplitudes (U, W)(r_k), already at unit remote far field.
    states_s, _ = _shear_state_seq(sphere, C₀)

    a_surf = zero(promote_type(T, typeof(κ₀), typeof(A_inf)))
    b_surf = zero(a_surf)
    for k in 1:N
        intf = layer_interface(sphere, k)
        intf isa MembraneInterface || continue
        κs = intf.κs; μs = intf.μs
        r = radii[k]
        u_r = inside_b[k][1] / A_inf
        a_surf += 4 * κs * u_r * r / R³
        U = states_s[k][1]; W = states_s[k][2]
        # (σzz−σxx)/V = C·r/R³;  𝕂-amplitude b = (σzz−σxx)/3.
        C = (-6 * κs * U + 18 * κs * W + 36 * μs * W) / 5
        b_surf += (C * r / R³) / 3
    end
    return a_surf, b_surf
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
    a_surf, b_surf = _membrane_surface_stress(sphere, C₀)
    return TensISO{3}(a + a_surf, b + b_surf)
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
    a_surf, b_surf = _membrane_surface_stress(sphere, C₀)
    return TensISO{3}(a + a_surf, b + b_surf)
end

"""
    flux_gradient_loc(sphere::LayeredSphere, K₁, K₀; kw...) -> TensISO{2,3}

Conductivity counterpart of [`stress_strain_loc`](@ref):
`⟨k∇T⟩_Ω = (Σ_k f_k k_k α_k) · ∇T∞`, plus the surface-conduction flux
[`_cond_surface_flux`](@ref) of any dual (surface-conductive) interface.
`K₁` is ignored.
"""
function flux_gradient_loc(
        sphere::LayeredSphere{T, N},
        ::TensND.AbstractTens{2, 3},
        K₀::TensND.TensISO{2, 3};
        kw...,
    ) where {T, N}
    k₀ = _iso_scalar(K₀)
    α = _cond_localization(sphere, k₀)
    k_layers = _cond_layer_moduli(sphere)
    f = _layer_fractions(sphere)
    return TensISO{3}(
        sum(f[k] * k_layers[k] * α[k] for k in 1:N) + _cond_surface_flux(sphere, k₀)
    )
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
