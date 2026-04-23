# =============================================================================
#  averages.jl — layer / sphere / cumulative averages of the strain (or
#  gradient) field inside an isotropic `LayeredSphere`.
#
#  For the bulk part (fully supported for any `N`), the volume
#  average in layer `k` of the strain tensor reduces to
#
#      <ε>_k = α_k · ε∞   for purely hydrostatic ε∞,
#
#  where `α_k` is the per-layer bulk localisation.  For a general remote
#  strain, the decomposition splits bulk + deviatoric and the shear
#  contribution is delegated to the multi-layer shear solver (single
#  layer only for now).
# =============================================================================

"""
    layer_strain_average(sphere, C₀, ε∞, layer) -> Tens{2,3}

Volume-averaged strain tensor `<ε>_layer` inside the `layer`-th layer
of a `LayeredSphere` embedded in an isotropic matrix `C₀`, under a
remote strain `ε∞`.  Returns a symmetric 2-tensor in the canonical
frame.  Combines the bulk localisation `α_k` (hydrostatic part) and
the shear localisation `β_k` (deviatoric part).
"""
function layer_strain_average(
        sphere::LayeredSphere{T, N},
        C₀::TensND.TensISO{4, 3},
        ε∞::TensND.AbstractTens{2, 3},
        layer::Int,
    ) where {T, N}
    1 ≤ layer ≤ N || throw(BoundsError(sphere, layer))
    κ₀, μ₀ = _iso_bulk_shear(C₀)
    α = _bulk_localization(sphere, κ₀, μ₀)[layer]
    β = _shear_localization(sphere, C₀)[layer]

    Tres = promote_type(T, eltype(C₀), eltype(ε∞))
    I2 = TensISO{3}(one(Tres))
    tr_ε∞ = sum(ε∞[i, i] for i in 1:3)
    ε_sph = (tr_ε∞ / 3) * I2
    ε_dev = ε∞ - ε_sph
    return α * ε_sph + β * ε_dev
end

"""
    sphere_strain_average(sphere, C₀, ε∞) -> Tens{2,3}

Volume-averaged strain over the whole composite sphere (all layers
combined): `<ε>_Ω = Σ_k f_k <ε>_k` where `f_k` is the volume fraction
of layer `k` inside the composite sphere.
"""
function sphere_strain_average(
        sphere::LayeredSphere{T, N},
        C₀::TensND.TensISO{4, 3},
        ε∞::TensND.AbstractTens{2, 3},
    ) where {T, N}
    f = ntuple(k -> layer_volume_fraction(sphere, k), N)
    avgs = ntuple(k -> layer_strain_average(sphere, C₀, ε∞, k), N)
    return sum(f[k] * avgs[k] for k in 1:N)
end

"""
    cumulative_strain_average(sphere, C₀, ε∞, r) -> Tens{2,3}

Volume-averaged strain over the ball of radius `r ∈ (0, r_N]` centred
on the composite sphere centre.  The ball may cross several layers;
the result is the volume-weighted average of the per-layer averages
truncated by the final partial layer.
"""
function cumulative_strain_average(
        sphere::LayeredSphere{T, N},
        C₀::TensND.TensISO{4, 3},
        ε∞::TensND.AbstractTens{2, 3},
        r,
    ) where {T, N}
    r > 0 || throw(ArgumentError("cumulative_strain_average radius must be > 0"))
    radii = sphere.radii

    # Accumulate the "volume × average" contribution layer by layer.
    Tres = promote_type(T, typeof(r), eltype(C₀), eltype(ε∞))
    acc_vol_times_avg = nothing
    total_vol = zero(Tres)

    for k in 1:N
        r_prev = k == 1 ? zero(Tres) : radii[k - 1]
        r_k = radii[k]
        if r ≤ r_prev
            break   # ball no longer reaches into this layer
        end
        r_upper = min(Tres(r), Tres(r_k))
        vol_k = (4 * π / 3) * (r_upper^3 - r_prev^3)
        avg_k = layer_strain_average(sphere, C₀, ε∞, k)
        acc_vol_times_avg = acc_vol_times_avg === nothing ?
            vol_k * avg_k :
            acc_vol_times_avg + vol_k * avg_k
        total_vol += vol_k
        if r ≤ r_k
            break   # ball does not extend beyond this layer
        end
    end

    acc_vol_times_avg === nothing &&
        throw(ArgumentError("cumulative_strain_average: ball is empty"))
    return (1 / total_vol) * acc_vol_times_avg
end
