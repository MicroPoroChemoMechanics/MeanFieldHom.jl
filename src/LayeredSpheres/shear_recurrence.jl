# =============================================================================
#  shear_recurrence.jl — deviatoric (Y₂-harmonic) recurrence for the
#  isotropic n-layer spherical inclusion, in the (κ, μ) parametrisation
#  that stays regular in the incompressibility limit (κ → ∞).
#
#  Under a remote pure-deviatoric loading the displacement field in an
#  isotropic layer has the axisymmetric form
#
#      u_r(r, θ) = U(r) · P₂(cos θ),
#      u_θ(r, θ) = V(r) · dP₂(cos θ)/dθ.
#
#  The four linearly-independent solutions of the coupled Navier ODEs on
#  (U, V) have radial dependencies (r, r³, 1/r⁴, 1/r²) with ratios that
#  involve only rational functions of (κ, μ) — every coefficient stays
#  finite as κ → ∞.  We use the "Christensen-Lo 1979 / Love" basis:
#
#      Mode 1 : U = r,         V = r
#      Mode 2 : U = D₁ r³,     V = D₃ r³
#      Mode 3 : U = 3/r⁴,      V = -2/r⁴
#      Mode 4 : U = D₂/r²,     V = D₄/r²
#
#  with D₁ = 2 - 3x, D₂ = 3(x + 1), D₃ = (15x + 11)/3, D₄ = 2(3x + 1)/3
#  and x = κ/μ.
#
#  State vector at radius r :    S(r) = (U, W, σ_rr, σ_rθ)
#  where (U, W) are the radial and tangential-amplitude displacement
#  components (u_r = U P₂, u_θ = W dP₂/dθ) and (σ_rr, σ_rθ) are the
#  physical traction amplitudes (not divided by μ) — they are
#  continuous across every perfect interface.  From Hooke
#  σ_ij = λ δ_ij tr(ε) + 2μ ε_ij applied to each mode:
#      σ_rr_k amp = (λ + 2μ) U_k' + (2λ/r) (U_k − 3 W_k)
#      σ_rθ_k amp = μ (W_k' + (U_k − W_k) / r)
#
#  giving the 4×4 fundamental matrix `_shear_M_matrix(r, κ, μ)` whose
#  columns are the four Navier modes evaluated at r in the basis
#  (U, W, σ_rr, σ_rθ).
#
#  Intra-layer field transfer:
#     S(r_out) = M(r_out; κ, μ) · M(r_in; κ, μ)⁻¹ · S(r_in).
#  Interface jumps are applied via `_shear_interface_T(intf, κ, μ, r)`
#  from `interface_transfer.jl`.
#
#  Seed: finiteness at the origin forces the two singular amplitudes
#  c₁ = d₁ = 0 in the core, leaving two free regular amplitudes (a, b).
#  We propagate two probe states (a, b) = (1, 0) and (0, 1) through the
#  whole stack, extract the matrix-side (a, b) amplitudes from the state
#  at r_N⁺, and form the linear combination that yields a unit remote
#  deviatoric far-field (a_{N+1}, b_{N+1}) = (1, 0).  The per-layer
#  shear localisation is β_k = a_k (volume-averaged Y₂ projection).
# =============================================================================

"""
    _shear_M_matrix(r, κ, μ) -> Matrix(4×4)

Fundamental 4×4 matrix of the Y₂-harmonic deviatoric problem.  Columns
1..4 are the four modes (r, r³, 1/r⁴, 1/r²) evaluated at `r` in an
isotropic layer of moduli `(κ, μ)`; rows are the state vector
`(U, V, τ_rr, τ_rθ)` with `τ = σ/μ`.

All entries are rational in `(κ, μ, r)`; no `1/(1−2ν)` remains, so the
matrix is finite in the incompressibility limit `κ → ∞`.
"""
@inline function _shear_M_matrix(r, κ, μ)
    T  = promote_type(typeof(r), typeof(κ), typeof(μ))
    Tκ = T(κ); Tμ = T(μ); Tr = T(r)
    x  = Tκ / Tμ

    r²     = Tr * Tr
    r³     = r² * Tr
    inv_r² = one(T) / r²
    inv_r³ = inv_r² / Tr
    inv_r⁴ = inv_r² * inv_r²
    inv_r⁵ = inv_r⁴ / Tr

    # Displacement-mode ratios derived directly from the Navier
    # characteristic equation for ℓ = 2:
    #   n = 1 :   U/W =  2                       (uniform deviatoric strain)
    #   n = 3 :   U/W =  6(3x − 2)/(15x + 11)    with x = κ/μ
    #   n = -4:   U/W = -3
    #   n = -2:   U/W =  3(x + 1)

    M = Matrix{T}(undef, 4, 4)

    # Mode 1 — (U, W) = (2r, r), uniform deviatoric strain.
    M[1, 1] = 2 * Tr
    M[2, 1] = Tr
    M[3, 1] = 4 * Tμ
    M[4, 1] = 2 * Tμ

    # Mode 2 — (U, W) = (6(3x−2) r³, (15x+11) r³).
    α₂ = 6 * (3 * x - 2)
    γ₂ = 15 * x + 11
    M[1, 2] = α₂ * r³
    M[2, 2] = γ₂ * r³
    M[3, 2] = 6 * (2 - 3 * x) * Tμ * r²
    M[4, 2] = 2 * (24 * x + 5) * Tμ * r²

    # Mode 3 — (U, W) = (3/r⁴, -1/r⁴).
    M[1, 3] =  3 * inv_r⁴
    M[2, 3] =     -inv_r⁴
    M[3, 3] = -24 * Tμ * inv_r⁵
    M[4, 3] =   8 * Tμ * inv_r⁵

    # Mode 4 — (U, W) = (3(x+1)/r², 1/r²).
    α₄ = 3 * (x + 1)
    M[1, 4] = α₄ * inv_r²
    M[2, 4] =      inv_r²
    M[3, 4] = -2 * (9 * x + 4) * Tμ * inv_r³
    M[4, 4] =  3 * x * Tμ * inv_r³

    return M
end

"""
    _shear_layer_transfer(r_out, r_in, κ, μ) -> Matrix(4×4)

Intra-layer field-to-field transfer `S(r_out) = T · S(r_in)` computed
as `M(r_out) · M(r_in)⁻¹`.
"""
@inline function _shear_layer_transfer(r_out, r_in, κ, μ)
    M_in  = _shear_M_matrix(r_in,  κ, μ)
    M_out = _shear_M_matrix(r_out, κ, μ)
    return M_out / M_in
end

"""
    _shear_seed_states(r_1, κ_1, μ_1) -> (probe_a, probe_b)

Two independent probe state vectors at `r = r_1⁻` corresponding to the
two regular amplitudes `(a₁, b₁) = (1, 0)` and `(0, 1)` (with
`c₁ = d₁ = 0` forced by finiteness at the origin).  Returned as the
matching columns of `M(r_1; κ_1, μ_1)`.
"""
@inline function _shear_seed_states(r_1, κ_1, μ_1)
    T = promote_type(typeof(r_1), typeof(κ_1), typeof(μ_1))
    M1 = _shear_M_matrix(T(r_1), T(κ_1), T(μ_1))
    probe_a = T[M1[1, 1], M1[2, 1], M1[3, 1], M1[4, 1]]
    probe_b = T[M1[1, 2], M1[2, 2], M1[3, 2], M1[4, 2]]
    return probe_a, probe_b
end

"""
    _shear_extract_amplitudes(r, κ, μ, state) -> (a, b, c, d)

Given the state `(U, V, τ_rr, τ_rθ)` at radius `r` in a layer of
moduli `(κ, μ)`, return the local mode amplitudes `(a, b, c, d)` by
solving `M(r; κ, μ) · x = state`.
"""
@inline function _shear_extract_amplitudes(r, κ, μ, state)
    T = promote_type(typeof(r), typeof(κ), typeof(μ), eltype(state))
    M = _shear_M_matrix(T(r), T(κ), T(μ))
    return M \ Vector{T}(state)
end

"""
    _shear_state_seq(sphere, C₀) -> NTuple{N, state⁻}, state⁺_N

Propagate two linearly-independent probe state vectors from the core
outward through every interface and every intermediate layer, then
form the linear combination that matches the remote far-field
`(a_{N+1}, b_{N+1}) = (1, 0)`.  Returns the composite state sequence
inside every layer (at its inner-interface radius `r_{k}⁻`) and the
state on the matrix side of the outer interface (`r_N⁺`).
"""
function _shear_state_seq(sphere::LayeredSphere{T, N}, C₀::TensND.TensISO{4, 3}) where {T, N}
    κμ = _bulk_layer_moduli(sphere)
    κ₀, μ₀ = _iso_bulk_shear(C₀)
    TP = _bulk_promote(sphere, κμ, κ₀, μ₀)
    radii = sphere.radii

    κ1, μ1 = κμ[1]
    sa, sb = _shear_seed_states(TP(radii[1]), TP(κ1), TP(μ1))

    inside_a = Vector{Vector{TP}}(undef, N)
    inside_b = Vector{Vector{TP}}(undef, N)
    for k in 1:N
        # `sa`, `sb` are rebound (not mutated) on each subsequent step,
        # so capturing the current vector pointer suffices — no `copy`.
        inside_a[k] = sa
        inside_b[k] = sb
        intf = layer_interface(sphere, k)
        (κk, μk) = κμ[k]
        Tint = _shear_interface_T(intf, TP(κk), TP(μk), TP(radii[k]))
        sa = Tint * sa
        sb = Tint * sb
        if k < N
            (κk1, μk1) = κμ[k + 1]
            Tlay = _shear_layer_transfer(TP(radii[k + 1]), TP(radii[k]),
                                         TP(κk1), TP(μk1))
            sa = Tlay * sa
            sb = Tlay * sb
        end
    end
    # Solve for the linear combination of the two probes that yields
    # (a_matrix, b_matrix) = (1, 0) in the outer matrix at r_N⁺.
    (aa, ab, _, _) = _shear_extract_amplitudes(TP(radii[N]), TP(κ₀), TP(μ₀), sa)
    (ba, bb, _, _) = _shear_extract_amplitudes(TP(radii[N]), TP(κ₀), TP(μ₀), sb)
    det_ab = aa * bb - ab * ba
    λa =  bb / det_ab
    λb = -ab / det_ab

    states = ntuple(k -> λa * inside_a[k] + λb * inside_b[k], N)
    s_matrix = λa * sa + λb * sb
    return states, s_matrix
end

"""
    _shear_localization_single_layer(sphere, C₀) -> β::T

For a single-layer (`N = 1`) composite sphere, delegate the deviatoric
localisation to the existing `Ellipsoid(r)` Eshelby machinery.
"""
function _shear_localization_single_layer(
        sphere::LayeredSphere{T, 1}, C₀::TensND.TensISO{4, 3}
    ) where {T}
    C₁ = layer_modulus(sphere, 1)
    ell = Elasticity.Ellipsoid(outer_radius(sphere))
    A = strain_strain_loc(ell, C₁, C₀)
    _, β_K = A.data
    return β_K
end

"""
    _layer_avg_dev_shear_factor(r_a, r_b, κ, μ) -> Number

Per-unit mode-2 amplitude `b` contribution to the layer-volume-averaged
deviatoric strain in a spherical shell `(r_a, r_b)` (with `r_a = 0` for
the innermost layer) of moduli `(κ, μ)`.  Equals
`(21/5) · (3κ + μ)/μ · (r_b⁵ - r_a⁵) / (r_b³ - r_a³)` (Christensen-Lo
mode-2 angular integral; modes 3 and 4 contribute zero to the dev β).

The full per-layer dev localisation is therefore
`β_k = a_k + b_k · _layer_avg_dev_shear_factor(r_a, r_b, κ_k, μ_k)`.
"""
@inline function _layer_avg_dev_shear_factor(r_a, r_b, κ, μ)
    T = promote_type(typeof(r_a), typeof(r_b), typeof(κ), typeof(μ))
    Tκ = T(κ); Tμ = T(μ); Tra = T(r_a); Trb = T(r_b)
    Trb3 = Trb^3; Tra3 = Tra^3
    Trb5 = Trb^5; Tra5 = Tra^5
    geom = (Trb5 - Tra5) / (Trb3 - Tra3)
    return T(21//5) * (3 * Tκ + Tμ) / Tμ * geom
end

"""
    _shear_localization_multi(sphere, C₀) -> NTuple{N}

Multi-layer (`N ≥ 2`) per-layer deviatoric localisation `β_k` from the
4×4 state-vector recurrence.  For a spherical shell layer, the
volume-averaged deviatoric strain involves both the mode-1 amplitude
(uniform deviatoric part) **and** the mode-2 amplitude (whose r³
displacement profile contributes a non-zero integrated dev strain
through the layer thickness).  Modes 3 (`1/r⁴`) and 4 (`1/r²`)
integrate to zero.  The returned per-layer `β_k` is therefore
`a_k + b_k · F_k` with `F_k = (21/5) (3κ_k + μ_k)/μ_k
(r_k⁵ - r_{k-1}⁵)/(r_k³ - r_{k-1}³)`.

Reference: Hervé-Zaoui 1993, Christensen-Lo 1979, ECHOES C++
`inclusion_sphere_nlayers.h::get_visco_layer_average_strain_Strain`.
"""
function _shear_localization_multi(
        sphere::LayeredSphere{T, N}, C₀::TensND.TensISO{4, 3}
    ) where {T, N}
    states, _ = _shear_state_seq(sphere, C₀)
    κμ = _bulk_layer_moduli(sphere)
    radii = sphere.radii
    return ntuple(N) do k
        (κk, μk) = κμ[k]
        (a_k, b_k, _, _) = _shear_extract_amplitudes(radii[k], κk, μk, states[k])
        r_a = (k == 1) ? zero(eltype(radii)) : radii[k - 1]
        r_b = radii[k]
        a_k + b_k * _layer_avg_dev_shear_factor(r_a, r_b, κk, μk)
    end
end

"""
    _shear_localization(sphere, C₀) -> NTuple{N}

Per-layer deviatoric localisation `β_k` under a remote unit deviatoric
far-field.  Dispatches to the single-layer Eshelby delegation for
`N = 1` and to the state-vector recurrence for `N ≥ 2`.
"""
function _shear_localization(
        sphere::LayeredSphere{T, N}, C₀::TensND.TensISO{4, 3}
    ) where {T, N}
    N == 1 && return (_shear_localization_single_layer(sphere, C₀),)
    return _shear_localization_multi(sphere, C₀)
end
