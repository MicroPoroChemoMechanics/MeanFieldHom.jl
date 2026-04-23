# =============================================================================
#  conductivity.jl — isotropic `n`-layer spherical inclusion, 2nd-order
#  transport problem (thermal / electric / Darcy).
#
#  Under a remote uniform gradient `∇T∞` (Y₁-harmonic direction field),
#  the temperature field inside a concentric-sphere laminate takes the
#  form
#
#      T(r, θ) = [A_k · r + B_k / r²] · (r̂ · ê∇),
#
#  where ê∇ is the direction of `∇T∞`.  The state vector at radius r is
#  `s(r) = (T̂(r), q̂_n(r))` with
#      T̂ = A r + B/r²                    (temperature amplitude),
#      q̂_n = -k (A - 2 B/r³)             (normal flux amplitude).
#
#  Intra-layer transfer  `s(r_out) = T_cond(r_out, r_in; k) · s(r_in)`
#  comes from `M(r_out) · M(r_in)⁻¹` with `M(r) = [r  1/r²; -k  2k/r³]`
#  and `det M(r) = 3k/r²`.  Substituting yields the closed form
#
#      T_cond = (1/3) · [ 2 ro/ri + (ri/ro)²            (ri³/ro² − ro)/k
#                         -2k (1/ri − ri²/ro³)           ro³/ri³ (wait…) ]
#
#  which expands to the four factored coefficients below (all finite
#  for any `k > 0`).
#
#  Interface jumps:
#     Perfect                       J = I
#     Kapitza(ρ)                    J = [1  ρ ; 0  1]
#     SurfaceConductive(ks)         J = [1  0 ; -n(n+1) ks/r²  1]  with n=1 ⇒ 2ks/r².
#
#  Per-layer gradient localisation `α_k = A_k / A_∞`.
# =============================================================================

"""
    _cond_layer_transfer(r_out, r_in, k_iso) -> Matrix(2×2)

Intra-layer field transfer for the Y₁-harmonic conductivity problem
between radii `r_in` and `r_out` in a layer of scalar conductivity
`k_iso`.  The closed form below shares `r²`, `r³`, and their
reciprocals; only one division by `k_iso` appears.
"""
@inline function _cond_layer_transfer(r_out, r_in, k_iso)
    T  = promote_type(typeof(r_out), typeof(r_in), typeof(k_iso))
    ro = T(r_out); ri = T(r_in); k = T(k_iso)

    inv_ri  = 1 / ri
    inv_ro² = 1 / (ro * ro)
    inv_ro³ = inv_ro² / ro
    ri²     = ri * ri
    ri³     = ri² * ri
    inv_k   = 1 / k

    # M(ri)⁻¹ = (ri²/(3k)) · [2k/ri³  -1/ri²; k  ri]
    #        = [2/(3 ri)  -1/(3k);  k ri²/(3k)  ri³/(3k)]
    # Let a = 2/(3ri), b = -1/(3k), c = ri²/3, d = ri³/(3k).
    # T_cond = M(ro) · M(ri)⁻¹ with M(ro) = [ro  1/ro²; -k  2k/ro³].
    a = 2 * inv_ri / 3
    b = -inv_k / 3
    c = ri² / 3
    d = ri³ * inv_k / 3

    T11 = ro * a + inv_ro² * c
    T12 = ro * b + inv_ro² * d
    T21 = -k * a + 2 * k * inv_ro³ * c
    T22 = -k * b + 2 * k * inv_ro³ * d

    return T[T11 T12; T21 T22]
end

"""
    _cond_seed_state(r_1, k_1) -> Vector(2)

Seed state at `r_1⁻` for a unit amplitude `A_1 = 1` in layer 1
(core).  `B_1 = 0` from finiteness at the origin, so `T̂(r_1) = r_1`,
`q̂_n(r_1) = -k_1`.
"""
@inline function _cond_seed_state(r_1, k_1)
    T = promote_type(typeof(r_1), typeof(k_1))
    return T[T(r_1), -T(k_1)]
end

"""
    _cond_extract_AB(r, k_iso, T̂, q̂n) -> (A, B)

Extract the local-expansion coefficients `(A, B)` in a layer of
conductivity `k_iso` given the state `(T̂, q̂_n)` at radius `r`.
"""
@inline function _cond_extract_AB(r, k_iso, T̂, q̂n)
    T   = promote_type(typeof(r), typeof(k_iso), typeof(T̂), typeof(q̂n))
    Tr  = T(r)
    Tr² = Tr * Tr
    Tr³ = Tr² * Tr
    Tk  = T(k_iso)
    inv_3 = one(T) / 3
    A = (2 * T(T̂) - Tr * T(q̂n) / Tk) * inv_3 / Tr
    B = (Tr² * T(T̂) + Tr³ * T(q̂n) / Tk) * inv_3
    return A, B
end

# ── Conductivity interface jump matrices ────────────────────────────────────

"""
    _cond_interface_T(intf, k_iso, r) -> Matrix(2×2)

Jump matrix for the Y₁-harmonic conductivity state vector
`(T̂, q̂_n)` at an interface of type `intf` located at radius `r`.
"""
function _cond_interface_T end

function _cond_interface_T(::PerfectInterface, k_iso, r)
    T = promote_type(typeof(k_iso), typeof(r))
    return T[one(T) zero(T); zero(T) one(T)]
end

function _cond_interface_T(intf::KapitzaInterface, k_iso, r)
    T = promote_type(eltype(intf), typeof(k_iso), typeof(r))
    return T[one(T) T(intf.resistance); zero(T) one(T)]
end

function _cond_interface_T(intf::SurfaceConductiveInterface, k_iso, r)
    T   = promote_type(eltype(intf), typeof(k_iso), typeof(r))
    Tr² = T(r) * T(r)
    return T[one(T) zero(T); -2 * T(intf.conductance) / Tr² one(T)]
end

# ── Propagation and localisation ────────────────────────────────────────────

# Cached tuple of per-layer scalar conductivities, built once from the
# `sphere.moduli` field.
@inline function _cond_layer_moduli(sphere::LayeredSphere{T, N}) where {T, N}
    return ntuple(k -> _iso_scalar(layer_modulus(sphere, k)), Val(N))
end

"""
    _cond_state_seq(sphere, k₀) -> NTuple{N, state⁻}, state⁺_N

Propagate the conductivity state vector from the core outward.
"""
function _cond_state_seq(sphere::LayeredSphere{T, N}, k₀) where {T, N}
    k_layers = _cond_layer_moduli(sphere)
    TP = promote_type(T, typeof(k₀),
                      ntuple(k -> typeof(k_layers[k]), N)...)
    radii = sphere.radii

    s = _cond_seed_state(TP(radii[1]), TP(k_layers[1]))

    inside_states = Vector{Vector{TP}}(undef, N)
    for k in 1:N
        inside_states[k] = copy(s)
        intf = layer_interface(sphere, k)
        Tint = _cond_interface_T(intf, TP(k_layers[k]), TP(radii[k]))
        s = Tint * s
        if k < N
            Tlay = _cond_layer_transfer(TP(radii[k + 1]), TP(radii[k]),
                                        TP(k_layers[k + 1]))
            s = Tlay * s
        end
    end
    return inside_states, s
end

"""
    _cond_localization(sphere, k₀) -> NTuple{N, TP}

Per-layer gradient localisation `α_k = A_k / A_∞` under a remote
uniform gradient.  Reduces, for `N = 1`, to the classical formula
`α_1 = 3 k_0 / (2 k_0 + k_1)` for a sphere inclusion.
"""
function _cond_localization(sphere::LayeredSphere{T, N}, k₀) where {T, N}
    k_layers = _cond_layer_moduli(sphere)
    TP = promote_type(T, typeof(k₀),
                      ntuple(k -> typeof(k_layers[k]), N)...)
    inside_states, s_matrix = _cond_state_seq(sphere, k₀)
    radii = sphere.radii

    A_inf, _ = _cond_extract_AB(TP(radii[N]), TP(k₀), s_matrix[1], s_matrix[2])
    inv_A_inf = one(TP) / A_inf

    return ntuple(N) do k
        (A_k, _) = _cond_extract_AB(TP(radii[k]), TP(k_layers[k]),
                                    inside_states[k][1], inside_states[k][2])
        A_k * inv_A_inf
    end
end

"""
    _effective_conductivity(sphere, k₀) -> k_eff

Effective conductivity of the composite sphere:
`k_eff = Σ_k f_k k_k α_k`.
"""
function _effective_conductivity(sphere::LayeredSphere{T, N}, k₀) where {T, N}
    α = _cond_localization(sphere, k₀)
    k_layers = _cond_layer_moduli(sphere)
    f = ntuple(k -> layer_volume_fraction(sphere, k), Val(N))
    return sum(f[k] * k_layers[k] * α[k] for k in 1:N)
end
