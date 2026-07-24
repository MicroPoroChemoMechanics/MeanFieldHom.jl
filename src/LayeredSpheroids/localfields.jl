# =============================================================================
#  localfields.jl — pointwise reconstruction of the temperature,
#  gradient and flux fields inside/around an `N`-layer confocal
#  spheroid, from the per-layer series coefficients
#  ([`spheroid_state_sequence`](@ref)).
#
#  Scope: fields are evaluated for the two CANONICAL remote loadings
#  the transfer-matrix recurrence directly solves — axial
#  (`H = H_axial·axis`) and transverse (`H = H_trans·ê₁`, `ê₁ ⟂ axis`
#  at `φ = 0`) — and their linear superposition (exact: the axial/
#  transverse harmonic subspaces are decoupled, eq:Tsol). A general
#  in-plane direction follows by evaluating at a shifted `φ`
#  (axisymmetry: `Hᵀ = H₁ê₁ + H₂ê₂` at azimuth `φ` behaves as a single
#  transverse magnitude `√(H₁²+H₂²)` at azimuth `φ - atan2(H₂,H₁)`).
#
#  Coordinates `(q, p, φ)` are in the spheroid's OWN principal frame
#  (revolution axis ≡ local `ê₃`); results are returned in that same
#  local frame. Rotate by `s.axis` (TensND) if the caller needs the
#  global frame and `s.axis ≠ (0,0,1)`.
# =============================================================================

"""
    get_layer(s::LayeredSpheroid, q) -> Int

Index of the layer whose confocal range `(q_{k-1}, q_k]` contains `q`
(`1, …, N`), or `N+1` if `q` lies in the surrounding matrix
(`|q| > |q_N|`).
"""
function get_layer(s::LayeredSpheroid{T, N}, q) where {T, N}
    k = 1
    while k ≤ N && abs(q) ≥ abs(s.q[k])
        k += 1
    end
    return k
end

@inline function _base_fond(q::Qx, p, φ) where {Qx}
    qb = sqrt(q^2 - 1)
    pb = sqrt(1 - p^2)
    qp = sqrt(q^2 - p^2)
    eR = (cos(φ), sin(φ), zero(φ))
    ephi = (-sin(φ), cos(φ), zero(φ))
    ez = (zero(φ), zero(φ), one(φ))
    e_p = ((-p * qb) .* eR .+ (q * pb) .* ez) ./ qp
    e_q = ((pb * q) .* eR .+ (p * qb) .* ez) ./ qp
    return e_q, e_p, ephi
end

@inline function _metric(c::Qx, q, p) where {Qx}
    qb = sqrt(q^2 - 1)
    pb = sqrt(1 - p^2)
    qp = sqrt(q^2 - p^2)
    return c * qp / qb, c * qp / pb, c * qb * pb   # h_q, h_p, h_φ
end

"""
    _spheroid_field_coeffs(s, k₀, q, p; N_only = nothing)
        -> (layer, k_layer, Aa, Ba, At, Bt)

The layer containing `q`, its conductivity (matrix `k₀` if `q` lies
outside the spheroid), and the axial/transverse coefficient
sub-vectors valid there.
"""
function _spheroid_field_coeffs(s::LayeredSpheroid{T, N}, k₀) where {T, N}
    Xa = spheroid_state_sequence(s, k₀, false)
    Xt = spheroid_state_sequence(s, k₀, true)
    return Xa, Xt
end

"""
    local_temperature(s, k₀, q, p, φ; H_axial = 1.0, H_trans = 0.0) -> T

Temperature at the spheroidal point `(q, p, φ)` (own frame) under a
remote gradient `H_axial·axis + H_trans·ê₁` (superposition of the
axial and transverse canonical problems, eq:Taxi/eq:Ttrans).
"""
function local_temperature(
        s::LayeredSpheroid{T, N}, k₀, q, p, φ; H_axial = 1.0, H_trans = 0.0,
    ) where {T, N}
    𝒩 = s.Nseries
    lay = get_layer(s, q)
    Xa, Xt = _spheroid_field_coeffs(s, k₀)
    Aa, Ba = Xa[lay][1:𝒩], Xa[lay][(𝒩 + 1):(2𝒩)]
    At, Bt = Xt[lay][1:𝒩], Xt[lay][(𝒩 + 1):(2𝒩)]

    P0p, _ = legendre_odd(:P0, p, 𝒩)
    P0q, _ = legendre_odd(:P0, q, 𝒩)
    Q0q, _ = legendre_odd(:Q0, q, 𝒩)
    P1p, _ = legendre_odd(:P1p, p, 𝒩)
    P1q, _ = legendre_odd(:P1, q, 𝒩)
    Q1q, _ = legendre_odd(:Q1, q, 𝒩)

    Ta = s.c * sum(P0p[r] * (Aa[r] * P0q[r] + Ba[r] * Q0q[r]) for r in 1:𝒩)
    Tt = s.c * sum(P1p[r] * (At[r] * P1q[r] + Bt[r] * Q1q[r]) for r in 1:𝒩)
    return real(H_axial * Ta + H_trans * Tt * cos(φ))
end

"""
    local_gradient(s, k₀, q, p, φ; H_axial = 1.0, H_trans = 0.0) -> (g₁, g₂, g₃)

Temperature gradient `∇T` at `(q, p, φ)` (own frame, real Cartesian
triple), under the same remote loading as [`local_temperature`](@ref).
"""
function local_gradient(
        s::LayeredSpheroid{T, N}, k₀, q, p, φ; H_axial = 1.0, H_trans = 0.0,
    ) where {T, N}
    𝒩 = s.Nseries
    lay = get_layer(s, q)
    Xa, Xt = _spheroid_field_coeffs(s, k₀)
    Aa, Ba = Xa[lay][1:𝒩], Xa[lay][(𝒩 + 1):(2𝒩)]
    At, Bt = Xt[lay][1:𝒩], Xt[lay][(𝒩 + 1):(2𝒩)]

    P0p, dP0p = legendre_odd(:P0, p, 𝒩)
    P0q, dP0q = legendre_odd(:P0, q, 𝒩)
    Q0q, dQ0q = legendre_odd(:Q0, q, 𝒩)
    P1p, dP1p = legendre_odd(:P1p, p, 𝒩)
    P1q, dP1q = legendre_odd(:P1, q, 𝒩)
    Q1q, dQ1q = legendre_odd(:Q1, q, 𝒩)

    dTa_dq = s.c * sum(P0p[r] * (Aa[r] * dP0q[r] + Ba[r] * dQ0q[r]) for r in 1:𝒩)
    dTa_dp = s.c * sum(dP0p[r] * (Aa[r] * P0q[r] + Ba[r] * Q0q[r]) for r in 1:𝒩)
    Tt = s.c * sum(P1p[r] * (At[r] * P1q[r] + Bt[r] * Q1q[r]) for r in 1:𝒩)
    dTt_dq = s.c * sum(P1p[r] * (At[r] * dP1q[r] + Bt[r] * dQ1q[r]) for r in 1:𝒩)
    dTt_dp = s.c * sum(dP1p[r] * (At[r] * P1q[r] + Bt[r] * Q1q[r]) for r in 1:𝒩)

    e_q, e_p, e_φ = _base_fond(q, p, φ)
    h_q, h_p, h_φ = _metric(s.c, q, p)

    g_axial = (dTa_dq / h_q) .* e_q .+ (dTa_dp / h_p) .* e_p
    g_trans_r = (dTt_dq / h_q) .* e_q .+ (dTt_dp / h_p) .* e_p   # radial-type part
    g_trans = cos(φ) .* g_trans_r .- (sin(φ) * Tt / h_φ) .* e_φ

    g = H_axial .* g_axial .+ H_trans .* g_trans
    return real(g[1]), real(g[2]), real(g[3])
end

"""
    local_flux(s, k₀, q, p, φ; H_axial = 1.0, H_trans = 0.0) -> (u₁, u₂, u₃)

Heat/mass flux `u = -k(x)·∇T` at `(q, p, φ)` (own frame), `k(x)` the
conductivity of the layer containing the point (or `k₀` outside the
spheroid).
"""
function local_flux(
        s::LayeredSpheroid{T, N}, k₀, q, p, φ; H_axial = 1.0, H_trans = 0.0,
    ) where {T, N}
    lay = get_layer(s, q)
    k_layers = _spheroid_layer_moduli(s)
    k_here = lay ≤ N ? k_layers[lay] : _as_scalar_k(k₀)
    g1, g2, g3 = local_gradient(s, k₀, q, p, φ; H_axial, H_trans)
    return -k_here * g1, -k_here * g2, -k_here * g3
end
