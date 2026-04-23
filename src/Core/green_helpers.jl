# =============================================================================
#  green_helpers.jl
#
#  Quadrature-agnostic building blocks of the crack-plane / half-sphere
#  Green-function evaluation.  Shared by every numerical back-end
#  (`green_nestedquadgk.jl`, `green_decuhr.jl`) in `Cracks/` and
#  `Elasticity/`.
#
#  `_A_and_Tn`          : φ-independent pre-computation (n̂-only).
#  `_phi_cache`         : φ-dependent, α-independent pre-computation.
#  `_qnn_pair_components!`: inner-α evaluation of
#           `[Q̂_{nn}(ζp) + Q̂_{nn}(ζm)] · ρ/sin²α`
#           from the pre-computed `(A, Vs, Ks, Kns, ca, sa)` quantities.
# =============================================================================

"""
    _A_and_Tn(C, n̂, ::Type{T}) -> (A, Tn)

Pre-compute the n̂-only quantities used by every `Q̂_{nn}` evaluation:

  * `Tn[i, p, q] = Σ_α C_{i α p q} n̂_α`          (3×3×3 array)
  * `A[i, p]    = Σ_q Tn[i, p, q] n̂_q`            (3×3 = V(n̂) = K(n̂))

Element type is the supplied `T = promote_type(...)`.
"""
function _A_and_Tn(C::AbstractArray, n̂::AbstractVector, ::Type{T}) where {T}
    Tn = Array{T}(undef, 3, 3, 3)
    A = Matrix{T}(undef, 3, 3)
    @inbounds for q in 1:3, p in 1:3, i in 1:3
        s = zero(T)
        for α in 1:3
            s += T(C[i, α, p, q]) * T(n̂[α])
        end
        Tn[i, p, q] = s
    end
    @inbounds for p in 1:3, i in 1:3
        s = zero(T)
        for q in 1:3
            s += Tn[i, p, q] * T(n̂[q])
        end
        A[i, p] = s
    end
    return A, Tn
end

"""
    _phi_cache(C, Tn, n̂, ξshat, ::Type{T}) -> (Vs, Ks, Kns)

Pre-compute the three 3×3 matrices that depend only on `ξshat` (the
in-plane unit direction):

  * `Vs[i, p]  = Σ_q Tn[i, p, q] ξshat_q                   = V(ξshat)`
  * `Ks[i, j]  = Σ_{k,l} C_{i k j l} ξshat_k ξshat_l        = K(ξshat)`
  * `Kns[i, j] = Σ_{k,l} C_{i k j l} (n̂_k ξshat_l + ξshat_k n̂_l)`
"""
function _phi_cache(
        C::AbstractArray, Tn::AbstractArray,
        n̂::AbstractVector, ξshat::AbstractVector,
        ::Type{T}
    ) where {T}
    Vs = Matrix{T}(undef, 3, 3)
    Ks = Matrix{T}(undef, 3, 3)
    Kns = Matrix{T}(undef, 3, 3)
    @inbounds for p in 1:3, i in 1:3
        s = zero(T)
        for q in 1:3
            s += Tn[i, p, q] * ξshat[q]
        end
        Vs[i, p] = s
    end
    @inbounds for j in 1:3, i in 1:3
        sk = zero(T); skn = zero(T)
        for k in 1:3, l in 1:3
            cc = T(C[i, k, j, l])
            sk += cc * ξshat[k] * ξshat[l]
            skn += cc * (n̂[k] * ξshat[l] + ξshat[k] * n̂[l])
        end
        Ks[i, j] = sk
        Kns[i, j] = skn
    end
    return Vs, Ks, Kns
end

"""
    _qnn_pair_components!(out, A, Vs, Ks, Kns, ca, sa, scale) -> out

Write `[Q̂_{nn}(ζp) + Q̂_{nn}(ζm)] · scale` into the 3×3 buffer `out`,
given the pre-computed φ-only quantities (`A, Vs, Ks, Kns`) and the
α-only trigs (`ca = cos α, sa = sin α`) plus a caller-supplied
`scale` factor that bundles the residual ρ / sin²α prefactor.

Operates in-place on `out`.
"""
@inline function _qnn_pair_components!(
        out::AbstractMatrix{T},
        A::AbstractMatrix{T},
        Vs::AbstractMatrix{T},
        Ks::AbstractMatrix{T},
        Kns::AbstractMatrix{T},
        ca::T, sa::T,
        scale::T
    ) where {T}
    cs = ca * sa
    ca² = ca * ca
    sa² = sa * sa
    Vp = ca .* A .+ sa .* Vs
    Vm = sa .* Vs .- ca .* A
    Kp = ca² .* A .+ cs .* Kns .+ sa² .* Ks
    Km = ca² .* A .- cs .* Kns .+ sa² .* Ks
    iKp = _inv3(Kp)
    iKm = _inv3(Km)
    Bp = (Vp * iKp) * transpose(Vp)
    Bm = (Vm * iKm) * transpose(Vm)
    @inbounds for j in 1:3, i in 1:3
        out[i, j] = (2 * A[i, j] - Bp[i, j] - Bm[i, j]) * scale
    end
    return out
end
