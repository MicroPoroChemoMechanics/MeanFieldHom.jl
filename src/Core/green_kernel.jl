# =============================================================================
#  green_kernel.jl
#
#  Elementary building blocks of the Fourier–space Green operator:
#
#    * `_acoustic_tensor(C, ξ)` : 3×3 acoustic tensor `Kᵢⱼ = Cᵢₖⱼₗ ξₖ ξₗ`
#    * `_Qnn_direct(C, ξ, n̂)`   : direct pointwise `Q̂_{nn}` evaluation
#    * `_inv3(K)`               : explicit cofactor-based 3×3 inverse
#                                 (ForwardDiff-safe, avoids LU factorisation)
#
#  Shared by the residue and DECUHR paths of the `Elasticity` and `Cracks`
#  sub-modules.
# =============================================================================

"""
    _acoustic_tensor(C, ξ) -> 3×3 matrix

3×3 acoustic tensor ``K_{ij} = C_{ikjl}\\,ξ_k\\,ξ_l`` for a 4th-order
stiffness `C` (component array) and a wave vector `ξ`.  Symmetric when
`C` has the minor / major symmetries.  Element type follows
`promote_type(eltype(C), eltype(ξ))`.
"""
@inline function _acoustic_tensor(C::AbstractArray{TC,4},
                                  ξ::AbstractVector{Tξ}) where {TC<:Number, Tξ<:Number}
    T = promote_type(TC, Tξ)
    K = zeros(T, 3, 3)
    @inbounds for i in 1:3, j in 1:3
        s = zero(T)
        for k in 1:3, l in 1:3
            s += C[i, k, j, l] * ξ[k] * ξ[l]
        end
        K[i, j] = s
    end
    return K
end

"""
    _inv3(K) -> Matrix{T}

Explicit closed-form inverse of a 3×3 matrix via the cofactor formula.
Avoids the overhead of `inv`/LU factorisation on tiny matrices and is
fully ForwardDiff-compatible (uses only `+`, `-`, `*`, `/`).
"""
@inline function _inv3(K::AbstractMatrix{T}) where {T}
    @inbounds begin
        a11 = K[1,1]; a12 = K[1,2]; a13 = K[1,3]
        a21 = K[2,1]; a22 = K[2,2]; a23 = K[2,3]
        a31 = K[3,1]; a32 = K[3,2]; a33 = K[3,3]
        c11 = a22*a33 - a23*a32
        c12 = a13*a32 - a12*a33
        c13 = a12*a23 - a13*a22
        c21 = a23*a31 - a21*a33
        c22 = a11*a33 - a13*a31
        c23 = a13*a21 - a11*a23
        c31 = a21*a32 - a22*a31
        c32 = a12*a31 - a11*a32
        c33 = a11*a22 - a12*a21
        det = a11*c11 + a12*c21 + a13*c31
        invd = inv(det)
        iK = Matrix{T}(undef, 3, 3)
        iK[1,1] = c11*invd; iK[1,2] = c12*invd; iK[1,3] = c13*invd
        iK[2,1] = c21*invd; iK[2,2] = c22*invd; iK[2,3] = c23*invd
        iK[3,1] = c31*invd; iK[3,2] = c32*invd; iK[3,3] = c33*invd
        iK
    end
end

"""
    _Qnn_direct(C, ξ, n̂) -> 3×3 matrix

Direct pointwise evaluation of the normal–normal projection
``\\hat Q_{nn}(\\vec ξ)`` of the Fourier Green operator for a
stiffness `C` (component array) and a wave vector `ξ`.

Element type is `promote_type(eltype(C), eltype(ξ), eltype(n̂))` —
ForwardDiff- and SymPy-compatible.  For `‖ξ‖ = 0` the acoustic tensor
becomes singular; the caller must handle this.
"""
function _Qnn_direct(C::AbstractArray{TC,4},
                     ξ::AbstractVector{Tξ},
                     n̂::AbstractVector{Tn}) where {TC<:Number, Tξ<:Number, Tn<:Number}
    T = promote_type(TC, Tξ, Tn)
    K    = _acoustic_tensor(C, ξ)
    Kinv = _inv3(K)
    A = zeros(T, 3, 3)
    @inbounds for i in 1:3, k in 1:3
        s = zero(T)
        for α in 1:3, β in 1:3
            s += C[i, α, k, β] * n̂[α] * n̂[β]
        end
        A[i, k] = s
    end
    B = zeros(T, 3, 3)
    @inbounds for i in 1:3, k in 1:3
        s = zero(T)
        for p in 1:3, q in 1:3, r in 1:3, s2 in 1:3, α in 1:3, β in 1:3
            Γ = (ξ[q] * Kinv[p, r] * ξ[s2] +
                 ξ[q] * Kinv[p, s2] * ξ[r] +
                 ξ[p] * Kinv[q, r] * ξ[s2] +
                 ξ[p] * Kinv[q, s2] * ξ[r]) / 4
            s += C[i, α, p, q] * Γ * C[r, s2, k, β] * n̂[α] * n̂[β]
        end
        B[i, k] = s
    end
    return A - B
end
