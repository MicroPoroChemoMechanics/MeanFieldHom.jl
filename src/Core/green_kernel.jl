# =============================================================================
#  green_kernel.jl
#
#    * `_inv3(K)` : explicit cofactor-based 3×3 inverse
#                   (ForwardDiff-safe, avoids LU factorization)
#
#  Shared by the residue and DECUHR paths of the `Elasticity` and `Cracks`
#  sub-modules.
#
#  This file used to also carry `_acoustic_tensor(C, ξ)` and a direct
#  pointwise `Q̂_{nn}` evaluation, `_Qnn_direct(C, ξ, n̂)`.  Both were dead
#  code, and `_Qnn_direct` contracted the stiffness with the Green kernel
#  through a six-deep `p,q,r,s,α,β` loop — 729 iterations where the
#  contraction factorizes into `U = (C·n̂)·ξ` then `B = U·K⁻¹·Uᵀ`, roughly
#  100 flops.  The factorized form is what the live code already uses, in
#  two places:
#
#    * `Cracks/green_residue.jl`  (`Tncon → V → M·Vᵀ`, polynomial path)
#    * `Core/green_helpers.jl`    (`_qnn_pair_components!`, quadrature path)
#
#  Look there rather than reinstating a third, slower copy.
# =============================================================================

"""
    _inv3(K) -> SMatrix{3,3,T}

Explicit closed-form inverse of a 3×3 matrix via the cofactor formula.
Avoids the overhead of `inv`/LU factorization on tiny matrices and is
fully ForwardDiff-compatible (uses only `+`, `-`, `*`, `/`).

Returns a stack-allocated `SMatrix`: this sits in the innermost loop of every
crack COD back-end, where the previous heap `Matrix{T}` was one of ~10
allocations per quadrature node.
"""
@inline function _inv3(K::AbstractMatrix{T}) where {T}
    return @inbounds begin
        a11 = K[1, 1]; a12 = K[1, 2]; a13 = K[1, 3]
        a21 = K[2, 1]; a22 = K[2, 2]; a23 = K[2, 3]
        a31 = K[3, 1]; a32 = K[3, 2]; a33 = K[3, 3]
        c11 = a22 * a33 - a23 * a32
        c12 = a13 * a32 - a12 * a33
        c13 = a12 * a23 - a13 * a22
        c21 = a23 * a31 - a21 * a33
        c22 = a11 * a33 - a13 * a31
        c23 = a13 * a21 - a11 * a23
        c31 = a21 * a32 - a22 * a31
        c32 = a12 * a31 - a11 * a32
        c33 = a11 * a22 - a12 * a21
        det = a11 * c11 + a12 * c21 + a13 * c31
        invd = inv(det)
        # column-major: (11,21,31, 12,22,32, 13,23,33)
        SMatrix{3, 3, T}(
            c11 * invd, c21 * invd, c31 * invd,
            c12 * invd, c22 * invd, c32 * invd,
            c13 * invd, c23 * invd, c33 * invd
        )
    end
end
