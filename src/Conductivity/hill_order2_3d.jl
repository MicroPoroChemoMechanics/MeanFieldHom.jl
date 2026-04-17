# =============================================================================
#  hill_order2_3d.jl — 2nd-order Hill tensor (3D, conductivity / diffusion).
# =============================================================================

"""
    _hill_order2_3d_iso(ell::Ellipsoid{3}, K₀) -> AbstractTens{2,3}

2nd-order Hill tensor `P` for an isotropic conductor `K₀ = k·δ`.
`P(A, k·δ) = I^A / k`.
"""
function _hill_order2_3d_iso(ell::Ellipsoid{3, Spherical}, K₀)
    T = promote_type(eltype(ell.semi_axes), eltype(K₀))
    k = K₀.data[1]
    IA = tens_IA(ell)
    return TensND.TensISO{3}(T(IA[1,1]) / k)
end

function _hill_order2_3d_iso(ell::Ellipsoid{3}, K₀)
    T = promote_type(eltype(ell.semi_axes), eltype(K₀))
    k = K₀.data[1]
    IA = tens_IA(ell)
    P_arr = zeros(T, 3, 3)
    for i in 1:3, j in 1:3
        P_arr[i,j] = T(IA[i,j]) / k
    end
    return TensND.Tens(P_arr)
end

"""
    _hill_order2_3d_aniso(ell::Ellipsoid{3}, K₀) -> AbstractTens{2,3}

2nd-order Hill tensor for a general anisotropic conductor via the
K^{-1/2} change-of-variable.
"""
function _hill_order2_3d_aniso(ell::Ellipsoid{3}, K₀)
    T_mat = eltype(K₀)

    K_arr = Matrix{T_mat}(undef, 3, 3)
    for i in 1:3, j in 1:3
        K_arr[i,j] = K₀[i,j]
    end

    F         = eigen(Symmetric(K_arr))
    invsqrt_K = F.vectors * Diagonal(1 ./ sqrt.(F.values)) * F.vectors'

    R_ell = [ell.basis[i,j] for i in 1:3, j in 1:3]
    A     = R_ell * Diagonal(collect(ell.semi_axes)) * R_ell'

    F2   = svd(A * invsqrt_K)
    perm = sortperm(F2.S, rev=true)
    s    = F2.S[perm]
    U    = F2.U[:, perm]

    Iv, _ = MFH_Core.newton_potential_3d(s[1], s[2], s[3])

    P₀_arr = U * Diagonal([Iv[1], Iv[2], Iv[3]] ./ (4π)) * U'
    P_arr  = invsqrt_K * P₀_arr * invsqrt_K

    return TensND.change_tens_canon(TensND.Tens(P_arr, TensND.CanonicalBasis{3,Float64}()))
end
