# =============================================================================
#  hill_order2_3d.jl — 2nd-order Hill tensor (3D, conductivity / diffusion).
# =============================================================================

"""
    _hill_order2_3d_iso(ell::Ellipsoid{3}, K₀) -> AbstractTens{2,3}

2nd-order Hill polarisation tensor of an ellipsoid in an isotropic
conductor ``\\mathbf K_0 = K\\,\\mathbf 1``:

```
P(A, K·1) = I^A / K ,
```

where ``\\mathbf I^{\\mathbf A}`` is the Newton-potential geometric
tensor ([`tens_IA`](@ref), [Willis 1977](@cite willis1977)).
"""
function _hill_order2_3d_iso(ell::Ellipsoid{3, Spherical}, K₀)
    T = promote_type(eltype(ell.semi_axes), eltype(K₀))
    k = K₀.data[1]
    IA = tens_IA(ell)
    return TensND.TensISO{3}(T(IA[1, 1]) / k)
end

function _hill_order2_3d_iso(ell::Ellipsoid{3}, K₀)
    T = promote_type(eltype(ell.semi_axes), eltype(K₀))
    k = K₀.data[1]
    IA = tens_IA(ell)
    P_arr = zeros(T, 3, 3)
    for i in 1:3, j in 1:3
        P_arr[i, j] = T(IA[i, j]) / k
    end
    return TensND.Tens(P_arr)
end

"""
    _hill_order2_3d_aniso(ell::Ellipsoid{3}, K₀) -> AbstractTens{2,3}

2nd-order Hill polarisation tensor of an ellipsoid in an arbitrarily
anisotropic conductor, via the closed-form square-root
change-of-variable of [Giraud & Gruescu 2019](@cite giraudMOM2019)
(equivalent derivation by Green's function in
[Barthélémy 2009](@cite barthelemyTIPM2009)):

```
P(A, K) = K⁻¹ᐟ² · I^(A·K⁻¹ᐟ²) · K⁻¹ᐟ² ,
```

where ``\\mathbf I^{\\mathbf A\\cdot\\mathbf K^{-1/2}}`` is the Newton
potential of the fictitious ellipsoid whose shape tensor is
``\\mathbf A\\cdot\\mathbf K^{-1/2}`` (semi-axes obtained by
diagonalising ``\\mathbf K^{-1/2}\\cdot\\mathbf A^{\\!T}\\!\\cdot
\\mathbf A\\cdot\\mathbf K^{-1/2}``).
"""
function _hill_order2_3d_aniso(ell::Ellipsoid{3}, K₀)
    T_mat = eltype(K₀)

    K_arr = Matrix{T_mat}(undef, 3, 3)
    for i in 1:3, j in 1:3
        K_arr[i, j] = K₀[i, j]
    end

    F = eigen(Symmetric(K_arr))
    invsqrt_K = F.vectors * Diagonal(1 ./ sqrt.(F.values)) * F.vectors'

    R_ell = [ell.basis[i, j] for i in 1:3, j in 1:3]
    A = R_ell * Diagonal(collect(ell.semi_axes)) * R_ell'

    F2 = svd(A * invsqrt_K)
    perm = sortperm(F2.S, rev = true)
    s = F2.S[perm]
    V = F2.V[:, perm]

    Iv, _ = MFH_Core.newton_potential_3d(s[1], s[2], s[3])

    P₀_arr = V * Diagonal([Iv[1], Iv[2], Iv[3]] ./ (4π)) * V'
    P_arr = invsqrt_K * P₀_arr * invsqrt_K

    return TensND.Tens(P_arr, TensND.CanonicalBasis{3, Float64}())
end

# ── Infinite cylinder (axis = e₁) ────────────────────────────────────────────

"""
    _hill_order2_3d_iso(cyl::Cylinder, K₀) -> AbstractTens{2,3}

2nd-order Hill polarisation tensor of an infinite cylinder (axis
``\\hat{\\mathbf e}_1``, transverse semi-axes ``b\\ge c>0``) in an
isotropic conductor ``\\mathbf K_0 = K\\,\\mathbf 1``, obtained from the
cylinder Newton-potential coefficients
([Mura 1987](@cite mura1987), §11.22):

```
P = I^cyl / K ,   with   I₁^cyl = 0,   I₂^cyl = c/(b+c),   I₃^cyl = b/(b+c) .
```

``P_{11} = 0`` expresses that no polarisation is transmitted along the
cylinder axis.
"""
function _hill_order2_3d_iso(cyl::Cylinder, K₀)
    T = promote_type(eltype(cyl.semi_axes), eltype(K₀))
    k = K₀.data[1]
    IA = tens_IA(cyl)
    P_arr = zeros(T, 3, 3)
    for i in 1:3, j in 1:3
        P_arr[i, j] = T(IA[i, j]) / k
    end
    return TensND.Tens(P_arr)
end

"""
    _hill_order2_3d_aniso(cyl::Cylinder, K₀) -> AbstractTens{2,3}

2nd-order Hill polarisation tensor of an infinite cylinder in an
arbitrarily anisotropic conductor.  The transverse plane
``(\\hat{\\mathbf e}_2,\\hat{\\mathbf e}_3)`` carries the full 2-D
Hill problem: the ``\\mathbf K^{-1/2}`` transformation of
[Giraud & Gruescu 2019](@cite giraudMOM2019) is applied to the
transverse 2×2 sub-matrix of ``\\mathbf K_0`` in the cylinder frame;
the 2-D Newton potentials produce the transverse block, and the axial
row/column is re-embedded as zero (``P_{1j}=0``).
"""
function _hill_order2_3d_aniso(cyl::Cylinder, K₀)
    T = promote_type(eltype(cyl.semi_axes), eltype(K₀))
    T_mat = eltype(K₀)

    K₀_princ = TensND.change_tens(K₀, cyl.basis)
    K_full = Matrix{T_mat}(undef, 3, 3)
    for i in 1:3, j in 1:3
        K_full[i, j] = K₀_princ[i, j]
    end

    K2 = Matrix{T_mat}(undef, 2, 2)
    for i in 1:2, j in 1:2
        K2[i, j] = K_full[i + 1, j + 1]
    end

    F = eigen(Symmetric(K2))
    invsqrt_K = F.vectors * Diagonal(1 ./ sqrt.(F.values)) * F.vectors'

    b, c = cyl.semi_axes
    A2 = Diagonal([b, c])
    F2 = svd(A2 * invsqrt_K)
    perm = sortperm(F2.S, rev = true)
    s = F2.S[perm]
    V = F2.V[:, perm]

    Iv2 = MFH_Core.newton_potential_2d(s[1], s[2])

    P_arr_2d_princ = V * Diagonal([Iv2[1], Iv2[2]] ./ (2π)) * V'
    P_arr_2d = invsqrt_K * P_arr_2d_princ * invsqrt_K

    P_arr_3d_princ = zeros(T, 3, 3)
    for i in 1:2, j in 1:2
        P_arr_3d_princ[i + 1, j + 1] = P_arr_2d[i, j]
    end

    return TensND.Tens(P_arr_3d_princ, cyl.basis)
end
