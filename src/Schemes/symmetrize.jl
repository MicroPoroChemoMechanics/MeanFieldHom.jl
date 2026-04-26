# =============================================================================
#  symmetrize.jl — orientation-distribution projection of a phase's tensors.
#
#  `_apply_symmetrize(t, sym)` applies the projection `sym` to a tensor `t`,
#  returning a tensor of the appropriate symmetry class :
#
#  * `NoSymmetrize`           — passthrough.
#  * `IsoSymmetrize`          — Reynolds average over the rotation group :
#      4th-order : (J::N) J + ((K_proj::N) / 5) K_proj          → TensISO{4}
#      2nd-order : (1/3) tr(N) δ                                → TensISO{2}
#  * `TISymmetrize(axis)`     — Reynolds average over rotations about `axis` :
#      4th-order : projection onto the 6-dim Walpole basis (axis-aligned)
#                   ; major-symmetric tensors collapse to 5-dim → TensTI{4}
#      2nd-order : a·(δ - n⊗n) + b·(n⊗n)                        → TensTI{2}
#
#  Mathematical references : the iso-projection coefficients follow from the
#  invariants of the rotation group ; the TI projection follows from the
#  Walpole basis ((1933, refined in Hill 1965 and Walpole 1981)).
# =============================================================================

"""
    _apply_symmetrize(t::AbstractTens, sym::AbstractSymmetrize) -> AbstractTens

Project `t` onto the symmetry class declared by `sym`. See
[`AbstractSymmetrize`](@ref) for the projection semantics.

`NoSymmetrize` is a passthrough ; `IsoSymmetrize` returns `TensISO` ;
`TISymmetrize(axis)` returns `TensTI` for 4th-order tensors and 2nd-order
tensors. Implemented for tensor orders 4 (elasticity) and 2 (conductivity).
"""
_apply_symmetrize(t::TensND.AbstractTens, ::NoSymmetrize) = t

# =============================================================================
#  Reference-medium projection for the localization-tensor computation.
#
#  When a phase declares an orientation-distribution projection `sym`, the
#  reference medium passed to `hill_tensor` / `strain_strain_loc` can be
#  pre-projected onto the same symmetry class. This keeps the localization
#  tensor in an analytical-friendly symmetry class (iso → iso analytical,
#  TI → TI-coaxial analytical) and avoids triggering the general-anisotropic
#  residue / DECUHR path which currently does not support `ForwardDiff.Dual`
#  coefficients.
#
#  Mathematically, this is consistent with the orientation-averaging
#  semantics of `sym` : the inclusion family sees an effectively
#  symmetry-projected matrix. For matrices that are already in the target
#  symmetry class, `_project_matrix` is the identity.
# =============================================================================

"""
    _project_matrix(P₀::AbstractTens, sym::AbstractSymmetrize) -> AbstractTens

Project `P₀` onto a symmetry class compatible with the analytical
localization-tensor branches, before computing the localization tensor
of a phase with the given symmetrize.

- `NoSymmetrize` : passthrough.
- `IsoSymmetrize` : project to the iso (J, K_proj) basis. The inclusion's
  hill tensor in an iso matrix is always TI (with the inclusion's axis),
  for which an analytical branch exists.
- `TISymmetrize` : also project the matrix to **iso** rather than to its
  TI form. Reason : the inclusion family at polar angle θ ≠ 0 from the
  symmetrize axis is *not* coaxial with the matrix's TI axis, and the
  analytical TI-coaxial localization branch does not apply ; routing
  through the general anisotropic branch is not currently
  ForwardDiff-compatible. The iso projection gives an analytical branch
  for every inclusion orientation, and the result is exact at the iso
  fixed-point of the SC iteration (where C₀ converges to iso anyway).
  The phase contribution is still projected onto TI(axis) by the
  outgoing `_apply_symmetrize`, so the outer symmetry semantics are
  preserved.
"""
_project_matrix(P₀::TensND.AbstractTens, ::NoSymmetrize) = P₀
_project_matrix(P₀::TensND.AbstractTens, ::IsoSymmetrize) =
    _apply_symmetrize(P₀, IsoSymmetrize())
_project_matrix(P₀::TensND.AbstractTens, ::TISymmetrize) =
    _apply_symmetrize(P₀, IsoSymmetrize())

# ── 4th-order tensors : iso projection ───────────────────────────────────────

function _apply_symmetrize(t::TensND.AbstractTens{4, 3}, ::IsoSymmetrize)
    arr = TensND.get_array(t)
    α = sum(arr[i, i, j, j] for i in 1:3, j in 1:3) / 3
    full_trace = sum(arr[i, j, i, j] for i in 1:3, j in 1:3)
    β = (full_trace - α) / 5
    return TensND.TensISO{3}(α, β)
end

# ── 2nd-order tensors : iso projection (spherical part) ──────────────────────

function _apply_symmetrize(t::TensND.AbstractTens{2, 3}, ::IsoSymmetrize)
    arr = TensND.get_array(t)
    λ = (arr[1, 1] + arr[2, 2] + arr[3, 3]) / 3
    return TensND.TensISO{3}(λ)
end

# ── 4th-order tensors : TI projection around axis n ──────────────────────────
#
# Walpole basis with axis n (n unit vector) :
#   nₙ = n⊗n,  nT = δ - n⊗n
#   W₁ = nₙ ⊗ nₙ
#   W₂ = (nT ⊗ nT) / 2
#   W₃ = (nₙ ⊗ nT) / √2
#   W₄ = (nT ⊗ nₙ) / √2
#   W₅ = nT ⊠ˢ nT − (nT ⊗ nT)/2
#   W₆ = nT ⊠ˢ nₙ + nₙ ⊠ˢ nT
#
# These six tensors are orthonormal in the Frobenius inner product. The TI
# projection is the orthogonal projection onto their span. For a major-
# symmetric input (W₃::N = W₄::N), we collapse to 5 components.

function _apply_symmetrize(t::TensND.AbstractTens{4, 3}, sym::TISymmetrize)
    arr = TensND.get_array(t)
    Tarr = eltype(arr)
    n = ntuple(i -> convert(Tarr, sym.axis[i]), 3)
    # Build W_k arrays for k=1..6 in cartesian indices, double-dot with `t`,
    # then reassemble. We rely on direct index summation (no TensND ops) to
    # stay generic w.r.t. Dual eltypes.
    nₙ = ntuple(i -> ntuple(j -> n[i] * n[j], 3), 3)        # 3×3 nested tuple
    δ  = ntuple(i -> ntuple(j -> i == j ? one(Tarr) : zero(Tarr), 3), 3)
    nT = ntuple(i -> ntuple(j -> δ[i][j] - nₙ[i][j], 3), 3)

    # W as 4-arrays of the same shape as `arr`
    W = Array{Tarr, 4}[zeros(Tarr, 3, 3, 3, 3) for _ in 1:6]
    @inbounds for i in 1:3, j in 1:3, k in 1:3, l in 1:3
        # symmetric helper : (a⊠ˢb)_ijkl = (a_ik b_jl + a_il b_jk) / 2
        nn_nn = nₙ[i][j] * nₙ[k][l]
        nT_nT = nT[i][j] * nT[k][l]
        nn_nT = nₙ[i][j] * nT[k][l]
        nT_nn = nT[i][j] * nₙ[k][l]
        nT_box_nT = (nT[i][k] * nT[j][l] + nT[i][l] * nT[j][k]) / 2
        nT_box_nn = (nT[i][k] * nₙ[j][l] + nT[i][l] * nₙ[j][k]) / 2
        nn_box_nT = (nₙ[i][k] * nT[j][l] + nₙ[i][l] * nT[j][k]) / 2

        W[1][i, j, k, l] = nn_nn
        W[2][i, j, k, l] = nT_nT / 2
        W[3][i, j, k, l] = nn_nT / sqrt(Tarr(2))
        W[4][i, j, k, l] = nT_nn / sqrt(Tarr(2))
        W[5][i, j, k, l] = nT_box_nT - nT_nT / 2
        W[6][i, j, k, l] = nT_box_nn + nn_box_nT
    end

    # Inner products ℓ_k = W_k :: t (Frobenius double-dot)
    ℓ = ntuple(k -> sum(W[k][i, j, p, q] * arr[i, j, p, q]
                         for i in 1:3, j in 1:3, p in 1:3, q in 1:3),
                6)

    # Force major-symmetry by averaging ℓ₃ and ℓ₄ (matches the W₃, W₄
    # convention — for a major-sym input they're already equal up to noise).
    ℓ34 = (ℓ[3] + ℓ[4]) / 2
    # Normalisation: W₁..W₄ are orthonormal in the Frobenius inner product
    # (each has ‖W‖² = 1) so the projection coefficient is the inner
    # product itself. W₅ and W₆ have ‖W‖² = 2 in 3D, so the basis
    # decomposition coefficient is the inner product divided by 2.
    data = (ℓ[1], ℓ[2], ℓ34, ℓ[5] / 2, ℓ[6] / 2)
    return TensND.TensTI{4, Tarr, 5}(data, n)
end

# ── 2nd-order tensors : TI projection around axis n ──────────────────────────
# A · (δ - n⊗n) + B · (n⊗n)  with  A = (tr T - n·T·n) / 2,  B = n·T·n.

function _apply_symmetrize(t::TensND.AbstractTens{2, 3}, sym::TISymmetrize)
    arr = TensND.get_array(t)
    Tarr = eltype(arr)
    n = ntuple(i -> convert(Tarr, sym.axis[i]), 3)
    # B = n · T · n
    B = zero(Tarr)
    @inbounds for i in 1:3, j in 1:3
        B += n[i] * arr[i, j] * n[j]
    end
    trT = arr[1, 1] + arr[2, 2] + arr[3, 3]
    A = (trT - B) / 2
    # TensTI{2, T, 2}(data = (a, b), n)  —  a transverse, b axial
    return TensND.TensTI{2, Tarr, 2}((A, B), n)
end
