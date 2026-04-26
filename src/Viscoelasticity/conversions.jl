# =============================================================================
#  conversions.jl — convert between a `(6n×6n)` block matrix and per-symmetry
#  scalar parameter matrices `n×n`.
#
#  An isotropic 4-tensor in Mandel form is described by two scalars
#  `(α, β)` (with `α = 3K`, `β = 2μ`).  An iso block matrix in ALV is
#  fully described by the two scalar Volterra matrices `α(t,t')` and
#  `β(t,t')` of size `n×n`.
#
#  Conversions are needed at two places :
#   * extracting the iso parameters from a 6n×6n viscoelastic block
#     matrix in order to apply the time-space decoupling formula
#     for the Hill kernel ;
#   * reassembling a 6n×6n matrix after applying a scalar Volterra
#     operation (e.g. a Volterra inverse) on each component.
# =============================================================================

"""
    iso_params_from_blocks(M) -> (α::Matrix, β::Matrix)

Decompose a `6n×6n` block matrix whose every 6×6 block is an
**isotropic** 4-tensor in Mandel form into the two scalar parameter
matrices `α` and `β`, both of size `n × n`.  The block layout is

```
  M_block(i,j) = α[i,j] · 𝕁_Mandel + β[i,j] · 𝕂_Mandel
```

where `𝕁_Mandel = (1/3) e₁ e₁ᵀ` (`e₁` the Mandel-1 unit) augmented
with the symmetric off-diagonal part `(1/3)` for the upper 3×3 block,
and `𝕂_Mandel = 𝕀_Mandel - 𝕁_Mandel`.

`α[i,j]` and `β[i,j]` are extracted from the diagonal entry `M_block[1,1]`
and the (4,4) Mandel-shear entry of each block:
   `α[i,j] = (M_block[1,1] + 2 M_block[1,2])`
   `β[i,j] = M_block[4,4]`.
"""
function iso_params_from_blocks(M::AbstractMatrix)
    sz = size(M, 1)
    sz == size(M, 2) || throw(ArgumentError("iso_params_from_blocks: M must be square"))
    sz % 6 == 0 || throw(ArgumentError("iso_params_from_blocks: size $(sz) not divisible by 6"))
    n = sz ÷ 6
    T = eltype(M)
    α = zeros(T, n, n)
    β = zeros(T, n, n)
    @inbounds for i in 1:n, j in 1:n
        r1 = 6 * (i - 1) + 1
        c1 = 6 * (j - 1) + 1
        c2 = 6 * (j - 1) + 2
        # α = M[1,1] + 2 M[1,2] (= 3K), β = M[4,4] (= 2μ in Mandel).
        α[i, j] = M[r1, c1] + 2 * M[r1, c2]
        β[i, j] = M[r1 + 3, c1 + 3]
    end
    return α, β
end

"""
    iso_blocks_from_params(α::Matrix, β::Matrix) -> Matrix

Inverse of [`iso_params_from_blocks`](@ref): build a `6n×6n` block
matrix whose every block is the iso 4-tensor `α[i,j] · 𝕁 + β[i,j] · 𝕂`
in Mandel form.

Both `α` and `β` must be `n × n` and have a common element type.
"""
function iso_blocks_from_params(α::AbstractMatrix, β::AbstractMatrix)
    size(α) == size(β) ||
        throw(ArgumentError("iso_blocks_from_params: α and β must have same size"))
    n = size(α, 1)
    n == size(α, 2) ||
        throw(ArgumentError("iso_blocks_from_params: α must be square"))
    T = promote_type(eltype(α), eltype(β))
    M = zeros(T, 6 * n, 6 * n)
    @inbounds for i in 1:n, j in 1:n
        a = T(α[i, j])
        b = T(β[i, j])
        rows = (6 * (i - 1) + 1):(6 * i)
        cols = (6 * (j - 1) + 1):(6 * j)
        # Mandel iso 4-tensor (α, β):
        #   diagonal 3×3 block:   diag = (α + 2β)/3,  off-diag = (α − β)/3
        #   diagonal shear block: 2μ = β on the diagonal of the 4..6 sub-block
        diag_top = (a + 2b) / 3
        offdiag_top = (a - b) / 3
        for k in 1:3
            for l in 1:3
                if k == l
                    M[rows[k], cols[l]] = diag_top
                else
                    M[rows[k], cols[l]] = offdiag_top
                end
            end
        end
        for k in 4:6
            M[rows[k], cols[k]] = b
        end
    end
    return M
end
