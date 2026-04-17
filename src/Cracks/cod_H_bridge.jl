# =============================================================================
#  cod_H_bridge.jl — bridge between compliance tensor H and COD tensor B.
# =============================================================================

"""
    cod_from_compliance(H, ℬ=getbasis(H)) -> Tens{2,3}

Extract the size-independent COD tensor ``\\mathbf B`` from a crack
compliance contribution tensor ``\\mathbb H`` expressed in an arbitrary
frame.
"""
function cod_from_compliance(H, ℬ = TensND.getbasis(H))
    T = eltype(H)
    newH = TensND.change_tens(H, ℬ)
    return TensND.Tens(
        Tensors.SymmetricTensor{2, 3}(
            (i, j) ->
            16 * newH[i, 3, j, 3] /
                (
                3 * (
                    one(T) + MFH_Core._δ(i, 3, T) + MFH_Core._δ(j, 3, T) +
                        MFH_Core._δ(i, 3, T) * MFH_Core._δ(j, 3, T)
                )
            )
        ),
        ℬ,
    )
end

const BfromH = cod_from_compliance

"""
    compliance_from_cod(B, ℬ=getbasis(B)) -> Tens{4,3}

Inverse of [`cod_from_compliance`](@ref).  Reconstruct ``\\mathbb H``
from ``\\mathbf B``.
"""
function compliance_from_cod(B, ℬ = TensND.getbasis(B))
    T = eltype(B)
    newB = TensND.change_tens(B, ℬ)
    data = zeros(T, 3, 3, 3, 3)
    @inbounds for i in 1:3, k in 1:3
        j = 3
        l = 3
        v = newB[i, k] * 3 *
            (
            one(T) + MFH_Core._δ(i, 3, T) + MFH_Core._δ(k, 3, T) +
                MFH_Core._δ(i, 3, T) * MFH_Core._δ(k, 3, T)
        ) /
            16
        data[i, j, k, l] = v
    end
    sym = zeros(T, 3, 3, 3, 3)
    @inbounds for i in 1:3, j in 1:3, k in 1:3, l in 1:3
        sym[i, j, k, l] = (
            data[i, j, k, l] + data[j, i, k, l] +
                data[i, j, l, k] + data[j, i, l, k] +
                data[k, l, i, j] + data[l, k, i, j] +
                data[k, l, j, i] + data[l, k, j, i]
        ) / 8
    end
    return TensND.Tens(sym, ℬ)
end

const HfromB = compliance_from_cod

"""
    cod_from_deltaS(ΔS, ε, ℬ=getbasis(ΔS)) -> Tens{2,3}

Back out ``\\mathbf B`` from an elliptic-crack compliance contribution
``Δ\\mathbb S = π ε (\\hat n ⊗ˢ B ⊗ˢ \\hat n)``.
"""
function cod_from_deltaS(ΔS, ε, ℬ = TensND.getbasis(ΔS))
    T = eltype(ΔS)
    newS = TensND.change_tens(ΔS / (T(π) * ε), ℬ)
    return TensND.Tens(
        Tensors.SymmetricTensor{2, 3}(
            (i, j) ->
            4 * newS[i, 3, j, 3] /
                (
                one(T) + MFH_Core._δ(i, 3, T) + MFH_Core._δ(j, 3, T) +
                    MFH_Core._δ(i, 3, T) * MFH_Core._δ(j, 3, T)
            )
        ),
        ℬ,
    )
end
