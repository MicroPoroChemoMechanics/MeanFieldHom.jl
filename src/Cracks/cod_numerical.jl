# =============================================================================
#  cod_numerical.jl — numerical COD tensor for anisotropic matrices.
# =============================================================================

"""
    _cod_elliptic_numerical(c, C₀, backend; abstol, reltol, maxiters)
"""
function _cod_elliptic_numerical(c::EllipticCrack{T}, C₀, backend;
                                  abstol::Real = 1e-8,
                                  reltol::Real = 1e-6,
                                  maxiters::Int = 100_000) where {T<:Number}
    η        = aspect_ratio(c)
    lhat, mhat, nhat = MFH_Core._frame_columns(crack_basis(c))
    basis    = crack_basis(c)

    C0_loc = TensND.change_tens(C₀, basis)
    Carr   = MFH_Core._C_array(C0_loc)

    Tp = promote_type(T, eltype(Carr), eltype(lhat))

    function M(φ)
        ξ_plane = similar(lhat, Tp)
        @inbounds for i in 1:3
            ξ_plane[i] = η * cos(φ) * lhat[i] + sin(φ) * mhat[i]
        end
        return backend(Carr, ξ_plane, nhat)
    end

    S = zeros(Tp, 3, 3)
    for i in 1:3, k in i:3
        val, _ = QuadGK.quadgk(φ -> M(φ)[i, k], 0.0, 2π;
                               atol=abstol, rtol=reltol, maxevals=maxiters)
        S[i, k] = val
        if k != i
            S[k, i] = val
        end
    end

    Smat = S / (Tp(8) / Tp(3))
    Bmat = inv(Smat)
    Bsym = (Bmat + transpose(Bmat)) / 2
    return TensND.Tens(Tensors.SymmetricTensor{2,3}((i, j) -> Tp(Bsym[i, j])), basis)
end

"""
    _cod_ribbon_numerical(c, C₀, backend; …)
"""
function _cod_ribbon_numerical(c::RibbonCrack{T}, C₀, backend;
                                abstol::Real = 1e-8,
                                reltol::Real = 1e-6,
                                maxiters::Int = 100_000) where {T<:Number}
    _, mhat, nhat = MFH_Core._frame_columns(crack_basis(c))
    basis   = crack_basis(c)
    C0_loc  = TensND.change_tens(C₀, basis)
    Carr    = MFH_Core._C_array(C0_loc)

    Tp = promote_type(T, eltype(Carr), eltype(mhat))

    mhat_p = Tp.(mhat)
    nhat_p = Tp.(nhat)

    Qstar = backend(Carr, mhat_p, nhat_p;
                    abstol=abstol, reltol=reltol, maxiters=maxiters)

    Bmat = inv(Qstar) * (Tp(π) / 4)
    Bsym = (Bmat + transpose(Bmat)) / 2
    return TensND.Tens(Tensors.SymmetricTensor{2,3}((i, j) -> Tp(Bsym[i, j])), basis)
end

# Backend dispatch
_residue_backend(Carr, ξ_plane, nhat; kw...) = _Qnn_star_residue(Carr, ξ_plane, nhat)

_decuhr_backend(Carr, ξ_plane, nhat; abstol, reltol, maxiters) =
    _Qnn_star_decuhr(Carr, ξ_plane, nhat; abstol=abstol, reltol=reltol, maxiters=maxiters)
