# =============================================================================
#  hill_3d_aniso_decuhr.jl
#
#  ForwardDiff-compatible 2D cubature for the 3D anisotropic Hill tensor.
#  Nested QuadGK over (z, φ) ∈ [0, 1] × [0, 2π].
#
#  The routine uses only `Core._acoustic_tensor` + `Core._inv3` from the
#  Core layer; the specific 2D parametrisation (`ζ = (√(1-z²)cosφ,
#  √(1-z²)sinφ/η, z/ω)`) is kept here because it is tied to the ellipsoid
#  semi-axes (η, ω).
# =============================================================================

function _hill_3d_aniso_decuhr(
        ell::Ellipsoid{3}, C₀;
        abstol::Float64 = 1.0e-8,
        reltol::Float64 = 1.0e-6,
        maxiters::Int = 1_000_000
    )
    T = promote_type(eltype(ell.semi_axes), eltype(C₀))

    η = ell.semi_axes[2] / ell.semi_axes[1]
    ω = ell.semi_axes[3] / ell.semi_axes[1]

    C₀_princ = TensND.change_tens(C₀, ell.basis)
    C₀_arr = Array{T, 4}(undef, 3, 3, 3, 3)
    for i in 1:3, j in 1:3, k in 1:3, l in 1:3
        C₀_arr[i, j, k, l] = C₀_princ[i, j, k, l]
    end

    voigt_ij = ((1, 1), (2, 2), (3, 3), (2, 3), (1, 3), (1, 2))

    function integrand_at(z, φ)
        uz = sqrt(one(T) - T(z) * T(z))
        ζ = (uz * cos(T(φ)), uz * sin(T(φ)) / η, T(z) / ω)

        K = Matrix{T}(undef, 3, 3)
        for i in 1:3, j in 1:3
            s = zero(T)
            for k in 1:3, l in 1:3
                s += ζ[k] * C₀_arr[k, i, j, l] * ζ[l]
            end
            K[i, j] = s
        end
        iK = inv(K)

        vals = Vector{T}(undef, 21)
        idx = 0
        for I in 1:6
            i, j = voigt_ij[I]
            for J in I:6
                k, l = voigt_ij[J]
                γ = T(0.25) * (
                    ζ[i] * (iK[j, k] * ζ[l] + iK[j, l] * ζ[k]) +
                        ζ[j] * (iK[i, k] * ζ[l] + iK[i, l] * ζ[k])
                )
                idx += 1
                vals[idx] = γ
            end
        end
        return vals
    end

    P_vals, _ = QuadGK.quadgk(
        0.0, 1.0;
        atol = abstol,
        rtol = reltol,
        maxevals = maxiters
    ) do z
        inner, _ = QuadGK.quadgk(
            0.0, 2π;
            atol = abstol / 10,
            rtol = reltol
        ) do φ
            integrand_at(z, φ)
        end
        inner
    end
    P_vals ./= T(2π)

    P_arr = zeros(T, 3, 3, 3, 3)
    idx = 0
    for I in 1:6, J in I:6
        i, j = voigt_ij[I]
        k, l = voigt_ij[J]
        idx += 1
        v = P_vals[idx]
        P_arr[i, j, k, l] = v;  P_arr[j, i, k, l] = v
        P_arr[i, j, l, k] = v;  P_arr[j, i, l, k] = v
        P_arr[k, l, i, j] = v;  P_arr[l, k, i, j] = v
        P_arr[k, l, j, i] = v;  P_arr[l, k, j, i] = v
    end

    return TensND.change_tens_canon(TensND.Tens(P_arr, ell.basis))
end
