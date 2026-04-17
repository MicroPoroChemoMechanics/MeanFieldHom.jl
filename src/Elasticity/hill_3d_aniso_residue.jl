# =============================================================================
#  hill_3d_aniso_residue.jl
#
#  Thin wrapper around the factored polynomial machinery of
#  `Core.green_residue.jl`.  The Hill residue algorithm specific pieces
#  — per-φ ζ(z) parametrisation, 21 Voigt-indexed numerators, Masson
#  log-factor post-processing — live here.  The generic polynomial
#  bookkeeping (acoustic tensor, adjugate, determinant, root finding) is
#  delegated to `Core._build_poly_system`.
#
#  Implementation of Masson (2008) Eq. (20). Float64 only.
# =============================================================================

function _hill_3d_aniso_residue(
        ell::Ellipsoid{3}, C₀;
        abstol::Float64 = 1.0e-8,
        reltol::Float64 = 1.0e-6,
        maxiters::Int = 100_000
    )
    η = Float64(ell.semi_axes[2] / ell.semi_axes[1])
    ω = Float64(ell.semi_axes[3] / ell.semi_axes[1])

    C₀_princ = TensND.change_tens(C₀, ell.basis)
    C = zeros(Float64, 3, 3, 3, 3)
    for i in 1:3, j in 1:3, k in 1:3, l in 1:3
        C[i, j, k, l] = Float64(C₀_princ[i, j, k, l])
    end

    voigt_ij = ((1, 1), (2, 2), (3, 3), (2, 3), (1, 3), (1, 2))

    const_I = ComplexF64(0.0, 1.0)

    function hill_at_phi(φ::Float64)
        m1 = cos(φ)
        m2 = sin(φ) / η
        # ζ(z) = (m1, m2, z/ω) = α₀ + α₁·z with α₀=(m1,m2,0), α₁=(0,0,1/ω).
        α₀ζ = [m1, m2, 0.0]
        α₁ζ = [0.0, 0.0, 1.0 / ω]

        sys = MFH_Core._build_poly_system(C, α₀ζ, α₁ζ)
        K_poly = sys.K_poly       # unused (kept for parity with Cracks residue)
        adj_poly = sys.adj_poly
        Q = sys.Q
        dQ = sys.dQ
        roots_uhp = sys.roots_uhp

        # ζ as a vector of Polynomial{ComplexF64,:z}
        ζ_poly = [Polynomial(ComplexF64[α₀ζ[i], α₁ζ[i]], :z) for i in 1:3]

        # 21 numerator polynomials (upper-Voigt) of degree ≤ 5
        P_num = Vector{Polynomial{ComplexF64, :z}}(undef, 21)
        idx = 0
        for I in 1:6
            vi, vj = voigt_ij[I]
            for J in I:6
                vk, vl = voigt_ij[J]
                idx += 1
                p = (
                    ζ_poly[vi] * (adj_poly[vj, vk] * ζ_poly[vl] + adj_poly[vj, vl] * ζ_poly[vk]) +
                        ζ_poly[vj] * (adj_poly[vi, vk] * ζ_poly[vl] + adj_poly[vi, vl] * ζ_poly[vk])
                )
                P_num[idx] = p * 0.25
            end
        end

        res = zeros(Float64, 21)

        # Term 1: z = i  (compute_residue_log_I at mult=0)
        Qi = Q(const_I)
        if abs2(Qi) > 1.0e-30
            for α in 1:21
                res[α] += real(P_num[α](const_I) / Qi)
            end
        else
            dQi = dQ(const_I)
            d2Qi = derivative(dQ)(const_I)
            for α in 1:21
                p0 = P_num[α](const_I)
                p1 = derivative(P_num[α])(const_I)
                r = (5im * p0 * dQi + 6 * p1 * dQi - 3 * p0 * d2Qi) / (6 * dQi^2)
                res[α] += real(r)
            end
        end

        # Term 2: other UHP roots of Q (compute_residue_log_z at mult=1)
        for zr in roots_uhp
            abs(zr - const_I) > 1.0e-6 || continue
            dQr = dQ(zr)
            abs(dQr) < 1.0e-30 && continue

            t1 = zr * zr
            t2 = 1.0 + t1
            t3 = sqrt(t2)
            L = MFH_Core._masson_log(zr)
            den = dQr * t3 * t2

            for α in 1:21
                res[α] += real(L * P_num[α](zr) / den)
            end
        end

        return res
    end

    P_vals, _ = QuadGK.quadgk(
        hill_at_phi, 0.0, π;
        atol = abstol, rtol = reltol, maxevals = maxiters
    )
    P_vals ./= π

    P_arr = zeros(Float64, 3, 3, 3, 3)
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
