# =============================================================================
#  green_residue.jl
#
#  Factored Masson / Cauchy residue algorithm used by both `Elasticity`
#  (Hill polarisation tensor — Masson 2008) and `Cracks` (line-integrated
#  Green kernel ``\\hat{\\mathbf Q}^{\\star}_{nn}`` — Cauchy residue).
#
#  Strategy
#  ========
#  Both algorithms reduce their integrand to the same rational-function
#  pattern :
#
#      f(z) = N(z) / Q(z)
#
#  where
#
#    * `Q(z) = det K(z)` — degree-6 polynomial in `z`
#    * `K(z)` — acoustic tensor evaluated on a caller-supplied parametric
#               wave-vector `ζ(z) = α₀ + α₁ z` (linear polynomial)
#
#  The common building blocks are factored in [`_build_poly_system`](@ref):
#  the coefficient matrices `α₀, α₁, α₂` such that
#  `K(z) = α₀ + α₁ z + α₂ z²`, the polynomial adjugate `adj(K)(z)` and the
#  determinant `Q(z)`, plus its derivative `Q'(z)` and the roots in the
#  upper half plane.
#
#  The algorithm-specific numerator and post-processing (Masson log
#  factor vs Cauchy residue without log) are delegated to a caller-supplied
#  closure in the downstream wrappers — see
#  `Elasticity/hill_3d_aniso_residue.jl` and
#  `Cracks/green_residue.jl`.
# =============================================================================

"""
    _build_acoustic_coeffs(C, α₀ζ, α₁ζ) -> (A₀, A₁, A₂)

Build the three 3×3 coefficient matrices of `K(z) = A₀ + A₁·z + A₂·z²`
given the stiffness `C` (Float64 `3×3×3×3` array) and the two linear
coefficients of the wave-vector parametrisation `ζ(z) = α₀ζ + α₁ζ · z`.

Each `Aₖ` is a 3×3 matrix of real scalars ; they are assembled below
into `ComplexF64` polynomial entries before root finding.
"""
function _build_acoustic_coeffs(
        C::AbstractArray{<:Real, 4},
        α₀ζ::AbstractVector{<:Real},
        α₁ζ::AbstractVector{<:Real}
    )
    A₀ = zeros(Float64, 3, 3)
    A₁ = zeros(Float64, 3, 3)
    A₂ = zeros(Float64, 3, 3)
    @inbounds for j in 1:3, i in 1:3
        s0 = 0.0; s1 = 0.0; s2 = 0.0
        for k in 1:3, l in 1:3
            cc = Float64(C[i, k, j, l])
            s0 += cc * α₀ζ[k] * α₀ζ[l]
            s1 += cc * (α₀ζ[k] * α₁ζ[l] + α₁ζ[k] * α₀ζ[l])
            s2 += cc * α₁ζ[k] * α₁ζ[l]
        end
        A₀[i, j] = s0; A₁[i, j] = s1; A₂[i, j] = s2
    end
    return A₀, A₁, A₂
end

"""
    _build_poly_system(C, α₀ζ, α₁ζ) -> (K_poly, adj_poly, Q, dQ, roots_uhp)

Build the polynomial acoustic tensor `K(z) = A₀ + A₁ z + A₂ z²`, its
polynomial adjugate `adj(K)(z)`, the determinant `Q(z) = det K(z)` and
its derivative `dQ(z)`, then compute the roots of `Q` lying in the
strict upper half plane.

All polynomial entries are `Polynomial{ComplexF64,:z}` (so that complex
roots can be handled uniformly). The roots filter keeps only points
with `imag(zr) > 1e-8`.
"""
function _build_poly_system(
        C::AbstractArray{<:Real, 4},
        α₀ζ::AbstractVector{<:Real},
        α₁ζ::AbstractVector{<:Real}
    )
    A₀, A₁, A₂ = _build_acoustic_coeffs(C, α₀ζ, α₁ζ)

    poly(coefs) = Polynomial(ComplexF64.(coefs), :z)
    K_poly = Matrix{Polynomial{ComplexF64, :z}}(undef, 3, 3)
    @inbounds for j in 1:3, i in 1:3
        K_poly[i, j] = poly([A₀[i, j], A₁[i, j], A₂[i, j]])
    end

    adj_poly = Matrix{Polynomial{ComplexF64, :z}}(undef, 3, 3)
    @inbounds for j in 1:3, i in 1:3
        rows = setdiff(1:3, [j])
        cols = setdiff(1:3, [i])
        minor = K_poly[rows[1], cols[1]] * K_poly[rows[2], cols[2]] -
            K_poly[rows[1], cols[2]] * K_poly[rows[2], cols[1]]
        adj_poly[i, j] = isodd(i + j) ? -minor : minor
    end

    Q = Polynomial(ComplexF64[0.0], :z)
    @inbounds for i in 1:3
        Q = Q + K_poly[1, i] * adj_poly[i, 1]
    end
    dQ = derivative(Q)

    Q_coeffs = coeffs(Q)
    while length(Q_coeffs) > 1 && abs(Q_coeffs[end]) < 1.0e-14 * abs(Q_coeffs[1])
        pop!(Q_coeffs)
    end
    roots_all = PolynomialRoots.roots(Q_coeffs)

    # Newton polish: `PolynomialRoots.roots` (Durand-Kerner) typically leaves
    # |Q(zᵣ)| ~ 1e-10…1e-12; 2 Newton iterations push this to O(eps).  Skip
    # roots where |Q'(zᵣ)| is too small (multiple roots).
    Q_scale = maximum(abs, Q_coeffs)
    dQ_thresh = 1.0e-14 * Q_scale
    @inbounds for _ in 1:2
        for i in eachindex(roots_all)
            zr = roots_all[i]
            Qv = Q(zr)
            dQv = dQ(zr)
            abs(dQv) < dQ_thresh && continue
            roots_all[i] = zr - Qv / dQv
        end
    end

    # After polish, genuine UHP roots have imag > O(eps); tighten the filter
    # to catch only real-axis roots.
    roots_uhp = ComplexF64[zr for zr in roots_all if imag(zr) > 1.0e-12]

    return (K_poly = K_poly, adj_poly = adj_poly, Q = Q, dQ = dQ, roots_uhp = roots_uhp)
end

# ─── Masson log factor (Hill-style) ──────────────────────────────────────────

"""
    _masson_log(z) -> ComplexF64

Masson (2008) log-factor ``L(z) = -(2\\log(z + \\sqrt{1+z²}) - iπ)``.
Used by the Hill residue algorithm to fold the analytic continuation of
the 1-D integral of the Green kernel back onto the real axis.
"""
@inline function _masson_log(z::ComplexF64)
    t2 = 1.0 + z * z
    t3 = sqrt(t2)
    t9 = log(z + t3)
    return -(2.0 * t9 - im * π)
end
