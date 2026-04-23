# =============================================================================
#  cod_analytical_thermal.jl — closed-form thermal COD scalars.
#
#  Analogue of `cod_analytical.jl` for the 2nd-order (conductivity /
#  diffusion) problem.  Returns a **scalar** `b` instead of a symmetric
#  2nd-order tensor: since the temperature jump across a flat crack is a
#  scalar driven only by the normal component of the heat flux, a single
#  scalar suffices.  The associated resistivity contribution is then
#  assembled as
#
#      ΔR = π ε · b · n̂ ⊗ n̂   (elliptic crack)
#      ΔR = (π/2) ε · b · n̂ ⊗ n̂   (ribbon crack)
#
#  in [`compliance.jl`](src/Cracks/compliance.jl).  Mathematical
#  derivation and reference conventions live in
#  [`docs/src/theory/thermal_cracks.md`](docs/src/theory/thermal_cracks.md).
# =============================================================================

# -----------------------------------------------------------------------------
#  ISO matrix
# -----------------------------------------------------------------------------

"""
    _cod_iso_ellipse_thermal(c::EllipticCrack, k₀) -> Real

Closed-form thermal COD scalar ``b`` of an elliptic crack of aspect
ratio ``\\eta = b/a`` in an isotropic conductor ``\\mathbf K_0 =
k_0\\,\\mathbf 1``:

```
b = η / (π · k₀ · 𝓔_η) ,   𝓔_η = 𝓔(√(1-η²))
```

with ``\\mathcal E_\\eta`` the complete elliptic integral of the
second kind ([Abramowitz & Stegun 1972](@cite abramowitz1972)).
Penny-crack limit ``\\eta=1``: ``b = 2/(π^{2} k_{0})``.
"""
function _cod_iso_ellipse_thermal(c::EllipticCrack{T}, k₀) where {T <: Number}
    η = aspect_ratio(c)
    _, _, ℰ = _elliptic_CS(η)
    return T(η) / (T(π) * k₀ * ℰ)
end

"""
    _cod_iso_ribbon_thermal(c::RibbonCrack, k₀) -> Real

Closed-form thermal COD scalar of a ribbon (tunnel) crack in an
isotropic conductor:

```
b = 2 / (π k₀) .
```

Ribbon limit of the 2-D thermal Eshelby problem
([Sevostianov & Kachanov 2002](@cite sevostianov2002),
 [Kachanov 2018](@cite kachanov2018)).
"""
function _cod_iso_ribbon_thermal(c::RibbonCrack{T}, k₀) where {T <: Number}
    return T(2) / (T(π) * k₀)
end

# -----------------------------------------------------------------------------
#  Anisotropic matrix — elliptic crack (K⁻¹ᐟ² transform, Giraud-Gruescu 2019)
# -----------------------------------------------------------------------------

"""
    _cod_aniso_ellipse_thermal(c::EllipticCrack, K₀) -> Real

Closed-form thermal COD scalar of an elliptic crack in an arbitrarily
anisotropic conductor, via the square-root change-of-variable of
[Giraud & Gruescu 2019](@cite giraudMOM2019). Let
``\\mathbf A = \\mathbf R_c \\cdot \\mathrm{diag}(a,b,0) \\cdot
\\mathbf R_c^{T}`` be the (flat) crack shape tensor and
``\\tilde{\\mathbf A} = \\mathbf A \\cdot \\mathbf K_0^{-1/2}``.  Its
singular values ``\\sigma_1 \\ge \\sigma_2 \\ge \\sigma_3 = 0`` give the
transformed in-plane aspect ratio
``\\eta_t = \\sigma_2/\\sigma_1 \\in (0,1]``.  Then

```
b = σ₂ · (n̂·K₀⁻¹·n̂) · √(n̂·K₀·n̂) / (π · aₘₐₓ · 𝓔(√(1-η_t²))) ,
```

with ``a_\\text{max} = \\max(a,b)`` and the resistivity contribution
rank-1 along ``\\hat{\\mathbf w} = \\mathbf K_0^{-1/2}\\hat{\\mathbf n}/
\\sqrt{\\hat{\\mathbf n}\\cdot\\mathbf K_0^{-1}\\hat{\\mathbf n}}``
(assembled by [`compliance_contribution`](@ref)).  Reduces to the iso
formula ``b = \\eta/(\\pi k_0 \\mathcal E_\\eta)`` for
``\\mathbf K_0 = k_0\\,\\mathbf 1`` (where ``\\hat{\\mathbf w} = \\hat
{\\mathbf n}``), and to ``b = 2/(\\pi^2 \\sqrt{k_t k_n})`` for a penny
crack in a TI matrix aligned with ``\\hat{\\mathbf n}``.
"""
function _cod_aniso_ellipse_thermal(c::EllipticCrack, K₀)
    T_mat = eltype(K₀)

    K_arr = Matrix{T_mat}(undef, 3, 3)
    for i in 1:3, j in 1:3
        K_arr[i, j] = K₀[i, j]
    end

    F = eigen(Symmetric(K_arr))
    invsqrt_K = F.vectors * Diagonal(1 ./ sqrt.(F.values)) * F.vectors'

    R_crack = [crack_basis(c)[i, j] for i in 1:3, j in 1:3]
    A = R_crack * Diagonal([c.a, c.b, zero(c.a)]) * R_crack'

    σ = sort(svdvals(A * invsqrt_K); rev = true)

    η_t = σ[2] / σ[1]
    _, _, ℰ = _elliptic_CS(η_t)

    n̂ = R_crack[:, 3]
    k_nn = dot(n̂, K_arr * n̂)
    k_nn_inv = dot(n̂, K_arr \ n̂)
    a_max = max(c.a, c.b)

    return σ[2] * k_nn_inv * sqrt(k_nn) / (π * a_max * ℰ)
end

# -----------------------------------------------------------------------------
#  Anisotropic matrix — ribbon crack (2-D K⁻¹ᐟ² on transverse plane)
# -----------------------------------------------------------------------------

"""
    _cod_aniso_ribbon_thermal(c::RibbonCrack, K₀) -> Real

Closed-form thermal COD scalar of a ribbon crack in an arbitrarily
anisotropic conductor.  Only the 2×2 block of ``\\mathbf K_0``
restricted to the plane ``(\\hat{\\mathbf m}, \\hat{\\mathbf n})``
(spanned by the in-plane crack direction and the crack normal) enters
the formula:

```
b = 2 / (π · √det(K₀|_{(m̂,n̂)})) .
```

For an isotropic conductor this reduces to ``b = 2/(\\pi k_0)``.
Derivation via the 2-D Giraud-Gruescu square-root transform on the
transverse plane of the cylinder
([Giraud & Gruescu 2019](@cite giraudMOM2019)).
"""
function _cod_aniso_ribbon_thermal(c::RibbonCrack, K₀)
    T_mat = eltype(K₀)

    # Express K₀ in the crack basis: indices 2, 3 give the (m̂, n̂) block.
    K₀_loc = TensND.change_tens(K₀, crack_basis(c))
    K_mm = T_mat(K₀_loc[2, 2])
    K_nn = T_mat(K₀_loc[3, 3])
    K_mn = T_mat(K₀_loc[2, 3])

    det_K⊥ = K_mm * K_nn - K_mn * K_mn
    return T_mat(2) / (T_mat(π) * sqrt(det_K⊥))
end
