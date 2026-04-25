# =============================================================================
#  hill_3d_ti_coaxial.jl — Analytical Hill polarisation tensor for a
#  spheroidal inclusion coaxial with a transversely isotropic matrix.
#
#  Implements the closed-form formula derived in
#  [Barthélémy 2020](@cite barthelemy2020) (eqs. 49–58 of the post-print)
#  for the Walpole-basis components `(P₁, P₂, P₃, P₅, P₆)` of the Hill
#  tensor of a spheroid (prolate or oblate) whose axis is parallel to
#  the symmetry axis of a TI elastic matrix.
#
#  The five elastic constants (C₁₁₁₁, C₁₁₂₂, C₁₁₃₃, C₃₃₃₃, C₂₃₂₃) of
#  the TI matrix and the spheroid aspect ratio `ω = (axial)/(transverse)`
#  enter through six elementary elliptic integrals
#  (I₀, I₂, I₄, J₀, J₂, J₄) that admit closed forms in terms of `acosh`
#  and complex square roots.
#
#  The output is a `TensTI{4, T, 5}` (major-symmetric Walpole tensor)
#  with axis equal to the common symmetry axis of `C₀` and the
#  spheroid.  When the matrix is in fact isotropic, the formula
#  recovers the classical Eshelby–Mura result
#  ([Mura 1987](@cite mura1987)).
# =============================================================================

# ── Elliptic-integral helpers — single argument η ----------------------------
#
#   I₀(η) = ∫₋₁¹ dz / [z² + η² (1−z²)]
#   I₂(η) = ∫₋₁¹ z² dz / [z² + η² (1−z²)]
#   I₄(η) = ∫₋₁¹ z⁴ dz / [z² + η² (1−z²)]
#
# Closed forms when η ≠ 1, with limit values 2, 2/3, 2/5 at η → 1.

# Tolerances chosen to absorb floating-point noise from acoustic-polynomial
# discriminants that should vanish exactly (e.g. when the TI matrix collapses
# to an isotropic one).  `_HILL_TI_EPS` is used both for `η ≈ 1` (single
# argument) and `η₁ ≈ η₂` (two arguments).
const _HILL_TI_EPS = 1.0e-6

@inline function _I0(η)
    k2 = η^2 - 1
    abs(k2) < _HILL_TI_EPS && return 2 * one(η)
    return 2 * acosh(η) / (η * sqrt(k2))
end

@inline function _I2(η)
    k2 = η^2 - 1
    abs(k2) < _HILL_TI_EPS && return 2 * one(η) / 3
    return 2 * (η * acosh(η) - sqrt(k2)) / k2^(3 // 2)
end

@inline function _I4(η)
    k2 = η^2 - 1
    abs(k2) < _HILL_TI_EPS && return 2 * one(η) / 5
    return 2 * (η^3 * acosh(η) - (1 + 4 * k2 / 3) * sqrt(k2)) / k2^(5 // 2)
end

# ── Elliptic-integral helpers — two arguments η₁, η₂ -------------------------
#
#   J_n(η₁,η₂) = ∫₋₁¹ z^n dz / { [z² + η₁²(1−z²)] [z² + η₂²(1−z²)] }
#
# Coincidence rule J_n(η,η) handled separately to avoid the 0/0 form.
# `η₁ ≈ η₂` is detected with a relative tolerance to absorb round-off
# noise from a near-degenerate acoustic polynomial (e.g. iso matrix
# expressed as a TI tensor).

@inline _close_eta(η1, η2) = abs(η1^2 - η2^2) < _HILL_TI_EPS * max(abs(η1^2), abs(η2^2), 1.0)

@inline function _J0(η1, η2)
    if _close_eta(η1, η2)
        η = (η1 + η2) / 2
        k2 = η^2 - 1; sk2 = sqrt(k2)
        abs(k2) < _HILL_TI_EPS && return 2 * one(η)
        return (acosh(η) + η * sk2) / (sk2 * η^3)
    end
    e12 = η1^2; e22 = η2^2
    return ((1 - e12) * _I0(η1) - (1 - e22) * _I0(η2)) / (e22 - e12)
end

@inline function _J2(η1, η2)
    if _close_eta(η1, η2)
        η = (η1 + η2) / 2
        k2 = η^2 - 1; sk2 = sqrt(k2)
        abs(k2) < _HILL_TI_EPS && return 2 * one(η) / 3
        return (-acosh(η) + η * sk2) / (sk2^3 * η)
    end
    e12 = η1^2; e22 = η2^2
    return ((1 - e12) * _I2(η1) - (1 - e22) * _I2(η2)) / (e22 - e12)
end

@inline function _J4(η1, η2)
    if _close_eta(η1, η2)
        η = (η1 + η2) / 2
        k2 = η^2 - 1; sk2 = sqrt(k2)
        abs(k2) < _HILL_TI_EPS && return 2 * one(η) / 5
        return (-3 * η * acosh(η) + (2 + η^2) * sk2) / sk2^5
    end
    e12 = η1^2; e22 = η2^2
    return ((1 - e12) * _I4(η1) - (1 - e22) * _I4(η2)) / (e22 - e12)
end

# ── Core routine --------------------------------------------------------------

"""
    _hill_ti_walpole(ω, C1111, C1122, C1133, C3333, C2323)
        -> (P1, P2, P3, P5, P6)

Closed-form Walpole-basis coefficients of the Hill tensor for a spheroid
of aspect ratio `ω = (axial)/(transverse)` coaxial with a transversely
isotropic matrix specified by the five independent elastic constants.

The arithmetic is carried out in `Complex{T}` where
`T = promote_type(typeof(ω), typeof(C1111), …, typeof(C2323))`. This
preserves compatibility with `ForwardDiff.Dual` numbers (so the
analytical Hill tensor is differentiable through the elastic constants
and the spheroid aspect ratio).  The imaginary parts cancel by
construction and only the real parts are returned. When the matrix is
in fact isotropic (`C1111 = C3333 = λ + 2μ`, `C1122 = C1133 = λ`,
`C2323 = μ`) the returned coefficients reduce to the classical Mura
formula.

!!! note "Symbolic numbers (SymPy `Sym`)"
    The function is **not** compatible with SymPy `Sym` inputs because
    Julia's `Complex{T}` parametric type requires `T <: Real`, which
    `Sym` does not satisfy. For symbolic exploration of the analytical
    formula, derive the components by hand from eqs. 49–58 of the
    paper, or use the residue/DECUHR backends with substituted
    numerical values.

Reference: [barthelemyIJES2020_hilltrans](@cite), eqs. 49–58.
"""
function _hill_ti_walpole(ω, C1111, C1122, C1133, C3333, C2323)
    T = promote_type(typeof(ω), typeof(C1111), typeof(C1122),
                     typeof(C1133), typeof(C3333), typeof(C2323))
    ωc  = Complex{T}(ω)
    C11 = Complex{T}(C1111)
    C12 = Complex{T}(C1122)
    C13 = Complex{T}(C1133)
    C33 = Complex{T}(C3333)
    C44 = Complex{T}(C2323)

    om2 = ωc^2
    om4 = om2^2

    # Roots γ₁, γ₂ of  a γ⁴ + b γ² + c = 0  (acoustic polynomial)
    a = C44 * C33
    b = C13^2 + 2 * C13 * C44 - C11 * C33
    c = C11 * C44
    sqd = sqrt(b^2 - 4 * a * c)
    γ1 = sqrt((-b + sqd) / (2 * a))
    γ2 = sqrt((-b - sqd) / (2 * a))

    η1 = ωc * γ1
    η2 = ωc * γ2
    η3 = ωc * sqrt((C11 - C12) / (2 * C44))

    coef = C44 * C33
    j0 = _J0(η1, η2) / coef
    j2 = _J2(η1, η2) / coef
    j4 = _J4(η1, η2) / coef
    i0 = _I0(η3) / C44
    i2 = _I2(η3) / C44

    # Walpole components (eqs. 53–58 of Barthélémy 2020)
    P1 = real(((C44 - om2 * C11) * j4 + om2 * C11 * j2) / 2)
    P2 = real(om2 * ((om2 * C44 - C33) * j4 + (C33 - 2 * om2 * C44) * j2 + om2 * C44 * j0) / 4)
    P3 = real(om2 / (2 * sqrt(T(2))) * (C44 + C13) * (j4 - j2))
    P5 = real(P2) / 2 + real(om2 * (i0 - i2) / 8)
    P6 = real(((om4 * C11 + C33 + 2 * om2 * C13) * j4
                - 2 * om2 * (om2 * C11 + C13) * j2
                + om4 * C11 * j0 + i2) / 8)

    return (P1, P2, P3, P5, P6)
end

# ── Public builders — coaxial spheroid in TI matrix --------------------------

"""
    _hill_3d_ti_coaxial(ell, C₀) -> TensTI{4, Float64, 5}

Hill polarisation tensor for a spheroidal inclusion `ell` coaxial with a
transversely isotropic matrix `C₀::TensTI{4, T, 5}`.  Both must share
the same symmetry axis; coaxiality is checked by the dispatcher
([`_ti_coaxial`](@ref)) and the `Analytical` algorithm is selected when
the test succeeds.  Otherwise the user is silently routed to the
generic residue/DECUHR backend.

# Arguments
- `ell::Ellipsoid{3, Oblate}`: oblate spheroid `a = b ≥ c`, axis `e₃`,
  aspect ratio `ω = c/a`.
- `ell::Ellipsoid{3, Prolate}`: prolate spheroid `a ≥ b = c`, axis `e₁`,
  aspect ratio `ω = a/b`.

# Returns
`TensTI{4, Float64, 5}` (major-symmetric Walpole tensor) with axis equal
to the matrix's TI axis.

Reference: [Barthélémy 2020](@cite barthelemy2020).
"""
function _hill_3d_ti_coaxial(ell::Ellipsoid{3, Oblate}, C₀::TensND.TensTI{4, T, 5}) where {T}
    a, _, c = ell.semi_axes
    ω = c / a   # generic division — preserves Dual derivatives
    C1111, C1122, C1133, C3333, C2323 = TensND.arg_TI(C₀)
    P1, P2, P3, P5, P6 = _hill_ti_walpole(ω, C1111, C1122, C1133, C3333, C2323)
    return TensND.TensTI{4}(P1, P2, P3, P5, P6, TensND.axis(C₀))
end

function _hill_3d_ti_coaxial(ell::Ellipsoid{3, Prolate}, C₀::TensND.TensTI{4, T, 5}) where {T}
    a, b, _ = ell.semi_axes
    ω = a / b   # generic division — preserves Dual derivatives
    C1111, C1122, C1133, C3333, C2323 = TensND.arg_TI(C₀)
    P1, P2, P3, P5, P6 = _hill_ti_walpole(ω, C1111, C1122, C1133, C3333, C2323)
    return TensND.TensTI{4}(P1, P2, P3, P5, P6, TensND.axis(C₀))
end

function _hill_3d_ti_coaxial(::Ellipsoid{3, Spherical}, C₀::TensND.TensTI{4, T, 5}) where {T}
    # Spheroid degenerates to sphere — ω = 1, formula limit values
    C1111, C1122, C1133, C3333, C2323 = TensND.arg_TI(C₀)
    P1, P2, P3, P5, P6 = _hill_ti_walpole(one(T), C1111, C1122, C1133, C3333, C2323)
    return TensND.TensTI{4}(P1, P2, P3, P5, P6, TensND.axis(C₀))
end

# ── Coaxiality test ----------------------------------------------------------

"""
    _ti_coaxial(C₀::TensTI{4}, ell::Ellipsoid{3}) -> Bool

`true` when the TI symmetry axis of `C₀` is parallel to the (unique)
spheroid axis of `ell` (column 3 for `Oblate`/`Spherical`, column 1 for
`Prolate`).  Used by the dispatcher to route coaxial-spheroid problems
to the analytical TI builder.
"""
function _ti_coaxial(C₀::TensND.TensTI{4}, ell::Ellipsoid{3, Oblate})
    axis_C = collect(TensND.axis(C₀))
    axis_e = collect(TensND.components_canon(TensND.tens_basis(ell.basis, 3)))
    return isapprox(abs(dot(axis_C, axis_e)), 1.0; atol = 1.0e-10)
end

function _ti_coaxial(C₀::TensND.TensTI{4}, ell::Ellipsoid{3, Prolate})
    axis_C = collect(TensND.axis(C₀))
    axis_e = collect(TensND.components_canon(TensND.tens_basis(ell.basis, 1)))
    return isapprox(abs(dot(axis_C, axis_e)), 1.0; atol = 1.0e-10)
end

# Sphere is coaxial with anything (no preferred axis on the inclusion side).
_ti_coaxial(::TensND.TensTI{4}, ::Ellipsoid{3, Spherical}) = true

# Triaxial ellipsoids cannot be coaxial with a TI matrix in general.
_ti_coaxial(::TensND.TensTI{4}, ::Ellipsoid{3, Triaxial}) = false
