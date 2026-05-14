# =============================================================================
#  moduli.jl
#
#  Extractors for the engineering moduli used by the analytical formulas in
#  the sub-modules `Elasticity`, `Cracks` and `Conductivity`.  The
#  routines only *read* components of the TensND tensor — they do not
#  perform any homogenization step.
# =============================================================================

"""
    extract_iso_moduli(C₀::TensISO{4,3}) -> (E, ν)

Extract Young's modulus `E` and Poisson's ratio `ν` from an isotropic
4th-order stiffness `TensISO{4,3}`.  The internal TensND convention is
``C_0 = 3k\\,\\mathbb J + 2μ\\,\\mathbb K`` i.e. `C₀.data = (3k, 2μ)`.
"""
function extract_iso_moduli(C₀::TensND.TensISO{4, 3})
    α, β = C₀.data           # α = 3k, β = 2μ
    k = α / 3
    μ = β / 2
    E = 9k * μ / (3k + μ)
    ν = (3k - 2μ) / (2 * (3k + μ))
    return E, ν
end

"""
    extract_iso_moduli(C₀::TensISO{4,2}) -> (E, ν)

2D plane-strain counterpart of [`extract_iso_moduli`](@ref) for
`TensISO{4,2}`.  The same 3D formulas are used because the TensND
storage is dimension-agnostic (`TensISO{4,d}` stores the same
`(α, β) = (3k, 2μ)` pair).
"""
function extract_iso_moduli(C₀::TensND.TensISO{4, 2})
    α, β = C₀.data
    k = α / 3
    μ = β / 2
    E = 9k * μ / (3k + μ)
    ν = (3k - 2μ) / (2 * (3k + μ))
    return E, ν
end

"""
    extract_ti_moduli(C₀, n̂) -> (E, H, ν₁, ν₂, Γ)

Read the five TI compliance moduli out of a stiffness tensor `C₀`
whose axis of symmetry is `n̂`.  Used by the closed-form COD formulas
of the `Cracks` sub-module.

The moduli are defined through the compliance tensor
``\\mathbb S = C_0^{-1}`` as:

* ``E = 1/S_{1111}``
* ``H = 1/(S_{3333}·E)``
* ``ν_1 = -E \\cdot S_{1122}``
* ``ν_2 = -E \\cdot S_{1133}``
* ``Γ = (1+ν_1)/(2·E·S_{2323})``

See the package documentation for the full derivation.
"""
function extract_ti_moduli(C₀, n̂)
    𝕊 = inv(C₀)
    ℬ_current = TensND.get_basis(𝕊)
    𝕊_rot = TensND.tens_basis(ℬ_current, 3) == n̂ ? 𝕊 :
        TensND.change_tens(
            𝕊,
            TensND.Basis(TensND.angles(TensND.components_canon(n̂))...)
        )
    E = inv(𝕊_rot[1, 1, 1, 1])
    H = inv(𝕊_rot[3, 3, 3, 3] * E)
    ν₁ = -E * 𝕊_rot[1, 1, 2, 2]
    ν₂ = -E * 𝕊_rot[1, 1, 3, 3]
    Γ = inv(𝕊_rot[2, 3, 2, 3]) * (one(E) + ν₁) / (2 * E)
    return (E = E, H = H, ν₁ = ν₁, ν₂ = ν₂, Γ = Γ)
end

"""
    extract_iso_conductivity(K₀::TensISO{2,d}) -> k

Extract the (scalar) conductivity coefficient of an isotropic 2nd-order
transport tensor `TensISO{2,d}` (``K_0 = k \\cdot \\delta``).
"""
extract_iso_conductivity(K₀::TensND.TensISO{2, 3}) = K₀.data[1]
extract_iso_conductivity(K₀::TensND.TensISO{2, 2}) = K₀.data[1]
