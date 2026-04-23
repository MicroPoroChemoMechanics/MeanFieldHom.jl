# =============================================================================
#  compliance.jl — crack compliance / resistivity contribution tensors.
#
#  Public API returns the size-independent contribution tensor H (or R in
#  the thermal case), consistent with the Echoes convention.  Helpers
#  `delta_compliance` / `delta_resistivity` apply the Budiansky density
#  factor to recover the full dilute compliance correction ΔS (or ΔR).
# =============================================================================

"""
    compliance_contribution(crack, C₀; method=:auto, kw...) -> Tens{4,3}

Size-independent **crack compliance contribution tensor**
``\\mathbb H`` (Echoes convention).  Assembled from the COD tensor
``\\mathbf B = `` [`cod_tensor`](@ref) through the factorisation

- Elliptic crack:  ``\\mathbb H = \\tfrac{3}{4}\\,\\hat{\\mathbf n}
  \\stackrel{s}{\\otimes}\\mathbf B\\stackrel{s}{\\otimes}\\hat{\\mathbf n}``.
- Ribbon crack:   ``\\mathbb H = \\tfrac{2}{\\pi}\\,\\hat{\\mathbf n}
  \\stackrel{s}{\\otimes}\\mathbf B\\stackrel{s}{\\otimes}\\hat{\\mathbf n}``.

The two geometric prefactors follow from the single definition
``\\mathbb H = \\lim_{c/b\\to 0}(c/b)\\,\\mathbb Q^{-1}
= (cS/V)\\,\\hat{\\mathbf n}\\stackrel{s}{\\otimes}\\mathbf B
\\stackrel{s}{\\otimes}\\hat{\\mathbf n}`` applied with
``S/V = 3/(4c)`` (elliptic, ``S=\\pi ab``, ``V=\\tfrac{4}{3}\\pi abc``)
or ``S/V = 2/(\\pi c)`` (ribbon, ``S=4ab``, ``V=2\\pi abc``,
``a\\to\\infty``).  The Kachanov–Echoes factorisation of the elliptic
case is recovered for ``\\eta=1`` (penny).

Apply [`delta_compliance`](@ref)`(crack, H, ε)` to obtain the dilute
compliance correction ``\\Delta\\mathbb S``:

```
ΔS = (4π/3) ε³ᵈ H   (elliptic, ε³ᵈ = N a b²)
ΔS =    π   ε²ᵈ H   (ribbon,   ε²ᵈ = N b²)
```

See [Kachanov 1992](@cite kachanov1992),
[Sevostianov & Kachanov 2002](@cite sevostianov2002),
[Barthélémy 2021](@cite barthelemyIJES2021).
"""
function compliance_contribution(crack::MFH_Core.AbstractCrack, C₀::TensND.AbstractTens{4, 3}; kw...)
    B = cod_tensor(crack, C₀; kw...)
    return _compliance_from_B(crack, B)
end

"""
    _compliance_from_B_elliptic(crack, B) -> Tens{4,3}

Elliptic: ``\\mathbb H = \\tfrac{3}{4}(\\hat n ⊗ˢ \\mathbf B ⊗ˢ \\hat n)``.
"""
function _compliance_from_B_elliptic(crack::EllipticCrack, B)
    n̂ = TensND.tensbasis(crack_basis(crack), 3)
    T = eltype(B)
    return (3 * one(T) / 4) * (n̂ ⊗ˢ B ⊗ˢ n̂)
end

"""
    _compliance_from_B_ribbon(crack, B) -> Tens{4,3}

Ribbon: ``\\mathbb H = \\tfrac{2}{\\pi}(\\hat n ⊗ˢ \\mathbf B ⊗ˢ \\hat n)``.
"""
function _compliance_from_B_ribbon(crack::RibbonCrack, B)
    n̂ = TensND.tensbasis(crack_basis(crack), 3)
    T = eltype(B)
    return (2 * one(T) / T(π)) * (n̂ ⊗ˢ B ⊗ˢ n̂)
end

_compliance_from_B(crack::EllipticCrack, B) = _compliance_from_B_elliptic(crack, B)
_compliance_from_B(crack::RibbonCrack, B) = _compliance_from_B_ribbon(crack, B)

# =============================================================================
#  Thermal (2nd-order) — R from the scalar COD b
# =============================================================================

"""
    _effective_crack_direction(crack, K₀) -> Tens{1,3}

Unit vector ``\\hat{\\mathbf w}`` giving the rank-1 direction of the
crack resistivity contribution.  For an isotropic matrix or a matrix
aligned with the crack frame this reduces to the geometric normal
``\\hat{\\mathbf n}``; for general anisotropy,

```
ŵ = K₀⁻¹ᐟ² · n̂ / ‖K₀⁻¹ᐟ² · n̂‖
```

(obtained as the null-space direction of ``\\mathbf K_0 - \\mathbf K_0
\\mathbf P(0) \\mathbf K_0`` in the crack limit — see
[`docs/src/theory/thermal_cracks.md`](docs/src/theory/thermal_cracks.md)).
"""
function _effective_crack_direction(crack::MFH_Core.AbstractCrack, K₀)
    T_mat = eltype(K₀)
    K_arr = Matrix{T_mat}(undef, 3, 3)
    for i in 1:3, j in 1:3
        K_arr[i, j] = K₀[i, j]
    end
    F = eigen(Symmetric(K_arr))
    invsqrt_K = F.vectors * Diagonal(1 ./ sqrt.(F.values)) * F.vectors'
    n̂_arr = [crack_basis(crack)[i, 3] for i in 1:3]
    w = invsqrt_K * n̂_arr
    w ./= norm(w)
    return TensND.Tens(w)
end

"""
    _resistivity_from_b_elliptic(crack, b, K₀) -> Tens{2,3}

Elliptic: ``\\mathbf R = \\tfrac{3}{4}\\,b\\,\\hat{\\mathbf w}
\\otimes\\hat{\\mathbf w}`` where ``\\hat{\\mathbf w}`` is the effective
rank-1 direction (see [`_effective_crack_direction`](@ref)).  Reduces
to ``\\mathbf R = \\tfrac{3}{4}\\,b\\,\\hat{\\mathbf n}\\otimes\\hat{\\mathbf n}``
for iso / aligned-TI matrices.
"""
function _resistivity_from_b_elliptic(crack::EllipticCrack, b, K₀)
    ŵ = _effective_crack_direction(crack, K₀)
    T = eltype(ŵ)
    return (T(3) / T(4) * b) * (ŵ ⊗ ŵ)
end

"""
    _resistivity_from_b_ribbon(crack, b, K₀) -> Tens{2,3}

Ribbon: ``\\mathbf R = \\tfrac{2}{\\pi}\\,b\\,\\hat{\\mathbf w}
\\otimes\\hat{\\mathbf w}``.  For iso / aligned-TI matrices
``\\hat{\\mathbf w} = \\hat{\\mathbf n}``.
"""
function _resistivity_from_b_ribbon(crack::RibbonCrack, b, K₀)
    ŵ = _effective_crack_direction(crack, K₀)
    T = eltype(ŵ)
    return (T(2) / T(π) * b) * (ŵ ⊗ ŵ)
end

_resistivity_from_b(crack::EllipticCrack, b, K₀) = _resistivity_from_b_elliptic(crack, b, K₀)
_resistivity_from_b(crack::RibbonCrack, b, K₀) = _resistivity_from_b_ribbon(crack, b, K₀)

# =============================================================================
#  Budiansky density helpers — dilute compliance / resistivity corrections.
# =============================================================================

"""
    delta_compliance(crack, H, ε) -> Tens{4,3}

Dilute compliance correction ``\\Delta\\mathbb S`` of a family of
identical parallel cracks of Budiansky density ``\\varepsilon`` from the
size-independent contribution tensor ``\\mathbb H``:

- Elliptic: ``\\Delta\\mathbb S = \\tfrac{4\\pi}{3}\\,\\varepsilon^{3\\mathrm d}\\,\\mathbb H``
  with ``\\varepsilon^{3\\mathrm d} = N a b^{2}``.
- Ribbon:   ``\\Delta\\mathbb S = \\pi\\,\\varepsilon^{2\\mathrm d}\\,\\mathbb H``
  with ``\\varepsilon^{2\\mathrm d} = N b^{2}``.

See [Budiansky & O'Connell 1976](@cite budiansky1976),
[Sevostianov & Kachanov 2002](@cite sevostianov2002).
"""
delta_compliance(crack::EllipticCrack, H, ε) = (4 * one(eltype(H)) * π / 3) * ε * H
delta_compliance(crack::RibbonCrack, H, ε) = (one(eltype(H)) * π) * ε * H

"""
    delta_resistivity(crack, R, ε) -> Tens{2,3}

Dilute resistivity correction ``\\Delta\\mathbf R`` of a family of
identical parallel cracks of Budiansky density ``\\varepsilon`` from the
size-independent contribution tensor ``\\mathbf R``:

- Elliptic: ``\\Delta\\mathbf R = \\tfrac{4\\pi}{3}\\,\\varepsilon^{3\\mathrm d}\\,\\mathbf R``.
- Ribbon:   ``\\Delta\\mathbf R = \\pi\\,\\varepsilon^{2\\mathrm d}\\,\\mathbf R``.
"""
delta_resistivity(crack::EllipticCrack, R, ε) = (4 * one(eltype(R)) * π / 3) * ε * R
delta_resistivity(crack::RibbonCrack, R, ε) = (one(eltype(R)) * π) * ε * R

# =============================================================================
#  Stiffness / conductivity contribution for cracks (API symmetry with
#  ellipsoids).  For a flat crack the dilute expansion of `inv(S₀+ΔS) - C₀`
#  at first order in the density ε gives N_crack = -C₀ : H : C₀, so both
#  contribution flavours are related by a simple ± C₀ : (·) : C₀ mapping.
# =============================================================================

"""
    stiffness_contribution(crack, C₀; kw...) -> Tens{4,3}

Size-independent **crack stiffness contribution tensor**
``\\mathbb N = -\\mathbb C_0 : \\mathbb H : \\mathbb C_0``, where
``\\mathbb H`` is the crack compliance contribution tensor
([`compliance_contribution`](@ref)).  Provided for API symmetry with
solid inclusions; the associated dilute correction is
``\\Delta\\mathbb C = (4\\pi/3)\\,\\varepsilon^{3\\mathrm d}\\,\\mathbb N``
(elliptic) or ``\\pi\\,\\varepsilon^{2\\mathrm d}\\,\\mathbb N`` (ribbon),
assembled by [`delta_stiffness`](@ref)`(crack, N, ε)`.
"""
function MFH_Core.stiffness_contribution(
        crack::MFH_Core.AbstractCrack,
        C₀::TensND.AbstractTens{4, 3};
        kw...
    )
    H = compliance_contribution(crack, C₀; kw...)
    return -(C₀ ⊡ H ⊡ C₀)
end

"""
    conductivity_contribution(crack, K₀; kw...) -> Tens{2,3}

Size-independent **crack conductivity contribution tensor**
``\\mathbf N_K = -\\mathbf K_0 \\cdot \\mathbf R \\cdot \\mathbf K_0``.
"""
function MFH_Core.conductivity_contribution(
        crack::MFH_Core.AbstractCrack,
        K₀::TensND.AbstractTens{2, 3};
        kw...
    )
    R = compliance_contribution(crack, K₀; kw...)
    return -(K₀ ⋅ R ⋅ K₀)
end

"""
    delta_stiffness(crack, N, ε) -> Tens{4,3}

Dilute stiffness correction ``\\Delta\\mathbb C`` from the
size-independent crack contribution tensor ``\\mathbb N`` and the
Budiansky density ``\\varepsilon``:

- Elliptic: ``\\Delta\\mathbb C = \\tfrac{4\\pi}{3}\\,\\varepsilon^{3\\mathrm d}\\,\\mathbb N``.
- Ribbon:   ``\\Delta\\mathbb C = \\pi\\,\\varepsilon^{2\\mathrm d}\\,\\mathbb N``.
"""
MFH_Core.delta_stiffness(crack::EllipticCrack, N, ε) = (4 * one(eltype(N)) * π / 3) * ε * N
MFH_Core.delta_stiffness(crack::RibbonCrack, N, ε) = (one(eltype(N)) * π) * ε * N

"""
    delta_conductivity(crack, N_K, ε) -> Tens{2,3}

Dilute conductivity correction from the crack contribution tensor
``\\mathbf N_K`` and the Budiansky density, with the same prefactors as
[`delta_stiffness`](@ref).
"""
MFH_Core.delta_conductivity(crack::EllipticCrack, N, ε) = (4 * one(eltype(N)) * π / 3) * ε * N
MFH_Core.delta_conductivity(crack::RibbonCrack, N, ε) = (one(eltype(N)) * π) * ε * N
