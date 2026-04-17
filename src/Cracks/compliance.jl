# =============================================================================
#  compliance.jl — crack compliance contribution.
# =============================================================================

"""
    compliance_contribution(crack, C₀, ε; method=:auto, kw...) -> Tens{4,3}

Contribution of a family of parallel, identical cracks to the effective
compliance tensor of the cracked medium.
"""
function compliance_contribution(crack::MFH_Core.AbstractCrack, C₀, ε; kw...)
    B = cod_tensor(crack, C₀; kw...)
    return _compliance_from_B(crack, B, ε)
end

"""
    _compliance_from_B_elliptic(crack, B, ε) -> Tens{4,3}
Elliptic: ``ΔS = π ε (\\hat n ⊗ˢ \\mathbf B ⊗ˢ \\hat n)``.
"""
function _compliance_from_B_elliptic(crack::EllipticCrack, B, ε)
    n̂ = TensND.tensbasis(crack_basis(crack), 3)
    return π * ε * (n̂ ⊗ˢ B ⊗ˢ n̂)
end

"""
    _compliance_from_B_ribbon(crack, B, ε) -> Tens{4,3}
Ribbon: ``ΔS = (π/2) ε (\\hat n ⊗ˢ \\mathbf B ⊗ˢ \\hat n)``.
"""
function _compliance_from_B_ribbon(crack::RibbonCrack, B, ε)
    n̂ = TensND.tensbasis(crack_basis(crack), 3)
    return (π / 2) * ε * (n̂ ⊗ˢ B ⊗ˢ n̂)
end

_compliance_from_B(crack::EllipticCrack, B, ε) = _compliance_from_B_elliptic(crack, B, ε)
_compliance_from_B(crack::RibbonCrack, B, ε) = _compliance_from_B_ribbon(crack, B, ε)
