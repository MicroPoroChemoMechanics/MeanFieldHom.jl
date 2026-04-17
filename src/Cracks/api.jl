# =============================================================================
#  api.jl — public entry point `cod_tensor` and `_kernel` method table
#  for flat cracks.  Dispatch via `Core._resolve_algo`.
# =============================================================================

"""
    _ti_aligned(C₀::TensWalpole, ℬ_crack) -> Bool

Return `true` when the axis of transverse isotropy stored in `C₀` is
parallel to the third axis (the crack normal) of the crack-local basis
`ℬ_crack`.
"""
function _ti_aligned(C₀, ℬ_crack::TensND.AbstractBasis)
    axis_C = TensND.components_canon(TensND.tensbasis(TensND.getbasis(C₀), 3))
    axis_n = TensND.components_canon(TensND.tensbasis(ℬ_crack, 3))
    d = abs(dot(axis_C, axis_n))
    return isapprox(d, 1.0; atol = 1.0e-10)
end

# TI-aligned dispatch rules — refine Core dispatch for `AbstractCrack` + TensWalpole
if isdefined(TensND, :TensWalpole)
    @eval function MFH_Core._resolve_algo(::Val{m}, crack::MFH_Core.AbstractCrack, C₀::TensND.TensWalpole) where {m}
        return _ti_aligned(C₀, crack_basis(crack)) ? MFH_Core.Analytical() :
            m === :decuhr ? MFH_Core.DECUHR() : MFH_Core.Residue()
    end
end

# ── Public API ──────────────────────────────────────────────────────────────

"""
    cod_tensor(crack, C₀; method=:auto, abstol=1e-8, reltol=1e-6, maxiters=100_000)
        -> Tens{2,3}

Size-independent crack opening displacement (COD) tensor ``\\mathbf B``.
"""
function cod_tensor(
        crack::MFH_Core.AbstractCrack,
        C₀::TensND.AbstractTens{4, 3};
        method::Symbol = :auto,
        abstol::Float64 = 1.0e-8,
        reltol::Float64 = 1.0e-6,
        maxiters::Int = 100_000
    )
    algo = MFH_Core._resolve_algo(Val(method), crack, C₀)
    return _kernel(crack, C₀, algo; abstol = abstol, reltol = reltol, maxiters = maxiters)
end

const B_tensor = cod_tensor

# Order-2 (conductivity) — stub
"""
    cod_tensor(crack, K₀::AbstractTens{2,3}; kw...)

Not implemented yet — see `docs/src/developer/roadmap.md`.
"""
function cod_tensor(
        crack::MFH_Core.AbstractCrack,
        K₀::TensND.AbstractTens{2, 3};
        method::Symbol = :auto,
        kw...
    )
    error(
        """
        cod_tensor for 2nd-order conductivity tensors is not implemented yet.
        """
    )
end

function compliance_contribution(
        crack::MFH_Core.AbstractCrack,
        K₀::TensND.AbstractTens{2, 3},
        ε; kw...
    )
    error(
        """
        compliance_contribution for 2nd-order conductivity tensors is not
        implemented yet.
        """
    )
end

# ── Level 2: _kernel methods ─────────────────────────────────────────────────

# Analytical — isotropic matrix
function _kernel(crack::EllipticCrack, C₀::TensND.TensISO{4, 3}, ::MFH_Core.Analytical; kw...)
    E, ν = MFH_Core.extract_iso_moduli(C₀)
    return _cod_iso_ellipse(crack, E, ν)
end

function _kernel(crack::RibbonCrack, C₀::TensND.TensISO{4, 3}, ::MFH_Core.Analytical; kw...)
    E, ν = MFH_Core.extract_iso_moduli(C₀)
    return _cod_iso_ribbon(crack, E, ν)
end

# Analytical — TI matrix aligned with n̂
if isdefined(TensND, :TensWalpole)
    @eval function _kernel(crack::EllipticCrack, C₀::TensND.TensWalpole, ::MFH_Core.Analytical; kw...)
        n̂ = TensND.tensbasis(crack_basis(crack), 3)
        nt = MFH_Core.extract_ti_moduli(C₀, n̂)
        return _cod_ti_ellipse(crack, nt.E, nt.H, nt.ν₁, nt.ν₂, nt.Γ)
    end

    @eval function _kernel(crack::RibbonCrack, C₀::TensND.TensWalpole, ::MFH_Core.Analytical; kw...)
        n̂ = TensND.tensbasis(crack_basis(crack), 3)
        nt = MFH_Core.extract_ti_moduli(C₀, n̂)
        return _cod_ti_ribbon(crack, nt.E, nt.H, nt.ν₁, nt.ν₂, nt.Γ)
    end
end

# Residue — anisotropic matrix
function _kernel(crack::EllipticCrack, C₀::TensND.AbstractTens{4, 3}, ::MFH_Core.Residue; kw...)
    return _cod_elliptic_numerical(
        crack, C₀, _residue_backend;
        abstol = get(kw, :abstol, 1.0e-8),
        reltol = get(kw, :reltol, 1.0e-6),
        maxiters = get(kw, :maxiters, 100_000)
    )
end

function _kernel(crack::RibbonCrack, C₀::TensND.AbstractTens{4, 3}, ::MFH_Core.Residue; kw...)
    return _cod_ribbon_numerical(
        crack, C₀,
        (Carr, ξ, n̂; kw2...) -> _Qnn_star_residue(Carr, ξ, n̂);
        abstol = get(kw, :abstol, 1.0e-8),
        reltol = get(kw, :reltol, 1.0e-6),
        maxiters = get(kw, :maxiters, 100_000)
    )
end

# DECUHR — anisotropic matrix
function _kernel(crack::EllipticCrack, C₀::TensND.AbstractTens{4, 3}, ::MFH_Core.DECUHR; kw...)
    return _cod_elliptic_decuhr_direct(
        crack, C₀;
        abstol = get(kw, :abstol, 1.0e-8),
        reltol = get(kw, :reltol, 1.0e-6),
        maxiters = get(kw, :maxiters, 100_000)
    )
end

function _kernel(crack::RibbonCrack, C₀::TensND.AbstractTens{4, 3}, ::MFH_Core.DECUHR; kw...)
    return _cod_ribbon_numerical(
        crack, C₀, _decuhr_backend;
        abstol = get(kw, :abstol, 1.0e-8),
        reltol = get(kw, :reltol, 1.0e-6),
        maxiters = get(kw, :maxiters, 100_000)
    )
end
