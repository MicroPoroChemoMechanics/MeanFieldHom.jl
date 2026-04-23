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
        _ti_aligned(C₀, crack_basis(crack)) && return MFH_Core.Analytical()
        m === :decuhr && return MFH_Core.DECUHR()
        m === :nestedquadgk && return MFH_Core.NestedQuadGK()
        return MFH_Core.Residue()
    end
end

# ── Public API ──────────────────────────────────────────────────────────────

"""
    cod_tensor(crack, C₀; method=:auto, abstol=1e-8, reltol=1e-6, maxiters=100_000)
        -> Tens{2,3}

Size-independent crack-opening-displacement (COD) tensor
``\\mathbf B`` defined from the average displacement jump on the crack
surface through

```
(1/S) ∫_S [u] dS = b · B · (Σ · n̂),
```

where ``b`` is the semi-minor in-plane semi-axis (``b\\ge c\\to 0``).
``\\mathbf B`` factors the crack compliance tensor as

- Elliptic: ``\\mathbb H = \\tfrac{3}{4}\\,\\hat{\\mathbf n}
  \\stackrel{s}{\\otimes}\\mathbf B\\stackrel{s}{\\otimes}\\hat{\\mathbf n}``;
- Ribbon:   ``\\mathbb H = \\tfrac{2}{\\pi}\\,\\hat{\\mathbf n}
  \\stackrel{s}{\\otimes}\\mathbf B\\stackrel{s}{\\otimes}\\hat{\\mathbf n}``;

with ``\\mathbb H = \\lim_{c/b\\to 0}(c/b)\\,\\mathbb Q^{-1}`` and
``\\mathbb Q = \\mathbb C - \\mathbb C:\\mathbb P:\\mathbb C`` the
second Hill tensor
([Kachanov 1992](@cite kachanov1992),
 [Sevostianov & Kachanov 2002](@cite sevostianov2002),
 [Barthélémy 2021](@cite barthelemyIJES2021)).  The elliptic and ribbon
factorisations are related by ``\\mathbf B^{\\mathcal R} =
\\tfrac{3\\pi}{8}\\,\\lim_{\\eta\\to 0}\\mathbf B^{\\mathcal E}``.

For isotropic or aligned-TI matrices the kernel is analytical
([Hoenig 1978](@cite hoenig1978),
 [Kanaun & Levin 2009](@cite kanaun2009));
for arbitrarily anisotropic matrices the limit ``c/b\\to 0`` is
resolved numerically through the first-order Taylor term of the Hill
tensor ([Barthélémy 2009](@cite barthelemyIJSS2009)).

Alias: [`B_tensor`](@ref).
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

# Order-2 (conductivity) — thermal COD scalar
"""
    cod_tensor(crack, K₀::AbstractTens{2,3}; method=:auto, kw...) -> Real

Size-independent **thermal crack-opening-displacement scalar** ``b``
for a flat crack in a conductor of 2nd-order conductivity tensor
``\\mathbf K_0``.  Analogue of the elasticity COD tensor: in the 2nd-
order problem, the temperature jump across the crack is scalar and
only the normal component of the heat flux drives it, so a single
scalar captures the full crack flexibility.  The associated
size-independent resistivity contribution ``\\mathbf R = `` [`compliance_contribution`](@ref)`(crack, K₀)` is

```
R = (3/4) b · ŵ ⊗ ŵ   (elliptic)
R = (2/π) b · ŵ ⊗ ŵ   (ribbon)
```

with ``\\hat{\\mathbf w}\\parallel\\mathbf K_0^{-1/2}\\hat{\\mathbf n}``
(reduces to ``\\hat{\\mathbf n}`` for iso / aligned-TI matrices).
Apply [`delta_resistivity`](@ref) to recover the dilute resistivity
correction ``\\Delta\\mathbf R = (4\\pi/3)\\varepsilon^{3\\mathrm d}\\mathbf R``
(elliptic) or ``\\Delta\\mathbf R = \\pi\\,\\varepsilon^{2\\mathrm d}\\mathbf R``
(ribbon).

Closed-form derivation via the square-root change-of-variable of
[Giraud & Gruescu 2019](@cite giraudMOM2019), in the framework of
[Sevostianov & Kachanov 2002](@cite sevostianov2002) for the rank-1
factorisation ``\\mathbf R \\propto \\hat{\\mathbf w}\\otimes\\hat{\\mathbf w}``.
See [`docs/src/theory/thermal_cracks.md`](docs/src/theory/thermal_cracks.md)
for the mathematical details and the elasticity ↔ conductivity
correspondence table.
"""
function cod_tensor(
        crack::MFH_Core.AbstractCrack,
        K₀::TensND.AbstractTens{2, 3};
        method::Symbol = :auto,
        kw...
    )
    return _cod_thermal(crack, K₀)
end

# Analytical dispatch: iso vs aniso.
_cod_thermal(crack::EllipticCrack, K₀::TensND.TensISO{2, 3}) =
    _cod_iso_ellipse_thermal(crack, MFH_Core.extract_iso_conductivity(K₀))

_cod_thermal(crack::RibbonCrack, K₀::TensND.TensISO{2, 3}) =
    _cod_iso_ribbon_thermal(crack, MFH_Core.extract_iso_conductivity(K₀))

_cod_thermal(crack::EllipticCrack, K₀::TensND.AbstractTens{2, 3}) =
    _cod_aniso_ellipse_thermal(crack, K₀)

_cod_thermal(crack::RibbonCrack, K₀::TensND.AbstractTens{2, 3}) =
    _cod_aniso_ribbon_thermal(crack, K₀)

"""
    compliance_contribution(crack, K₀::AbstractTens{2,3}; kw...) -> Tens{2,3}

Size-independent **crack resistivity contribution tensor** ``\\mathbf R``
(thermal analogue of the elasticity [`compliance_contribution`](@ref)):

- Elliptic crack:  ``\\mathbf R = \\tfrac{3}{4}\\,b\\,
  \\hat{\\mathbf w}\\otimes\\hat{\\mathbf w}``.
- Ribbon crack:    ``\\mathbf R = \\tfrac{2}{\\pi}\\,b\\,
  \\hat{\\mathbf w}\\otimes\\hat{\\mathbf w}``.

with ``b = `` [`cod_tensor`](@ref)`(crack, K₀)` the thermal COD scalar
and ``\\hat{\\mathbf w}\\parallel\\mathbf K_0^{-1/2}\\hat{\\mathbf n}``
(reduces to ``\\hat{\\mathbf n}`` for iso / aligned-TI matrices).
Apply [`delta_resistivity`](@ref)`(crack, R, ε)` to obtain the dilute
resistivity correction ``\\Delta\\mathbf R``.
"""
function compliance_contribution(
        crack::MFH_Core.AbstractCrack,
        K₀::TensND.AbstractTens{2, 3};
        kw...
    )
    b = cod_tensor(crack, K₀; kw...)
    return _resistivity_from_b(crack, b, K₀)
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

# NestedQuadGK — anisotropic matrix (historical nested-1D-QuadGK path,
# formerly shipped under the DECUHR name).
function _kernel(crack::EllipticCrack, C₀::TensND.AbstractTens{4, 3}, ::MFH_Core.NestedQuadGK; kw...)
    return _cod_elliptic_nestedquadgk_direct(
        crack, C₀;
        abstol = get(kw, :abstol, 1.0e-8),
        reltol = get(kw, :reltol, 1.0e-6),
        maxiters = get(kw, :maxiters, 100_000)
    )
end

function _kernel(crack::RibbonCrack, C₀::TensND.AbstractTens{4, 3}, ::MFH_Core.NestedQuadGK; kw...)
    return _cod_ribbon_numerical(
        crack, C₀, _nestedquadgk_backend;
        abstol = get(kw, :abstol, 1.0e-8),
        reltol = get(kw, :reltol, 1.0e-6),
        maxiters = get(kw, :maxiters, 100_000)
    )
end
