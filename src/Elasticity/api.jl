# =============================================================================
#  api.jl
#
#  Public entry point `hill_tensor(ell, C₀; method, kw...)` and the
#  `_kernel` method table for elasticity and 2nd-order Hill tensors
#  (4th-order case handled in this module; 2nd-order case handled in
#  `Conductivity`).  Dispatch is delegated to `Core._resolve_algo`.
# =============================================================================

"""
    hill_tensor(ell, C₀; method=:auto, abstol=1e-8, reltol=1e-6, maxiters=1_000_000)
        → AbstractTens

Hill polarisation tensor **P** for an ellipsoidal inclusion `ell`
embedded in a reference medium `C₀`.  `C₀` can be a 4th-order stiffness
(elasticity) or a 2nd-order conductivity tensor — dispatch selects
the appropriate formulation automatically.

The general expression of the elastic polarisation tensor is
([Willis 1977](@cite willis1977), [Mura 1987](@cite mura1987)):

```
P(A, C) = (det A)/(4π) ∫_{|ξ|=1} ξ ⊗ˢ (ξ·C·ξ)⁻¹ ⊗ˢ ξ / ‖A·ξ‖³ dS_ξ
```

The isotropic case (`C₀::TensISO`) is evaluated analytically; the
anisotropic case uses the Cauchy-residue reduction of
[Masson 2008](@cite masson2008) (trait `Residue`) or the DECUHR
adaptive cubature of [Espelid & Genz 1994](@cite espelid1994)
(trait `DECUHR`). See the `Hill polarisation tensors` theory page
for the full dispatch table and return types.
"""
function hill_tensor(
        ell::AbstractEllipsoidalInclusion,
        C₀::TensND.AbstractTens;
        method::Symbol = :auto,
        abstol::Float64 = 1.0e-8,
        reltol::Float64 = 1.0e-6,
        maxiters::Int = 1_000_000
    )
    algo = MFH_Core._resolve_algo(Val(method), ell, C₀)
    return _kernel(ell, C₀, algo; abstol = abstol, reltol = reltol, maxiters = maxiters)
end

# ── 4th-order, 3D ────────────────────────────────────────────────────────────

_kernel(ell::Ellipsoid{3}, C₀::TensND.TensISO{4, 3}, ::MFH_Core.Analytical; kw...) =
    _hill_3d_iso(ell, C₀)

_kernel(ell::Ellipsoid{3}, C₀::TensND.AbstractTens{4, 3}, ::MFH_Core.Residue; kw...) =
    _hill_3d_aniso_residue(
    ell, C₀;
    abstol = get(kw, :abstol, 1.0e-8),
    reltol = get(kw, :reltol, 1.0e-6),
    maxiters = get(kw, :maxiters, 1_000_000)
)

_kernel(ell::Ellipsoid{3}, C₀::TensND.AbstractTens{4, 3}, ::MFH_Core.DECUHR; kw...) =
    _hill_3d_aniso_decuhr(
    ell, C₀;
    abstol = get(kw, :abstol, 1.0e-8),
    reltol = get(kw, :reltol, 1.0e-6),
    maxiters = get(kw, :maxiters, 1_000_000)
)

_kernel(ell::Ellipsoid{3}, C₀::TensND.AbstractTens{4, 3}, ::MFH_Core.NestedQuadGK; kw...) =
    _hill_3d_aniso_nestedquadgk(
    ell, C₀;
    abstol = get(kw, :abstol, 1.0e-8),
    reltol = get(kw, :reltol, 1.0e-6),
    maxiters = get(kw, :maxiters, 1_000_000)
)

# ── 4th-order, 3D — infinite cylinder ────────────────────────────────────────

_kernel(cyl::Cylinder, C₀::TensND.TensISO{4, 3}, ::MFH_Core.Analytical; kw...) =
    _hill_3d_iso(cyl, C₀)

_kernel(cyl::Cylinder, C₀::TensND.AbstractTens{4, 3}, ::MFH_Core.CylinderQuadrature; kw...) =
    _hill_3d_cylinder_aniso(
    cyl, C₀;
    abstol = get(kw, :abstol, 1.0e-8),
    reltol = get(kw, :reltol, 1.0e-6),
    maxiters = get(kw, :maxiters, 1_000_000)
)

# ── 4th-order, 2D ────────────────────────────────────────────────────────────

_kernel(ell::Ellipsoid{2}, C₀::TensND.TensISO{4, 2}, ::MFH_Core.Analytical; kw...) =
    _hill_2d_iso(ell, C₀)

_kernel(ell::Ellipsoid{2}, C₀::TensND.AbstractTens{4, 2}, ::MFH_Core.Analytical; kw...) =
    _hill_2d_aniso(
    ell, C₀;
    abstol = get(kw, :abstol, 1.0e-8),
    reltol = get(kw, :reltol, 1.0e-6),
    maxiters = get(kw, :maxiters, 1_000_000)
)

# ── Eshelby tensor (4th order) — S = P : C₀ ──────────────────────────────────

"""
    eshelby_tensor(incl::AbstractEllipsoidalInclusion, C₀::TensND.AbstractTens{4}; kw...)

4th-order Eshelby tensor ``\\mathbb S = \\mathbb P : \\mathbb C_0``
of an ellipsoidal inclusion `incl` embedded in a matrix of stiffness
`C₀`. Thin wrapper around [`hill_tensor`](@ref) followed by the double
contraction with `C₀`.
"""
MFH_Core.eshelby_tensor(
    incl::AbstractEllipsoidalInclusion, C₀::TensND.AbstractTens{4, 3};
    kw...
) = hill_tensor(incl, C₀; kw...) ⊡ C₀

MFH_Core.eshelby_tensor(
    incl::AbstractEllipsoidalInclusion, C₀::TensND.AbstractTens{4, 2};
    kw...
) = hill_tensor(incl, C₀; kw...) ⊡ C₀
