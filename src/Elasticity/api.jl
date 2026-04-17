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
(elasticity) or a 2nd-order conductivity tensor — dispatch handles
both.

See the package documentation for the full dispatch table and a
discussion of return types vs matrix symmetry class.
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
