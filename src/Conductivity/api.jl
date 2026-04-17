# =============================================================================
#  api.jl — add `_kernel` methods for 2nd-order Hill tensor (conductivity).
#
#  `hill_tensor` itself is defined in the `Elasticity` sub-module — the
#  dispatch on the 2nd-order case goes through the same generic entry
#  point, so we only need to register new `_kernel` methods here.
# =============================================================================

# ── 2nd-order, 3D (conductivity) ─────────────────────────────────────────────

function Elasticity._kernel(ell::Elasticity.Ellipsoid{3},
                            K₀::TensND.TensISO{2,3},
                            ::MFH_Core.Analytical; kw...)
    return _hill_order2_3d_iso(ell, K₀)
end

function Elasticity._kernel(ell::Elasticity.Ellipsoid{3},
                            K₀::TensND.AbstractTens{2,3},
                            ::MFH_Core.Analytical; kw...)
    return _hill_order2_3d_aniso(ell, K₀)
end

# ── 2nd-order, 2D (conductivity) ─────────────────────────────────────────────

function Elasticity._kernel(ell::Elasticity.Ellipsoid{2},
                            K₀::TensND.TensISO{2,2},
                            ::MFH_Core.Analytical; kw...)
    return _hill_order2_2d_iso(ell, K₀)
end

function Elasticity._kernel(ell::Elasticity.Ellipsoid{2},
                            K₀::TensND.AbstractTens{2,2},
                            ::MFH_Core.Analytical; kw...)
    return _hill_order2_2d(ell, K₀)
end
