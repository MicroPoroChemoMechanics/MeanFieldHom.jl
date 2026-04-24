# =============================================================================
#  dispatch.jl
#
#  Centralised `_resolve_algo(Val(method), incl, C₀)` dispatch — single
#  place where the symbol `:auto` / `:residues` / `:decuhr` is translated
#  into an [`AbstractAlgorithm`](@ref) instance, taking into account the
#  material symmetry class (TensND type of `C₀`) *and* the inclusion
#  class (`AbstractEllipsoidalInclusion` vs `AbstractCrack`).
#
#  Rules
#  -----
#  - Isotropic matrix (`TensISO`) of any order, any dimension  → Analytical
#  - 2D stiffness / conductivity (any inclusion)               → Analytical
#  - 3D conductivity (2nd-order tensor)                        → Analytical
#  - 3D anisotropic elasticity (4th-order `AbstractTens{4,3}`) :
#      * `AbstractEllipsoidalInclusion`, `:auto` or `:residues` → Residue
#      * `AbstractEllipsoidalInclusion`, `:decuhr`             → DECUHR
#      * `AbstractCrack`, TI + aligned with n̂                  → Analytical
#      * `AbstractCrack`, `:auto` or `:residues`                → Residue
#      * `AbstractCrack`, `:decuhr`                            → DECUHR
#
#  The rules that depend on the *inclusion* class are injected at the end
#  of the file so that they can refer to the abstract types defined in
#  `abstractions.jl` without introducing cycles.
#
#  Every dispatch method lists BOTH the inclusion and the TensND type as
#  method parameters.  Matching all three axes (method symbol × inclusion
#  class × material symmetry) in a single signature avoids the
#  classical "diamond" ambiguity between the inclusion-refined and the
#  symmetry-refined specialisations.
# =============================================================================

# ─── TensISO (any kind) → Analytical — most specific on symmetry ───────────

_resolve_algo(::Val, ::AbstractInclusion, ::TensND.TensISO) = Analytical()
_resolve_algo(::Val{:auto}, ::AbstractInclusion, ::TensND.TensISO) = Analytical()
_resolve_algo(::Val{:residues}, ::AbstractInclusion, ::TensND.TensISO) = Analytical()
_resolve_algo(::Val{:decuhr}, ::AbstractInclusion, ::TensND.TensISO) = Analytical()

# Specific inclusion + TensISO — keep the same rule (needed to avoid
# ambiguity with the inclusion-refined 3D-aniso methods below).
_resolve_algo(::Val, ::AbstractEllipsoidalInclusion, ::TensND.TensISO) = Analytical()
_resolve_algo(::Val{:auto}, ::AbstractEllipsoidalInclusion, ::TensND.TensISO) = Analytical()
_resolve_algo(::Val{:residues}, ::AbstractEllipsoidalInclusion, ::TensND.TensISO) = Analytical()
_resolve_algo(::Val{:decuhr}, ::AbstractEllipsoidalInclusion, ::TensND.TensISO) = Analytical()

_resolve_algo(::Val, ::AbstractCrack, ::TensND.TensISO) = Analytical()
_resolve_algo(::Val{:auto}, ::AbstractCrack, ::TensND.TensISO) = Analytical()
_resolve_algo(::Val{:residues}, ::AbstractCrack, ::TensND.TensISO) = Analytical()
_resolve_algo(::Val{:decuhr}, ::AbstractCrack, ::TensND.TensISO) = Analytical()

# ─── 2D elasticity → Analytical ──────────────────────────────────────────────

_resolve_algo(::Val, ::AbstractInclusion, ::TensND.AbstractTens{4, 2}) = Analytical()
_resolve_algo(::Val, ::AbstractEllipsoidalInclusion, ::TensND.AbstractTens{4, 2}) = Analytical()

# ─── Conductivity (2nd-order) → Analytical ───────────────────────────────────

_resolve_algo(::Val, ::AbstractInclusion, ::TensND.AbstractTens{2, 3}) = Analytical()
_resolve_algo(::Val, ::AbstractEllipsoidalInclusion, ::TensND.AbstractTens{2, 3}) = Analytical()
_resolve_algo(::Val, ::AbstractInclusion, ::TensND.AbstractTens{2, 2}) = Analytical()
_resolve_algo(::Val, ::AbstractEllipsoidalInclusion, ::TensND.AbstractTens{2, 2}) = Analytical()

# ─── 3D anisotropic elasticity ───────────────────────────────────────────────

# Ellipsoidal inclusions — 3D anisotropic default is the residue algorithm
_resolve_algo(::Val{:auto}, ::AbstractEllipsoidalInclusion, ::TensND.AbstractTens{4, 3}) = Residue()
_resolve_algo(::Val{:residues}, ::AbstractEllipsoidalInclusion, ::TensND.AbstractTens{4, 3}) = Residue()
_resolve_algo(::Val{:decuhr}, ::AbstractEllipsoidalInclusion, ::TensND.AbstractTens{4, 3}) = DECUHR()
_resolve_algo(::Val{:nestedquadgk}, ::AbstractEllipsoidalInclusion, ::TensND.AbstractTens{4, 3}) = NestedQuadGK()

# ─── Cylinder dispatch — anisotropic elasticity (3D) routes to the dedicated
# 1D transverse-plane quadrature.  The residue algorithm is not applicable
# to a cylinder (acoustic polynomial degenerates), so `:residues` silently
# falls back to `CylinderQuadrature`.  The rules are injected by the
# Elasticity sub-module (see `Elasticity.jl`) to avoid a Core→Elasticity
# dependency — this file only declares the infrastructure.

# Generic inclusion fallback (also used by `AbstractCrack` before the
# TI-aligned refinement injected from the `Cracks` sub-module).
_resolve_algo(::Val{:auto}, ::AbstractInclusion, ::TensND.AbstractTens{4, 3}) = Residue()
_resolve_algo(::Val{:residues}, ::AbstractInclusion, ::TensND.AbstractTens{4, 3}) = Residue()
_resolve_algo(::Val{:decuhr}, ::AbstractInclusion, ::TensND.AbstractTens{4, 3}) = DECUHR()
_resolve_algo(::Val{:nestedquadgk}, ::AbstractInclusion, ::TensND.AbstractTens{4, 3}) = NestedQuadGK()

# Plain catch-all for the cases where a sub-module passes `method` symbols
# we don't know about (e.g. a user extension): default to Analytical.
_resolve_algo(::Val, ::AbstractInclusion, ::TensND.AbstractTens) = Analytical()

# ─── Public helper ───────────────────────────────────────────────────────────

"""
    _resolve_algo(Val(method), incl, C₀) -> AbstractAlgorithm

Translate a method symbol (`:auto`, `:residues`, `:decuhr`) into the
concrete [`AbstractAlgorithm`](@ref) instance, taking the inclusion
class and the TensND symmetry class of `C₀` into account.  Centralised
here so that every high-level entry point (`hill_tensor`,
`cod_tensor`, `compliance_contribution`, `sif`, `dif`) shares a single
truth table — new dispatch rules are added as new methods of this
function.
"""
_resolve_algo
