# =============================================================================
#  symmetrize.jl — orientation-distribution treatment of a phase's tensors.
#
#  TWO distinct mechanisms, mirroring echoes (tensor_symmetry.h vs tensor_ti.h);
#  they must never be conflated :
#
#  1. EXACT rotation-group averaging (runtime, inside scheme kernels) —
#     `_apply_symmetrize(t, sym)` delegates to `Core.isotropify` /
#     `Core.transverse_isotropify`.  Valid for minor-symmetric tensors with
#     or without major symmetry ; the TI average returns the full 8-parameter
#     axially-invariant tensor (`TensND.TensTI{4,T,8}`), preserving ℓ₃ ≠ ℓ₄
#     and the antisymmetric azimuthal couplings (ℓ₇, ℓ₈) that concentration
#     tensors generally carry.  At 2nd order the antisymmetric in-plane part
#     is preserved (`TensTI{2,T,3}`).
#
#  2. BEST-FIT projection (reporting / parameter extraction ONLY) —
#     `best_fit_ti(t, axis)` (→ major-symmetric `TensTI{4,T,5}`) and
#     `best_fit_iso(t)`.  This is the orthogonal projection onto the
#     symmetric Walpole span, the analog of echoes' `.paramsym(sym=TI)`.
#     Do NOT use it inside scheme kernels : it silently drops the
#     non-major-symmetric content of concentration tensors.
#
#  References : Walpole (1981) for the TI basis ; echoes
#  `tensor_symmetry.h` for the exact azimuthal-average closed form.
# =============================================================================

"""
    _apply_symmetrize(t::AbstractTens, sym::AbstractSymmetrize) -> AbstractTens

Apply the **exact** orientation average declared by `sym` to `t`.

`NoSymmetrize` is a passthrough ; `IsoSymmetrize` returns the SO(3) average
(`TensISO`) ; `TISymmetrize(axis)` returns the azimuthal average about
`axis` (`TensTI{4,T,8}` / `TensTI{2,T,3}`, non-major-symmetric content
preserved). Implemented for tensor orders 4 (elasticity) and 2
(conductivity).
"""
_apply_symmetrize(t::TensND.AbstractTens, ::NoSymmetrize) = t

_apply_symmetrize(t::TensND.AbstractTens{4, 3}, ::IsoSymmetrize) =
    MFH_Core.isotropify(t)
_apply_symmetrize(t::TensND.AbstractTens{2, 3}, ::IsoSymmetrize) =
    MFH_Core.isotropify(t)

_apply_symmetrize(t::TensND.AbstractTens{4, 3}, sym::TISymmetrize) =
    MFH_Core.transverse_isotropify(t, sym.axis)
_apply_symmetrize(t::TensND.AbstractTens{2, 3}, sym::TISymmetrize) =
    MFH_Core.transverse_isotropify(t, sym.axis)

# =============================================================================
#  Reference-medium projection for the localization-tensor computation.
# =============================================================================

"""
    _project_matrix(P₀::AbstractTens, sym::AbstractSymmetrize) -> AbstractTens

Project the reference medium `P₀` before computing the localization tensor
of a phase declaring the orientation distribution `sym`.

- `NoSymmetrize` : passthrough.
- `IsoSymmetrize` : isotropic average of `P₀` (the inclusion's Hill tensor
  in an isotropic matrix always has an analytical branch).
- `TISymmetrize` : controlled by `sym.matrix_projection` :
    * `:iso` (default) — isotropic average of `P₀`.  Approximation whenever
      `P₀` is not isotropic ; exact at the isotropic fixed point of the SC
      iteration (where the reference converges to its isotropic average).
      Rationale : an inclusion family at polar angle θ ≠ 0 from the
      symmetrize axis is *not* coaxial with a TI reference, so the
      TI-coaxial analytical Hill branch does not apply ; the isotropic
      projection guarantees an analytical, ForwardDiff-compatible branch
      for every orientation.
    * `:none` — no projection ; non-coaxial anisotropic references route
      through the general-anisotropy `NestedQuadGK` Hill branch
      (ForwardDiff-compatible, quadrature-priced).
    * `:ti` — best-fit TI projection of `P₀` about `sym.axis`
      (`TensTI{4,T,5}`) ; only meaningful when the phase's inclusions are
      coaxial with the axis.

The phase contribution is in all cases still exactly averaged by the
outgoing `_apply_symmetrize`, so the outer symmetry semantics are
preserved.
"""
_project_matrix(P₀::TensND.AbstractTens, ::NoSymmetrize) = P₀
_project_matrix(P₀::TensND.AbstractTens, ::IsoSymmetrize) = MFH_Core.isotropify(P₀)
function _project_matrix(P₀::TensND.AbstractTens, sym::TISymmetrize)
    mp = sym.matrix_projection
    mp === :none && return P₀
    mp === :ti && return best_fit_ti(P₀, sym.axis)
    return MFH_Core.isotropify(P₀)
end

# =============================================================================
#  Best-fit projections (reporting / parameter extraction only)
# =============================================================================

"""
    best_fit_iso(t::AbstractTens{4,3}) -> TensND.TensISO{4}
    best_fit_iso(t::AbstractTens{2,3}) -> TensND.TensISO{2}

Orthogonal (Frobenius) projection of `t` onto the isotropic basis — thin
wrapper over [`TensND.proj_tens`](@ref)`(Val(:ISO), t)`, TensND's canonical
"paramsym"-style extraction. For minor-symmetric tensors this coincides with
the exact SO(3) average [`Core.isotropify`](@ref) (the isotropic subspace is
`{𝕁, 𝕂}` either way).
"""
best_fit_iso(t::TensND.AbstractTens) = TensND.proj_tens(Val(:ISO), t)[1]

"""
    best_fit_ti(t::AbstractTens{4,3}, axis) -> TensND.TensTI{4,T,5}
    best_fit_ti(t::AbstractTens{2,3}, axis) -> TensND.TensTI{2,T,2}

Orthogonal (Frobenius) projection of `t` onto the **major-symmetric** TI
(Walpole) span about `axis` — thin wrapper over [`TensND.proj_tens`](@ref)`(
Val(:TI), t, axis)`, the analog of echoes' `.paramsym(sym=TI)` parameter
extraction. Numerically identical to the previous in-house implementation
(exact azimuthal average then forced major symmetry), verified to ~1e-11.

!!! warning
    This is a reporting utility, NOT the orientation average : it forces
    major symmetry (ℓ₃+ℓ₄)/2 and drops the antisymmetric azimuthal
    couplings.  Inside scheme kernels use `_apply_symmetrize` /
    [`Core.transverse_isotropify`](@ref) instead.
"""
best_fit_ti(t::TensND.AbstractTens{4, 3}, axis) = TensND.proj_tens(Val(:TI), t, axis)[1]
best_fit_ti(t::TensND.AbstractTens{2, 3}, axis) = TensND.proj_tens(Val(:TI), t, axis)[1]

"""
    best_fit_ortho(t::AbstractTens{4,3}, frame) -> TensND.TensOrtho
    best_fit_ortho(t::AbstractTens{2,3}, frame) -> Matrix

Orthogonal (Frobenius) projection of `t` onto the orthotropic span in the
given material `frame` — thin wrapper over [`TensND.proj_tens`](@ref)`(
Val(:ORTHO), t, frame)`, the analog of echoes' `.paramsym(sym=ORTHO)`.
There was previously no orthotropic parameter extraction in MeanFieldHom.jl;
this closes that gap using the TI/ORTHO projection machinery already tested
in TensND (`test/test_tens_projection.jl`).
"""
best_fit_ortho(t::TensND.AbstractTens{4, 3}, frame) = TensND.proj_tens(Val(:ORTHO), t, frame)[1]
best_fit_ortho(t::TensND.AbstractTens{2, 3}, frame) = TensND.proj_tens(Val(:ORTHO), t, frame)[1]
