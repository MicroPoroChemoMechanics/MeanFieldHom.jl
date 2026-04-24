# =============================================================================
#  traits.jl
#
#  Holy-style traits shared by every high-level entry point of
#  `MeanFieldHom`:
#
#    - Algorithm traits     (`AbstractAlgorithm`, `Analytical`, `Residue`,
#                            `DECUHR`, `Auto`) — drive `_resolve_algo`.
#    - Material-symmetry    (`MaterialSymmetry`, `IsotropicSym`,
#      traits                `TransverselyIsotropicSym`, …) — classify a
#                            TensND stiffness / conductivity tensor.
#
#  The traits are intentionally plain empty structs: they carry no data,
#  only a type identity used as a dispatch tag.
# =============================================================================

# ─── Algorithm traits ────────────────────────────────────────────────────────

"""
    AbstractAlgorithm

Abstract supertype for computation algorithms exposed through
`MeanFieldHom`. The subtypes are plain singleton structs used as
dispatch tags by the internal `_kernel(...)` machinery.
"""
abstract type AbstractAlgorithm end

"Closed-form / analytical algorithm — used for every isotropic matrix and for every 2D / conductivity case."
struct Analytical <: AbstractAlgorithm end

"Residue-theorem algorithm — polynomial-root based, `Float64` only (incompatible with `ForwardDiff` and `SymPy`)."
struct Residue <: AbstractAlgorithm end

"DECUHR 2D hyper-cubature (Espelid & Genz 1994) via Integrals.jl `DecuhrAlgorithm`. Suitable for general anisotropy in 3D."
struct DECUHR <: AbstractAlgorithm end

"Nested 1D QuadGK cubature (fallback to DECUHR, ForwardDiff-compatible). Historically shipped as `DECUHR` before the split."
struct NestedQuadGK <: AbstractAlgorithm end

"1D QuadGK quadrature dedicated to infinite cylinders (transverse-plane parametrisation ζ(φ) = (0, cos φ / b, sin φ / c)). ForwardDiff-compatible, selected whenever an [`AbstractEllipsoidalInclusion`](@ref) with a cylindrical-shape trait meets a general-anisotropic 3D stiffness."
struct CylinderQuadrature <: AbstractAlgorithm end

"Placeholder singleton used by `_resolve_algo(Val(:auto), …)` as an explicit *automatic* selection request."
struct Auto <: AbstractAlgorithm end

# ─── Material-symmetry traits ───────────────────────────────────────────────

"""
    MaterialSymmetry

Abstract supertype classifying the symmetry class of a TensND tensor.
Used by dispatch rules to select the most specific algorithm (e.g.
closed-form formulas for isotropic and transversely isotropic matrices).
"""
abstract type MaterialSymmetry end

"Isotropic material (TensND `TensISO`)."
struct IsotropicSym <: MaterialSymmetry end

"Transversely isotropic material (TensND `TensTI{4}` / `TensTI`)."
struct TransverselyIsotropicSym <: MaterialSymmetry end

"Orthotropic material (TensND `TensOrtho`)."
struct OrthotropicSym <: MaterialSymmetry end

"General anisotropic material (generic `Tens` subtype with no structured symmetry)."
struct GeneralAnisotropicSym <: MaterialSymmetry end

# ─── Classification of TensND tensors ───────────────────────────────────────

"""
    material_symmetry(C₀)

Return the [`MaterialSymmetry`](@ref) trait corresponding to the TensND
tensor `C₀`.  Dispatches on the concrete TensND type: `TensISO` →
`IsotropicSym`, `TensTI{4}` → `TransverselyIsotropicSym`, `TensOrtho`
→ `OrthotropicSym`, anything else → `GeneralAnisotropicSym`.
"""
function material_symmetry(::TensND.TensISO)
    return IsotropicSym()
end

# `TensTI{4}` and `TensOrtho` may or may not exist in any given TensND
# release — guard the method definitions with `isdefined` so the package
# still loads if those symbols have been renamed upstream.

if isdefined(TensND, :TensTI)
    @eval material_symmetry(::TensND.TensTI{4}) = TransverselyIsotropicSym()
end

if isdefined(TensND, :TensOrtho)
    @eval material_symmetry(::TensND.TensOrtho) = OrthotropicSym()
end

# Fallback — generic `Tens` (or any non-classified subtype)
material_symmetry(::TensND.AbstractTens) = GeneralAnisotropicSym()
