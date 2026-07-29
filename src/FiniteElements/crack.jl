# =============================================================================
#  crack.jl — the flat elliptical crack solved by 3-D finite elements.
#
#  Method: the "finite Eshelby cell with a first-order corrected boundary
#  condition" of Adessina, Barthélémy, Lavergne & Ben Fraj, *Int. J. Eng. Sci.*
#  119 (2017) 1-15, in the crack declination used in the `SifAniso` study —
#  the infinite matrix is truncated to a ball of radius `R = radius_ratio · a`
#  and the truncation bias is removed by adding the dipole far field of the
#  open crack to the imposed boundary displacement.  See
#  `docs/src/manual/fe_inclusions.md`.
# =============================================================================

"""
    FEMeshOptions(; radius_ratio = 5.0, htipdiv = 12.0, order = 2)

Discretization settings of a finite-element inclusion.

| Field | Default | Meaning |
|---|---|---|
| `radius_ratio` | `5.0` | radius of the surrounding ball of matrix, in units of the largest semi-axis. `5` is enough *because* the boundary condition is corrected — the uncorrected problem needs 10 to 40. |
| `htipdiv` | `12.0` | element size at the crack front, as `b / htipdiv`. The mesh coarsens to `R/3` at the outer boundary. |
| `order` | `2` | polynomial order of the displacement interpolation (1 or 2). Tetrahedra are geometrically straight; `order = 2` is therefore subparametric. |

`order = 2` is strongly recommended: linear tetrahedra badly under-resolve the
square-root field at the crack front.
"""
struct FEMeshOptions
    radius_ratio::Float64
    htipdiv::Float64
    order::Int
    function FEMeshOptions(; radius_ratio = 5.0, htipdiv = 12.0, order = 2)
        order in (1, 2) || throw(
            ArgumentError(
                "interpolation `order` must be 1 or 2 (Ferrite provides no " *
                    "higher-order Lagrange element on tetrahedra), got $order"
            )
        )
        radius_ratio > 1 ||
            throw(ArgumentError("`radius_ratio` must exceed 1, got $radius_ratio"))
        htipdiv > 0 || throw(ArgumentError("`htipdiv` must be positive, got $htipdiv"))
        return new(Float64(radius_ratio), Float64(htipdiv), Int(order))
    end
end

"""
    FEEllipticCrack(a, b; euler_angles = (), radius_ratio = 5.0,
                    htipdiv = 12.0, order = 2)

Flat elliptical crack whose **crack-opening-displacement tensor is computed by
finite elements** instead of the closed form of
[`EllipticCrack`](@ref MeanFieldHom.Cracks.EllipticCrack).

It subtypes [`AbstractCrack`](@ref) and declares the standard
[`shape_trait`](@ref MeanFieldHom.Core.shape_trait), so implementing
[`cod_tensor`](@ref MeanFieldHom.Cracks.cod_tensor) is *all* it takes: ℍ, ℕ,
𝐑, 𝐍_K, the bundled pair and the four `delta_*` with the Budiansky `4π/3`
prefactor are inherited. It is a drop-in replacement for `EllipticCrack` in
every scheme — the point of the exercise being that the same machinery accepts
a morphology for which no closed form exists.

Requires `Ferrite`, `FerriteGmsh` and `Gmsh` to be loaded.

# Scope

Isotropic matrix only. An anisotropic reference medium needs the anisotropic
Green-function gradient (Pan-Chou, or the Barnett-Willis line integral), which
is not implemented — see `docs/src/developer/roadmap.md`. Isotropy is tested on
the tensor's *content*, not its TensND type.

That restriction bites in the iterative schemes: `SelfConsistent` and
`DifferentialScheme` re-evaluate the crack in the current estimate, which for a
family of **parallel** cracks is transversely isotropic. Add
`symmetrize = IsoSymmetrize()` to the phase and the scheme hands the kernel an
isotropic reference at every iteration — which also collapses the whole
orientation family onto one cached solve.

`ForwardDiff` cannot be propagated through the solve (the sparse factorization
is `Float64`-only), so use finite differences for sensitivities.

# Example

```julia
using MeanFieldHom, Ferrite, FerriteGmsh, Gmsh

crack = FEEllipticCrack(1.0, 0.25; htipdiv = 12.0)
C₀ = iso_stiffness(0.8333, 0.3846)          # E = 1, ν = 0.3

B_fe = cod_tensor(crack, C₀)                 # finite elements
B_an = cod_tensor(EllipticCrack(1.0, 0.25), C₀)   # closed form

rve = RVE(:M)
add_matrix!(rve, Ellipsoid(1.0), Dict(:C => C₀))
add_phase!(rve, :cracks, crack, Dict(:C => C₀); density = 0.05,
           symmetrize = IsoSymmetrize())
homogenize(rve, MoriTanaka(), :C)
```

See also [`FEMeshOptions`](@ref), [`FECache`](@ref).
"""
struct FEEllipticCrack{T <: Number, S, B <: TensND.AbstractBasis} <: Core.AbstractCrack{T}
    a::T
    b::T
    basis::B
    mesh::FEMeshOptions
    cache::FECache
end

function FEEllipticCrack(
        a::Ta, b::Tb;
        euler_angles::Tuple{Vararg{Real}} = (),
        basis::Union{Nothing, TensND.AbstractBasis} = nothing,
        radius_ratio = 5.0, htipdiv = 12.0, order = 2
    ) where {Ta <: Real, Tb <: Real}
    T = Core._floatlike(promote_type(Ta, Tb))
    bas = basis === nothing ? Core._default_basis(T, euler_angles) : basis
    a_, b_ = T(a), T(b)
    a_ >= b_ > 0 || throw(
        ArgumentError("semi-axes must satisfy a ≥ b > 0, got a = $a_, b = $b_")
    )
    S = Cracks._classify_crack(T, a_, b_)
    opts = FEMeshOptions(; radius_ratio, htipdiv, order)
    return FEEllipticCrack{T, S, typeof(bas)}(a_, b_, bas, opts, FECache())
end

Core.shape_trait(::FEEllipticCrack{T, S}) where {T, S} = S
_fe_cache(c::FEEllipticCrack) = c.cache

function Core.shape_tensor(c::FEEllipticCrack{T}) where {T}
    D = zeros(T, 3, 3)
    D[1, 1] = c.a
    D[2, 2] = c.b
    return TensND.Tens(D, c.basis)
end

# ─── Extension seam ──────────────────────────────────────────────────────────

"""
    _fe_cod_tensor(crack, C₀; kw...)

Backend seam of the finite-element COD tensor. The implementation lives in the
package extension `MeanFieldHomFerriteExt`; this fallback is hit when the
extension is not loaded and raises an informative error.
"""
_fe_cod_tensor(args...; kwargs...) = error(
    "`FEEllipticCrack` needs the Ferrite extension: run " *
        "`import Ferrite, FerriteGmsh, Gmsh` first."
)

"""
    fe_cod_breakdown(crack, C₀) -> NamedTuple

Diagnostic view of the corrected finite-element solve — the COD tensor of the
truncated cell, the response to the crack's own dipole far field, and the
infinite-medium result. Implemented in `MeanFieldHomFerriteExt`; see its
docstring once the extension is loaded.
"""
fe_cod_breakdown(args...; kwargs...) = error(
    "`fe_cod_breakdown` needs the Ferrite extension: run " *
        "`import Ferrite, FerriteGmsh, Gmsh` first."
)

"""
    fe_mesh_report(crack) -> NamedTuple

Mesh diagnostics of a finite-element inclusion (cell/dof counts, welded
crack-front pairs, lip areas against the exact `πab`). Implemented in
`MeanFieldHomFerriteExt`.
"""
fe_mesh_report(args...; kwargs...) = error(
    "`fe_mesh_report` needs the Ferrite extension: run " *
        "`import Ferrite, FerriteGmsh, Gmsh` first."
)

Cracks.cod_tensor(
    crack::FEEllipticCrack, C₀::TensND.AbstractTens{4, 3};
    K_interface::Union{Nothing, TensND.AbstractTens{2, 3}} = nothing, kw...
) = let B = _fe_cod_tensor(crack, C₀; kw...)
    K_interface === nothing ? B :
        Cracks._apply_interface_stiffness(B, K_interface, crack.b)
end

Cracks.cod_tensor(::FEEllipticCrack, ::TensND.AbstractTens{2, 3}; kw...) = error(
    "`FEEllipticCrack` implements the elastic COD tensor only. The transport " *
        "problem would need its own finite-element resolution (see " *
        "`docs/src/developer/roadmap.md`); use `EllipticCrack` for conduction."
)
