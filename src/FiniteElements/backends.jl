# =============================================================================
#  backends.jl — which finite-element library performs the solve.
#
#  The scientific content of the axisymmetric solver — the Fourier operators,
#  the boundary data, the fixed point of the corrected boundary condition, the
#  memoization — lives in `src/`.  What a backend supplies is only the
#  discretization: a mesh, scalar Lagrange spaces, an assembly, a Dirichlet
#  lift and a quadrature.  That is the nine functions declared below.
#
#  A backend never sees a Fourier mode or a physics: the driver closes `Bop`
#  and `proj` over the mode before handing them over, so the backend is told
#  only "this many scalar fields, this operator, this projection".
# =============================================================================

"""
    FEBackend

Which finite-element library performs the solve. Concrete singletons:
[`FerriteBackend`](@ref), [`GridapBackend`](@ref) and the default
[`AutoBackend`](@ref).

The three are always defined, whether or not the corresponding package is
loaded — naming a backend costs nothing, only *solving* with it requires the
matching extension.
"""
abstract type FEBackend end

"""
    AutoBackend()

Default backend of a finite-element inclusion: pick whichever backend is
loaded, **at the first solve** rather than at construction.

Deferring the choice is what lets an inclusion be built, printed, stored in an
`RVE` and passed around in a session where no finite-element package has been
imported; the informative error arrives only when a scheme actually asks for a
localization tensor.

Priority is `FerriteBackend` then `GridapBackend`. Loading both is not an
error — it just means the first is chosen. To pick deliberately, pass
`backend = GridapBackend()` to the constructor.
"""
struct AutoBackend <: FEBackend end

"""
    FerriteBackend()

Solve with [Ferrite.jl](https://ferrite-fem.github.io/Ferrite.jl); needs
`import Ferrite, FerriteGmsh, Gmsh`. The reference implementation, and the
only backend for [`FEEllipticCrack`](@ref).
"""
struct FerriteBackend <: FEBackend end

"""
    GridapBackend()

Solve with [Gridap.jl](https://gridap.github.io/Gridap.jl); needs
`import Gridap, GridapGmsh` (GridapGmsh carries its own `gmsh`, so `Gmsh.jl`
is not required on this path).

Available for [`FEExcenteredSphere`](@ref) only. Gridap states the weak form
directly — `∫( Bᵐ(v)' * D * Bᵐ(u) * ρ )dΩ` — which makes it the easier of the
two to read and to modify; Ferrite's explicit assembly loop is the faster of
the two to run.
"""
struct GridapBackend <: FEBackend end

_backend_extension(::FerriteBackend) = :MeanFieldHomFerriteExt
_backend_extension(::GridapBackend) = :MeanFieldHomGridapExt

_backend_import(::FerriteBackend) = "import Ferrite, FerriteGmsh, Gmsh"
_backend_import(::GridapBackend) = "import Gridap, GridapGmsh"

_backend_loaded(b::FEBackend) =
    Base.get_extension(parentmodule(@__MODULE__), _backend_extension(b)) !== nothing

const _AXI_BACKENDS = (FerriteBackend(), GridapBackend())

"""
    _resolve_backend(b) -> FEBackend

Turn [`AutoBackend`](@ref) into a concrete backend, and check that a concrete
one is actually available. Called once per inclusion, at the first solve; the
result is pinned in the cache, so a single inclusion never mixes two backends'
grids.
"""
function _resolve_backend(b::FEBackend)
    _backend_loaded(b) || error(
        "the $(nameof(typeof(b))) finite-element backend is not loaded: " *
            "run `$(_backend_import(b))` first."
    )
    return b
end

function _resolve_backend(::AutoBackend)
    for b in _AXI_BACKENDS
        _backend_loaded(b) && return b
    end
    return error(
        "no finite-element backend is loaded: run `" *
            join(map(_backend_import, _AXI_BACKENDS), "` or `") * "` first."
    )
end

# ─── The backend contract ────────────────────────────────────────────────────
#
#  Nine functions, each with a fallback that names the offending backend.  A
#  new backend is exactly this list; nothing else in the package needs to know
#  that it exists.

_no_backend_method(f, b) = error(
    "`$f` is not implemented for $(nameof(typeof(b))). Either the extension " *
        "failed to load (`$(_backend_import(b))`), or this backend does not " *
        "support this inclusion."
)

"""
    fe_axi_grid(backend, incl)

Backend-native mesh of the meridian half-plane of `incl`, carrying the cell
sets `"core"`, `"shell"`, `"matrix"` and the boundary sets `"outer"`, `"axis"`
of [`_build_gmsh_axi_model`](@ref). The first node coordinate is the
cylindrical radius `ρ`, the second the axial coordinate `z`.
"""
fe_axi_grid(b::FEBackend, incl) = _no_backend_method("fe_axi_grid", b)

"""
    fe_axi_grid_counts(backend, grid) -> (; ncells, nnodes, ncells_by_set)

Cell and node counts, `ncells_by_set` being a `Dict{String,Int}` over the three
regions. Diagnostics only.
"""
fe_axi_grid_counts(b::FEBackend, grid) = _no_backend_method("fe_axi_grid_counts", b)

"""
    fe_axi_region_volume(backend, grid, set) -> Float64

Volume of revolution `2π ∫_set ρ dρ dz` of one region, on the geometric
interpolation of the mesh. Diagnostics only — the driver measures its own
volume with the mode's quadrature, and the two need not agree exactly.
"""
fe_axi_region_volume(b::FEBackend, grid, set) =
    _no_backend_method("fe_axi_region_volume", b)

"""
    fe_axi_mode(backend, grid, order, ncomp, axis_zeros) -> mode

Discretize one Fourier mode: `ncomp` scalar Lagrange fields of degree `order`
on `grid`, with a quadrature exact to degree `2 * order + 1`.

The dof numbering must span the **whole** space, with no Dirichlet elimination
— the driver does the free/prescribed split itself, so that one factorization
serves every right-hand side.

`axis_zeros` lists the component indices that the axis regularity of this mode
forces to vanish on the set `"axis"`.
"""
fe_axi_mode(b::FEBackend, grid, order, ncomp, axis_zeros) =
    _no_backend_method("fe_axi_mode", b)

"""
    fe_axi_dof_split(backend, mode) -> (ndofs, free, presc)

Total dof count and the two index vectors, in the numbering of
[`fe_axi_mode`](@ref). `presc` is the sorted union of the outer-boundary dofs
and the axis-pinned dofs; `free` is its complement.
"""
fe_axi_dof_split(b::FEBackend, mode) = _no_backend_method("fe_axi_dof_split", b)

"""
    fe_axi_set_dirichlet!(backend, mode, u, f) -> u

Write the Dirichlet data of one right-hand side into `u`: `u[d] = f(ρ_d, z_d)[k]`
for every dof `d` of the outer boundary, `k` being its component, **then**
`u[d] = 0` for every axis-pinned dof.

Two things make this the delicate function of the contract.

*The order matters.* The poles `(0, ±R)` belong to both the `"outer"` and the
`"axis"` sets. The axis must win, so the zeroing comes second.

*It must not touch the matrix.* The driver assembles and factorizes once, then
calls this once per load case. Re-assembling here would multiply the cost of a
solve by the number of loads.
"""
fe_axi_set_dirichlet!(b::FEBackend, mode, u, f) =
    _no_backend_method("fe_axi_set_dirichlet!", b)

"""
    fe_axi_stiffness(backend, mode, Dmap, Bop) -> AbstractMatrix

Stiffness of one mode over the whole dof numbering,

```
K = Σ_regions ∫_region Bop(v)' * D_region * Bop(u) * ρ dρ dz .
```

`Dmap` is a `Vector{Pair{String,Matrix{Float64}}}` — a vector, not a `Dict`, so
that the assembly order is reproducible — mapping a cell-set name to that
region's material matrix in the cylindrical `(ρ, θ, z)` basis.

`Bop(N, dNρ, dNz, ρ) -> Matrix{Float64}` of size `nrow × ncomp` is the
generalized-strain operator of one scalar shape function, already closed over
the Fourier mode. It is **R-linear in `(N, dNρ, dNz)`**, so a backend that
manipulates whole trial functions rather than shape functions may apply it to
`(u_c, ∂ρu_c, ∂zu_c)` directly instead of building an element `B` matrix.

The `ρ` in the measure is the single factor that turns a plane problem into a
solid of revolution.
"""
fe_axi_stiffness(b::FEBackend, mode, Dmap, Bop) =
    _no_backend_method("fe_axi_stiffness", b)

"""
    fe_axi_average(backend, mode, Dmap, u, Bop, proj, sets) -> (prim, dual, V)

Volume averages over `∪ sets` of the generalized strain and of the associated
generalized stress, both projected by `proj` onto the Kelvin basis of the mode,
plus the volume of revolution `V = 2π ∫ ρ dρ dz` of that union:

```
prim = ∫ proj(B u) ρ / ∫ ρ ,      dual = ∫ proj(D · B u) ρ / ∫ ρ .
```

The azimuthal integration has already been performed analytically inside
`proj`; what remains is the meridian quadrature.
"""
fe_axi_average(b::FEBackend, mode, Dmap, u, Bop, proj, sets) =
    _no_backend_method("fe_axi_average", b)
