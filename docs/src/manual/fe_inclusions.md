# [Finite-element inclusions](@id man-fe-inclusions)

When a morphology has no closed-form Eshelby solution, the response can be
computed numerically and fed to the schemes through the
[custom-inclusion contract](@ref man-custom-inclusions). `MeanFieldHom` ships
one such inclusion: [`FEEllipticCrack`](@ref MeanFieldHom.FEEllipticCrack),
whose crack-opening-displacement tensor comes out of a `Ferrite.jl` solve.

```julia
using MeanFieldHom
import Ferrite, FerriteGmsh, Gmsh      # activates MeanFieldHomFerriteExt
```

Without those three packages the type still exists but every solve raises an
informative error.

## The method

A finite-element Eshelby problem has to truncate the infinite matrix. Doing so
naively — imposing the remote field `u = E·x` on the boundary of a ball of
radius `R` — leaves an ``O((a/R)^3)`` bias, which is why an uncorrected
computation needs `R/a` of 10 to 40.

The fix, from
[Adessina et al. 2017](https://doi.org/10.1016/j.ijengsci.2017.03.015), is to
put the **dipole far field of the inclusion itself** into the boundary
condition. For a flat crack this takes the "3 + 3" form used here:

```math
\mathbf u(\mathbf x)\;\underset{\|\mathbf x\|\to\infty}{\approx}\;
  (\mathbb S_0:\boldsymbol\Sigma)\cdot\mathbf x
  \;-\; b\,S_f\bigl(\nabla\mathbf G(\mathbf x):\mathbb C_0\cdot\hat{\mathbf n}\bigr)\cdot\mathbf U,
\qquad \mathbf U=\frac{\langle[\![\mathbf u]\!]\rangle}{b},\quad S_f=\pi ab .
```

Per evaluation:

1. one assembly and **one** Cholesky factorization, shared by all six solves;
2. three *traction* solves, `u|∂Ω = (𝕊₀:Σ⁽ⁱ⁾)·x` with `Σ⁽ⁱ⁾·n̂ = eᵢ`, giving the
   COD tensor `B_s` of the truncated cell;
3. three *dipole* solves, `u|∂Ω = −∇G(x):Πₘ` with `Πₘ = b·S_f·ℂ₀:(eₘ⊗ˢn̂)`,
   giving its response `B_u` to the crack's own far field;
4. the fixed point, in closed form: `B_∞ = (𝟏 − B_u)⁻¹·B_s`.

`∇G` is the real-space Kelvin Green gradient,
[`green_gradient_iso`](@ref MeanFieldHom.Core.green_gradient_iso) /
[`dipole_displacement_iso`](@ref MeanFieldHom.Core.dipole_displacement_iso),
which is why the matrix must currently be isotropic.

The crack itself is a **zero-thickness discontinuity** — duplicated nodes —
whose lips are traction-free naturally: no interface term, no multiplier, no
contact condition. The mean opening is measured as a surface integral of the
jump over each lip, with no assumption on the opening profile.

## Using it

```julia
crack = FEEllipticCrack(1.0, 0.25; htipdiv = 12.0)

C₀ = iso_stiffness(0.8333, 0.3846)                 # E = 1, ν = 0.3
B_fe = cod_tensor(crack, C₀)
B_an = cod_tensor(EllipticCrack(1.0, 0.25), C₀)    # closed form, for comparison

rve = RVE(:M)
add_matrix!(rve, Ellipsoid(1.0), Dict(:C => C₀))
add_phase!(rve, :cracks, crack, Dict(:C => C₀); density = 0.05,
           symmetrize = IsoSymmetrize())
homogenize(rve, MoriTanaka(), :C)
```

That is the whole point of the contract: only `cod_tensor` is implemented, and
because the type subtypes `AbstractCrack` with a standard
[`shape_trait`](@ref MeanFieldHom.Core.shape_trait), ℍ, ℕ, 𝐑, 𝐍_K, the bundled
pair and the four `delta_*` with their Budiansky ``4\pi/3`` prefactor all
follow. Orientation averaging is applied by the scheme afterwards.

### Discretization

[`FEMeshOptions`](@ref MeanFieldHom.FEMeshOptions) — passed as keywords to the
constructor:

| Keyword | Default | Meaning |
|---|---|---|
| `radius_ratio` | `5.0` | `R/a`. Five is enough *because* the boundary condition is corrected. |
| `htipdiv` | `12.0` | element size at the crack front, `b/htipdiv`; the mesh coarsens to `R/3` at the outer boundary. |
| `order` | `2` | displacement interpolation order. |

Straight tetrahedra with a P2 displacement field. Curving the geometry is not
an option: `setOrder(2)` runs *after* gmsh's `Crack` plugin and curves the two
lips' front edges differently, which would pull the welded crack front apart
again at the mid-side nodes.

### Diagnostics

```julia
fe_mesh_report(crack)     # cells, dofs, welded front pairs, lip areas vs πab
fe_cod_breakdown(crack, C₀)   # B_s, B_u, B_inf — bypasses the cache
fe_assembly_count(crack)  # factorizations actually performed
fe_reset!(crack)          # drop the mesh and the memoized tensors
```

`norm(B_u)` is the size of the truncation bias being removed; it falls like
``(a/R)^3``, and checking that ratio is the sharpest test that the correction
is wired correctly.

### Cost and memoization

One evaluation is one assembly, one factorization and six solves. Results are
**memoized on the reference medium expressed in the crack's own frame**, which
has a useful consequence: a whole family of orientations of the same crack in
the same isotropic matrix shares a single finite-element resolution.

One-shot schemes (`Dilute`, `MoriTanaka`, `Maxwell`, `PonteCastanedaWillis`)
therefore cost exactly one solve. Iterative schemes (`SelfConsistent`,
`DifferentialScheme`) change the reference medium at every step and pay one
solve per distinct `C₀` — usable, but expensive.

!!! note "Iterative schemes need `IsoSymmetrize`"
    An iterative scheme re-evaluates the crack in the *current estimate*. For a
    family of **parallel** cracks that estimate is transversely isotropic, and
    the isotropic-only boundary correction refuses it — with an error saying
    so, rather than silently substituting an isotropic projection. Adding
    `symmetrize = IsoSymmetrize()` to the phase makes the scheme hand the
    kernel an isotropic reference at every iteration, and has the pleasant side
    effect of turning the whole orientation family into a single cached solve.

    Isotropy is tested on the tensor's *content*, not its TensND type, so an
    iterate arriving as a `TensCanonical` with isotropic content is accepted.

## Accuracy

The crack front carries a square-root displacement field, so a fixed-order
element converges slowly — in practice **first order in the element size**.
On the meshes above the finite-element opening sits a few percent *below* the
closed form (a truncated cell is stiffer, and the front is under-resolved), and
Richardson extrapolation in `h` lands within about a percent. Script
`81_fe_crack_eshelby.jl` runs that convergence study.

For reference, the FEniCSx implementation of the same scheme in the `SifAniso`
study reports ±5 % on ``\mathbf B_\infty`` at `htipdiv = 12` with P3 elements.

## Limitations

- **Isotropic matrix only.** The corrected boundary condition uses the
  closed-form Kelvin dipole; the anisotropic Green gradient (Pan-Chou, or the
  Barnett-Willis line integral) is not implemented. An anisotropic `C₀` raises
  an `ArgumentError`.
- **Elasticity only.** The transport problem would need its own finite-element
  resolution.
- **No automatic differentiation** through the solve — the sparse factorization
  is `Float64`-only. Use finite differences for sensitivities.
- Ferrite provides no cubic Lagrange element on tetrahedra, so `order` is
  capped at 2.

These, along with the general 6 + 6 scheme for solid inclusions (the
excentered-core sphere of the 2017 paper) and Fourier axisymmetric elements,
are tracked in `docs/src/developer/roadmap.md`.

## See also

- [Custom inclusions](@ref man-custom-inclusions) — the contract this
  implements.
- [Adding a new inclusion](@ref) — the full developer contract.
- `scripts/81_fe_crack_eshelby.jl`, `scripts/82_fe_crack_schemes.jl`.
