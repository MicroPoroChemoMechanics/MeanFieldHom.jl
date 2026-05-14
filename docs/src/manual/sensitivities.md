# Sensitivities — autodiff via ForwardDiff

Since v0.4.0 `MeanFieldHom` exposes a small but expressive API for computing
derivatives of `homogenize(rve, scheme)` outputs with respect to any scalar
input parameter — physical (stiffness coefficients, conductivities) or
geometric (radii, aspect ratios, volume fractions, crack densities, distribution
shape envelopes) — *and* with respect to scalar fields of inclusion types
defined later by the user.

The whole machinery is a thin convenience layer on top of [ForwardDiff.jl];
ForwardDiff is shipped as a [weak dependency](https://pkgdocs.julialang.org/v1/creating-packages/#Conditional-loading-of-code-in-packages-(Extensions))
so the API only activates when you `using ForwardDiff` alongside `MeanFieldHom`.

## Why autodiff

Every kernel of the package (`hill_tensor`, `eshelby_tensor`, the ten
schemes, the SC/ASC/Differential solvers) is engineered to be
`ForwardDiff.Dual`-friendly: cross-checked by `test/Schemes/test_dual_compat.jl`
and validated via direct propagation through every numerical branch. The
sensitivity API just lifts that property into a friendly user-facing form,
so you do not have to rebuild the RVE manually with `Dual` values, plumb
the right `Tag`, or extract scalars by hand.

Practical consequences:

- **Geometry parameters are first-class.** Differentiate w.r.t. semi-axes,
  crack opening, distribution-shape envelopes, anything stored as a `Number`
  field on an inclusion.
- **User-defined inclusions just work.** Define your own
  `<: AbstractEllipsoidalInclusion` (or any subtype of `AbstractInclusion`)
  with `Number` fields, register `hill_tensor` and friends — and the
  sensitivity API differentiates them with no further code change.
- **Multi-scale chain rule is automatic.** Compose two or three calls to
  `homogenize` in a closure; ForwardDiff propagates partial derivatives
  through all of them in a single Dual sweep — no manual chain rule
  required.

## Quick start

```julia
using MeanFieldHom, ForwardDiff, TensND

rve = RVE(:M)
add_matrix!(rve, Ellipsoid(1.0), Dict(:C => TensISO{3}(30.0, 10.0)))
add_phase!(rve, :I, Ellipsoid(1.0),
            Dict(:C => TensISO{3}(60.0, 20.0)); fraction = 0.2)

# Sensitivity of C_eff[1,1,1,1] w.r.t. the inclusion volume fraction
∂_f = derivative(rve, MoriTanaka(), amount(:I);
                 indexer = C -> get_array(C)[1, 1, 1, 1])
```

That's it — no closure to write, no manual `set_param`/`homogenize`
plumbing.

## Parameter lenses

The first argument to `derivative`/`gradient`/`jacobian` is an
`AbstractParameter` *lens* describing the scalar input you want to
differentiate against. Four concrete lens kinds are shipped, plus a
`sensitivity(f, x₀)` closure fallback for everything else.

| Helper                          | Kind                                | Underlying type             |
| ------------------------------- | ----------------------------------- | --------------------------- |
| `amount(:I)`                    | volume fraction or crack density    | `AmountParameter`           |
| `property(:I, :C, :bulk)`       | scalar coefficient of a tensor      | `PropertyParameter`         |
| `geometry(:I, :semi_axes, 3)`   | scalar geometry field               | `GeometryParameter`         |
| `shape_param(:semi_axes, 1)`    | distribution-shape geometry field   | `DistributionShapeParameter`|

Named selectors recognised by `property` (other symbols fall back to a
positional `Int` index into `get_data(tensor)`):

| Tensor type    | Named selectors                                      |
| -------------- | ---------------------------------------------------- |
| `TensISO{2}`   | `:scalar`, `:λ`                                      |
| `TensISO{4,3}` | `:bulk`, `:K`, `:α` ; `:shear`, `:μ`, `:β`           |
| `TensTI{2}`    | `:transverse`, `:a` ; `:axial`, `:b`                 |
| `TensTI{4}`    | `:ℓ₁` … `:ℓ₆` (with `ℓ₃ = ℓ₄` in the major-symmetric case) |

## Single derivative, gradient, full Jacobian

```julia
# scalar in / scalar out
∂_K = derivative(rve, MoriTanaka(), property(:I, :C, :bulk);
                 indexer = C -> get_array(C)[1, 1, 1, 1])

# vector in / scalar out
ps = [amount(:I), property(:I, :C, :bulk), property(:M, :C, :shear)]
∇  = gradient(rve, MoriTanaka(), ps;
              indexer = C -> get_array(C)[1, 1, 1, 1])

# vector in / tensor out → full Jacobian (flattened to 81 × N for a 4-tensor)
J = jacobian(rve, MoriTanaka(), ps)         # size(J) == (81, 3)
```

`gradient` and `jacobian` accept an optional `chunk = ForwardDiff.Chunk(N)`
kwarg; without it ForwardDiff picks a chunk size automatically.

## Closure fallback for arbitrary parameterisations

Anything that can't be expressed as a single lens (composite parameters,
parameters of a user inclusion that don't map to a `Number` field, etc.)
can be differentiated through a user-supplied closure:

```julia
∂α = sensitivity(0.3) do α
    rve = RVE(:M)
    add_matrix!(rve, Ellipsoid(1.0), Dict(:C => TensISO{3}(30.0, 10.0)))
    add_phase!(rve, :I, Ellipsoid(1.0, α, α^2),
                Dict(:C => TensISO{3}(60.0, 20.0)); fraction = 0.2)
    return get_array(homogenize(rve, MoriTanaka(), :C))[1, 1, 1, 1]
end
```

`sensitivity(f, x₀)` auto-detects the derivative / gradient / Jacobian flavor
from the shape of `x₀` and the return type of `f(x₀)`. Pass
`kind = :derivative | :gradient | :jacobian` to force a specific mode.

## User-defined inclusions

```julia
struct MyBlob{T <: Number, B <: TensND.AbstractBasis} <:
       MeanFieldHom.AbstractEllipsoidalInclusion{3, T}
    radius::T
    eccentricity::T
    basis::B
end
# Register hill_tensor / eshelby_tensor (delegate to an equivalent Ellipsoid)
MeanFieldHom.hill_tensor(b::MyBlob, C₀::TensND.AbstractTens; kw...) =
    MeanFieldHom.hill_tensor(Ellipsoid(b.radius, b.radius*(1-b.eccentricity),
                                        b.radius*(1-b.eccentricity)^2), C₀; kw...)

# Differentiate w.r.t. any scalar field — no library change required.
∂_e = derivative(rve, Dilute(), geometry(:B, :eccentricity);
                 indexer = C -> get_array(C)[1,1,1,1])
```

The generic `_replace_geom_field` reflects on the struct's fieldnames and
promotes any sibling `<:Number` field to the new `Dual` element type, so
the parametric inner constructor of `MyBlob{T,B}` resolves cleanly.

## Multi-scale chain rule

Composing several `homogenize` calls in a single closure means the chain
rule is taken care of automatically. The script
[`scripts/28_multiscale_strength.jl`](https://codeberg.org/MicroPoroChemoMechanics/MeanFieldHom.jl/src/branch/main/scripts/28_multiscale_strength.jl)
is a complete worked example following Pichler & Hellmich (CCR 2011): a
three-scale upscaling of cement-paste / mortar elasticity and quasi-brittle
strength built from a self-consistent hydrate foam plus two Mori-Tanaka
stages. The strength criterion requires `∂C_mortar / ∂μ_hyd` through the
full chain — exposed in two complementary styles:

- **Manual chain rule** — one `jacobian` call per scale, then an explicit
  tensor product of partial Jacobians. Useful when the intermediate
  partials are themselves of interest.
- **End-to-end autodiff** — a single `ForwardDiff.derivative` (or
  `sensitivity`) on a closure that runs the three nested schemes in one
  pass. Shorter, scheme-agnostic, and naturally extends to longer chains.

Both approaches agree to the floor of ForwardDiff itself, and the
end-to-end approach is the recommended default unless you specifically
need the intermediate Jacobians.

## Symmetrize and orientation distributions

Inclusions are often modelled as oriented distributions rather than a
single oriented inclusion: a thin oblate spheroid with a *uniform spatial
distribution of orientations* is, on average, isotropic. The
`symmetrize` keyword on `add_matrix!` / `add_phase!` declares such a
distribution at the RVE level so the homogenization kernel projects the
phase localisation tensor onto the corresponding symmetry class:

| Symmetrize value           | Meaning                                                     | Result class |
| -------------------------- | ----------------------------------------------------------- | ------------ |
| `:none` (default)          | inclusion at its declared orientation                       | as input     |
| `:iso`                     | uniform spatial distribution (Reynolds avg over `SO(3)`)    | `TensISO`    |
| `:ti` / `TISymmetrize(n)`  | uniform azimuthal distribution around axis `n` (default ez) | `TensTI(n)`  |

```julia
rve = RVE(:M)
add_matrix!(rve, Ellipsoid(1.0), Dict(:C => TensISO{3}(216.0, 64.0)))
# Oblate inclusions with a uniform-in-orientation distribution: the
# effective stiffness is iso even though each individual inclusion is TI.
add_phase!(rve, :I, Ellipsoid(1.0, 1.0, 0.2),
            Dict(:C => TensISO{3}(3.0e-6, 2.0e-6));
            fraction = 0.3, symmetrize = :iso)
```

Sensitivities on RVEs that carry a `symmetrize` keyword work the same way:
`derivative` / `gradient` / `jacobian` propagate Duals through the
projection automatically.

## Limitations

- **No `Complex{T}` autodiff.** ForwardDiff doesn't mix Dual + Complex
  cleanly. Use closure-style `sensitivity` only when the input is real.
  Frequency-domain viscoelastic computations remain available via
  `Complex{Float64}` moduli (independent of the autodiff API).
- **Symbolic differentiation is not exposed via this API.** Use SymPy
  directly through the closed-form schemes (Voigt, Reuss, Dilute, MT,
  Maxwell, PCW) when symbolic derivatives are needed.
- **Geometry derivative across shape categories.** Perturbing axes of an
  `Ellipsoid{Spherical}` inclusion preserves the `Spherical` shape trait
  in the parametric type, so the symmetry-imposed derivative is `0`. To
  get a non-trivial geometric sensitivity, perturb around an already
  non-degenerate (`Triaxial` / `Prolate` / `Oblate`) configuration.
- **TI symmetrize on non-coaxial inclusions** routes the matrix used for
  the localisation-tensor computation through an isotropic projection —
  the result still satisfies the outer `TI(axis)` projection, and is
  exact at the iso fixed-point of the SC iteration. This is documented
  in [`src/Schemes/symmetrize.jl`](https://codeberg.org/MicroPoroChemoMechanics/MeanFieldHom.jl/src/branch/main/src/Schemes/symmetrize.jl).

## Validation

The package ships a cross-cutting test in
[`test/Schemes/test_sensitivities.jl`](https://codeberg.org/MicroPoroChemoMechanics/MeanFieldHom.jl/src/branch/main/test/Schemes/test_sensitivities.jl)
that compares every sensitivity against centred finite differences on every
scheme, plus an exact match against the Christensen 1990 closed form for
`∂k_MT/∂f`. Expect agreement to `rtol ≈ 1e-6` for closed-form schemes and
`rtol ≈ 1e-4` for iterative schemes (limited by the fixed-point tolerance,
not the autodiff itself).
