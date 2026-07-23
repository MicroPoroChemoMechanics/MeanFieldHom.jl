# Localization and contribution tensors

MeanFieldHom exposes the four **dilute localization tensors** of the
Eshelby problem, together with the size-independent **stiffness and
compliance contribution tensors** of Kachanov–Sevostianov.

## Pivot formula

For an inclusion `I` of stiffness ``\mathbb C_1`` embedded in an
infinite matrix of stiffness ``\mathbb C_0``, the strain–strain
localization tensor is

```math
\mathbb A_{\varepsilon\varepsilon} = \bigl[\,\mathbb I +
\mathbb P(\mathrm I, \mathbb C_0) :
(\mathbb C_1 - \mathbb C_0)\,\bigr]^{-1},
```

where ``\mathbb P`` is the Hill polarization tensor
([`hill_tensor`](@ref)) and ``\mathbb I`` is the symmetric identity
4-tensor.  The three other localization tensors follow algebraically:

```math
\mathbb A_{\sigma\varepsilon} = \mathbb C_1 : \mathbb A_{\varepsilon\varepsilon},\qquad
\mathbb A_{\varepsilon\sigma} = \mathbb A_{\varepsilon\varepsilon} : \mathbb S_0,\qquad
\mathbb A_{\sigma\sigma} = \mathbb C_1 : \mathbb A_{\varepsilon\varepsilon} : \mathbb S_0,
```

with ``\mathbb S_0 = \mathbb C_0^{-1}``.  The four functions exposed by
MeanFieldHom are:

| Function                                     | Return value                               |
| -------------------------------------------- | ------------------------------------------ |
| [`strain_strain_loc`](@ref)`(incl, C₁, C₀)`  | ``\mathbb A_{\varepsilon\varepsilon}``     |
| [`stress_strain_loc`](@ref)`(incl, C₁, C₀)`  | ``\mathbb A_{\sigma\varepsilon}``          |
| [`strain_stress_loc`](@ref)`(incl, C₁, C₀)`  | ``\mathbb A_{\varepsilon\sigma}``          |
| [`stress_stress_loc`](@ref)`(incl, C₁, C₀)`  | ``\mathbb A_{\sigma\sigma}``               |

## Contribution tensors

The **stiffness contribution tensor** (Kachanov–Sevostianov 2018) is

```math
\mathbb N = (\mathbb C_1 - \mathbb C_0) : \mathbb A_{\varepsilon\varepsilon},
```

and its dilute-scheme volume average is

```math
\Delta\mathbb C_\mathrm{eff} = f \,\mathbb N,
```

for a dilute family of volume fraction ``f``.  The dual **compliance
contribution tensor** is

```math
\mathbb H = (\mathbb S_1 - \mathbb S_0) : \mathbb A_{\sigma\sigma},\qquad
\Delta\mathbb S_\mathrm{eff} = f\,\mathbb H.
```

Functions: [`stiffness_contribution`](@ref),
[`compliance_contribution`](@ref), with density helpers
[`delta_stiffness`](@ref) and [`delta_compliance`](@ref).

## Cracks (Kachanov convention)

For flat cracks the Budiansky density convention is used instead of a
volume fraction.  The same entry points apply, with the density
``\varepsilon`` replacing ``f``:

- `compliance_contribution(crack, C₀)` returns the size-independent
  `H = (3/4) n̂ ⊗ˢ B ⊗ˢ n̂` (elliptic) or `(2/π) n̂ ⊗ˢ B ⊗ˢ n̂` (ribbon);
- `stiffness_contribution(crack, C₀)` returns `N = -C₀ : H : C₀`
  (first order, provided for API symmetry);
- `delta_compliance(crack, H, ε)` and `delta_stiffness(crack, N, ε)`
  apply the appropriate `(4π/3)` or `π` geometric prefactor.

## Conductivity (2nd-order transport)

Every routine above has a 2-tensor analog, triggered by dispatch on
`::AbstractTens{2,3}` matrices:

| Elasticity                           | Conductivity                           |
| ------------------------------------ | -------------------------------------- |
| `strain_strain_loc`                  | [`gradient_gradient_loc`](@ref)        |
| `stress_strain_loc`                  | [`flux_gradient_loc`](@ref)            |
| `strain_stress_loc`                  | [`gradient_flux_loc`](@ref)            |
| `stress_stress_loc`                  | [`flux_flux_loc`](@ref)                |
| `stiffness_contribution`             | [`conductivity_contribution`](@ref)    |
| `compliance_contribution` (ellipsoid)| [`resistivity_contribution`](@ref)     |
| `delta_stiffness`                    | [`delta_conductivity`](@ref)           |
| `delta_compliance` (ellipsoid)       | [`delta_resistivity`](@ref)            |

## Type-genericity

All four localization and both contribution tensors are generic in the
element type — Float64, BigFloat, `ForwardDiff.Dual`, `SymPy.Sym`, and
`Symbolics.Num` all flow through the computation.  The only
requirement is that [`hill_tensor`](@ref) supports the chosen element
type.

## Extending to user-defined inclusions

Any concrete subtype of `AbstractInclusion` inherits the four
localization and the contribution tensors automatically once it
provides a method for [`hill_tensor`](@ref).  For inclusions where the
Hill polarization tensor has no convenient closed form (e.g.
`LayeredSphere`), the user may directly override
[`strain_strain_loc`](@ref) (and, if needed, its variants) — the three
remaining localization tensors and the contribution tensors are
derived algebraically and do not require additional methods.

See the developer guide [Adding a new inclusion](../developer/adding_inclusion.md)
for a step-by-step recipe.
