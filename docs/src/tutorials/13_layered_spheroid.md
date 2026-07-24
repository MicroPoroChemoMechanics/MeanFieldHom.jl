# The n-layer confocal spheroid (conduction)

Every RVE phase seen so far has been a single ellipsoid or a layered
*sphere*. This page introduces
[`LayeredSpheroid`](@ref) — an `N`-layer confocal spheroidal composite
inclusion for **conduction** (thermal / electric / Darcy) — following
[Barthélémy & Bignonnet (2020)](@cite barthelemyBignonnetIJES2020). It
plugs into the same `RVE`/[`homogenize`](@ref) workflow as every other
phase, so if you have read
[the first estimate](01_first_estimate.md) tutorial, most of this page
will look familiar. For the underlying confocal-harmonic series and the
numerical-precision choices behind it, see
[the theory page](../theory/layered_spheroid.md); for four worked
application scripts (interface-parameter sweeps, series convergence,
equivalent-particle conductivity, local fields), see the
[Gallery](../gallery/index.md).

!!! note "Conduction only"
    Unlike [`LayeredSphere`](@ref), `LayeredSpheroid` has **no elastic
    counterpart** — the harmonic decomposition it relies on is specific
    to the scalar Laplace equation.

## Building a single-layer spheroid

A spheroid is specified by its per-layer **axis** (revolution) and
**disk** (transverse) semi-axes — ascending, and confocal (constant
`axis² - disk²`, up to sign):

```@example spheroid
using MeanFieldHom
using TensND
using LinearAlgebra

K1 = TensISO{3}(20.0)   # inclusion conductivity
K0 = TensISO{3}(2.0)    # matrix conductivity

s = LayeredSpheroid((3.0,), (1.0,), (K1,))   # a prolate spheroid, ω = 3
```

With a single layer and a perfect interface, `LayeredSpheroid` must
reduce exactly to the classical Eshelby result for a homogeneous
spheroidal inclusion — a useful sanity check, and the first regression
test in `test/LayeredSpheroids/test_conductivity.jl`:

```@example spheroid
A_layered = gradient_gradient_loc(s, K1, K0)

ell = Spheroid(3.0)                 # equivalent homogeneous Ellipsoid
P = hill_tensor(ell, K0)
A_classic = inv(TensISO{3}(1.0) + P ⋅ (K1 - K0))

get_array(A_layered) ≈ get_array(A_classic)
```

`Spheroid(ω)`'s revolution axis is `ê₁` for prolate (`ω > 1`) and `ê₃`
for oblate (`ω < 1`); `LayeredSpheroid`'s own default axis is `ê₃`,
overridable with the `axis` keyword — pass `axis = (1.0, 0.0, 0.0)`
when comparing against a prolate `Spheroid`, as above.

## Imperfect interfaces

`LayeredSpheroid` reuses [`LayeredSphere`](@ref)'s interface types:
[`KapitzaInterface`](@ref)`(ρ)` (a temperature-jump resistance — "LC",
low-conducting) and
[`SurfaceConductiveInterface`](@ref)`(β)` (a flux-jump surface
conductance — "HC", highly-conducting), on top of the default
`PerfectInterface`. Two exact limits anchor them:

```@example spheroid
s_lc = LayeredSpheroid((3.0,), (1.0,), (K1,); interfaces = (KapitzaInterface(1.0e10),))
s_insulated = LayeredSpheroid((3.0,), (1.0,), (TensISO{3}(1.0e-12),))

# an infinitely resistive interface behaves like a fully insulated core
get_array(gradient_gradient_loc(s_lc, K1, K0)) ≈
    get_array(gradient_gradient_loc(s_insulated, TensISO{3}(1.0e-12), K0))
```

## Multiple layers

Add layers by extending both semi-axis tuples (ascending, confocal)
and the per-layer conductivity/interface tuples:

```@example spheroid
focal2 = 3.0^2 - 1.0^2                  # shared focal distance² (outer layer)
a_in = 2.9
b_in = sqrt(a_in^2 - focal2)             # confocal inner layer

s2 = LayeredSpheroid(
    (a_in, 3.0), (b_in, 1.0), (TensISO{3}(5.0), K1);
    interfaces = (KapitzaInterface(0.1), PerfectInterface()),
)
layer_volume_fraction(s2, 1), layer_volume_fraction(s2, 2)
```

The convenience constructor
[`layered_spheroid_from_fractions`](@ref) builds the confocal
parameters for you from an outer aspect ratio, an outer size, and
target volume fractions:

```@example spheroid
s3 = layered_spheroid_from_fractions(3.0, 3.0, (0.3, 0.7), (TensISO{3}(5.0), K1))
layer_volume_fraction(s3, 1), layer_volume_fraction(s3, 2)
```

## In an RVE

`LayeredSpheroid` is a first-class phase geometry, exactly like
`Ellipsoid` or `LayeredSphere` — it carries **no** Hill tensor
(`is_homogeneous_inclusion(s) == false`), so it plugs into the schemes
through its own volume-averaged concentration tensors instead. The
declared phase property is accepted for signature compatibility but
ignored (the real moduli live in the layers):

```@example spheroid
rve = RVE(:M)
add_matrix!(rve, Ellipsoid(1.0, 1.0, 1.0), Dict(:K => K0))
add_phase!(rve, :I, s, Dict(:K => K1); fraction = 0.15)

Keff = homogenize(rve, MoriTanaka(), :K)
get_array(Keff)
```

## Local fields

The pointwise temperature, gradient and flux fields inside/around the
spheroid are available from
[`local_temperature`](@ref MeanFieldHom.LayeredSpheroids.local_temperature),
[`local_gradient`](@ref MeanFieldHom.LayeredSpheroids.local_gradient), and
[`local_flux`](@ref MeanFieldHom.LayeredSpheroids.local_flux), evaluated
at the spheroid's own confocal coordinates `(q, p, φ)` under a remote
axial and/or transverse loading:

```@example spheroid
using MeanFieldHom.LayeredSpheroids: local_temperature

qN = s.q[1]
local_temperature(s, K0, real(qN) + 1.0e-6, 0.5, 0.0; H_axial = 1.0)
```

See `scripts/35_spheroid_local_fields.jl` in the [Gallery](../gallery/index.md)
for a full meridian-plane temperature/flux map.
