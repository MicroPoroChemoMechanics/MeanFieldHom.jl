# Elliptic integrals — examples

`MeanFieldHom` re-exports the five entry points of the
[`MeanFieldHom.Elliptic`](@ref) submodule: `ell_K`, `ell_E`, `ell_F`,
`ell_RF`, `ell_RD`.

## Basic calls

```@example ell_basic
using MeanFieldHom

# Complete integrals, Float64 fast path
m = 0.5
K = ell_K(m)
E = ell_E(m)
(K, E)
```

```@example ell_basic
# The generic AGM path agrees with the Float64 fast path when promoted
K_big = ell_K(BigFloat(m))
abs(Float64(K_big) - K)
```

```@example ell_basic
# Incomplete integrals
φ = π/4
F  = ell_F(φ, m)
E2 = ell_E(φ, m)
(F, E2)
```

## Carlson symmetric forms

```@example ell_basic
# Unit-integration identity: R_F(1, 1, 1) = 1
ell_RF(1.0, 1.0, 1.0)
```

```@example ell_basic
# Homogeneity: R_F(λx, λy, λz) = λ^{-1/2} R_F(x, y, z)
λ = 3.0
λ^(-0.5) * ell_RF(1.0, 2.0, 3.0), ell_RF(λ, 2λ, 3λ)
```

## Automatic differentiation

The elliptic submodule is designed to flow through `ForwardDiff`. The
`Float64` fast path is automatically replaced by the generic
AGM / Carlson code whenever the input is a `ForwardDiff.Dual`:

```@example ell_ad
using MeanFieldHom, ForwardDiff

# Derivative of K(m) at m = 0.3
dKdm = ForwardDiff.derivative(ell_K, 0.3)

# Reference value: dK/dm = (E - (1-m)K) / (2m(1-m))
m = 0.3
ref = (ell_E(m) - (1 - m) * ell_K(m)) / (2m * (1 - m))
(dKdm, ref, abs(dKdm - ref))
```

```@example ell_ad
# Gradient of a scalar functional involving F(φ, m)
using ForwardDiff: gradient
f(x) = ell_F(x[1], x[2])^2 + ell_E(x[1], x[2])
g = gradient(f, [0.4, 0.5])
```

## Arbitrary precision

Using `BigFloat` gives access to the full precision of the AGM
recursion:

```@example ell_big
using MeanFieldHom

setprecision(BigFloat, 128)    # 128 mantissa bits
m = BigFloat("0.25")

K_big = ell_K(m)
K_f64 = ell_K(Float64(m))

# Round-trip error (BigFloat → Float64)
abs(Float64(K_big) - K_f64)
```

## Symbolic computation

With the optional `SymPy` extension loaded (the package declares
`SymPy` as a `[weakdeps]`), the complete integrals become symbolic
closed-form expressions:

```julia
using MeanFieldHom, SymPy

@syms m
ell_K(m)       # → K(m) as a SymPy expression
ell_E(π/4, m)  # → SymPy elliptic_e(π/4, m)
```

With `Symbolics` (no extension required — the pure-arithmetic AGM path
is used):

```julia
using MeanFieldHom, Symbolics

@variables m
expr = ell_K(m)
# Symbolics expression — use `simplify(expr)` to compact the output.
```

## Custom scalar backends

Adding a new backend for a user-defined `Number` type amounts to
adding methods on the public API:

```julia
struct MyScalar <: Number
    x::Float64
end

Base.float(s::MyScalar) = s.x

# Route K and E through the fast Float64 path
MeanFieldHom.Elliptic.ell_K(s::MyScalar) = ell_K(s.x)
MeanFieldHom.Elliptic.ell_E(s::MyScalar) = ell_E(s.x)
```

Downstream code that calls `ell_K`, `ell_E`, `ell_F`, `ell_RF`,
`ell_RD` will pick up the new methods automatically — multiple
dispatch does the rest.
