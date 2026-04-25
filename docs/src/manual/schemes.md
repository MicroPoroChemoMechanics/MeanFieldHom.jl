# Homogenisation schemes — user manual

The `MeanFieldHom.Schemes` module provides ten classical mean-field
homogenisation schemes plus a [`RVE`](@ref) container holding the
matrix and inclusion phases with their geometries, properties and
volume fractions or crack densities.

## Building an RVE

```julia
using MeanFieldHom, TensND

rve = RVE(:M)                                            # matrix is named :M
add_matrix!(rve, Ellipsoid(1.0), Dict(:C => TensISO{3}(30.0, 10.0)))
add_phase!(rve, :I, Ellipsoid(1.0, 1.0, 0.5),
           Dict(:C => TensISO{3}(60.0, 20.0)); fraction = 0.2)
add_phase!(rve, :CRACK, PennyCrack(1.0),
           Dict(:C => TensISO{3}(30.0, 10.0)); density  = 0.05)
```

Volume fractions are stored at the RVE level (not on the inclusions),
so a single inclusion remains usable for localisation-tensor
calculations without any RVE machinery
([`hill_tensor`](@ref), [`strain_strain_loc`](@ref), …). The matrix
volume fraction is implicit (`1 - Σ f_inc`); crack densities are
excluded from that sum.

## Calling a scheme

```julia
C_voigt = homogenize(rve, Voigt())         # type-instance
C_mt    = homogenize(rve, :mt)             # Symbol shortcut (lowercase canonical)
C_sc    = homogenize(rve, SelfConsistent(; abstol = 1e-12, maxiters = 200))
```

Every scheme takes the optional kwarg `property = :C` (default,
elasticity) or `property = :K` (conductivity). Iterative schemes also
accept `abstol`, `maxiters`, `damping`, `verbose`.

| Long form | Short / ECHOES code |
| --- | --- |
| `:voigt` | `:v`, `:V`, `:Voigt`, `:VOIGT` |
| `:reuss` | `:r`, `:R` … |
| `:dilute` | `:dil`, `:DIL` |
| `:dilute_dual` | `:dild`, `:DILD` |
| `:mori_tanaka` | `:mt`, `:MT` |
| `:maxwell` | `:max`, `:MAX` |
| `:ponte_castaneda_willis` | `:pcw`, `:PCW` |
| `:self_consistent` | `:sc`, `:SC` |
| `:asymmetric_self_consistent` | `:asc`, `:ASC` |
| `:differential` | `:diff`, `:DIFF` |

## Distribution shape (Maxwell, PCW)

```julia
rve = RVE(:M; distribution_shape = Ellipsoid(1.0, 1.0, 0.3))   # oblate outer
add_matrix!(rve, Ellipsoid(1.0), Dict(:C => TensISO{3}(30.0, 10.0)))
add_phase!(rve, :I, Ellipsoid(1.0), Dict(:C => TensISO{3}(60.0, 20.0));
           fraction = 0.3)
homogenize(rve, Maxwell())
```

The `distribution_shape` field is wrapped in `UniformDistribution`. A
future `PairwiseDistribution` (Willis 1982) can be added without
breaking the public API — see
[`AbstractDistributionShape`](@ref).

## Iterative solvers

```julia
homogenize(rve, SelfConsistent())                            # built-in damped Picard

# With NonlinearSolve.jl loaded:
using NonlinearSolve
homogenize(rve, SelfConsistent(; algorithm = NewtonRaphson(),
                                abstol = 1e-12, maxiters = 200))
```

`NewtonDefault()` reuses the SciML default Newton-Raphson once
`NonlinearSolve.jl` is loaded; without that load, it raises an explicit
error.

## Differential trajectories

```julia
homogenize(rve, DifferentialScheme(; nsteps = 200))                  # Proportional (default)
homogenize(rve, DifferentialScheme(; trajectory = Sequential([:I1, :I2])))
custom = CustomPath(Dict(:I => collect(range(0.0, 1.0; length = 101))))
homogenize(rve, DifferentialScheme(; trajectory = custom, nsteps = 100))
```

For multi-phase RVEs the trajectory choice is *physical* — the schemes
agree in the dilute limit and diverge at finite fractions. Cracks
(`CrackDensity`) are added at constant target density per step
regardless of the trajectory.

## Frequency-domain viscoelasticity

```julia
δ = 0.05
rve = RVE(:M)
add_matrix!(rve, Ellipsoid(1.0),
            Dict(:C => TensISO{3}(30.0 + δ * im, 10.0 + 0.5δ * im)))
add_phase!(rve, :I, Ellipsoid(1.0),
           Dict(:C => TensISO{3}(60.0 + δ * im, 20.0 + 0.5δ * im));
           fraction = 0.3)

C_eff = homogenize(rve, MoriTanaka())   # eltype(C_eff) == ComplexF64
```

All schemes propagate `Complex{Float64}` through their tensor algebra.
The `Im → 0` limit consistently recovers the real-modulus result.

## Sensitivity (ForwardDiff)

```julia
using ForwardDiff
df = ForwardDiff.derivative(0.3) do f
    rve = RVE(:M; T = typeof(f))
    add_matrix!(rve, Ellipsoid(1.0), Dict(:C => TensISO{3}(30.0, 10.0)))
    add_phase!(rve, :I, Ellipsoid(1.0), Dict(:C => TensISO{3}(60.0, 20.0));
               fraction = f)
    KM(homogenize(rve, MoriTanaka()))[1, 1]
end
```

Every scheme is differentiable through the fractions, moduli, and
inclusion geometry. The volume-fraction eltype `T` is fixed at RVE
construction (default `Float64`); pass `T = ForwardDiff.Dual{...}` or
`T = Complex{Float64}` explicitly when needed.
