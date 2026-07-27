# The differential scheme and path dependence

The differential (or incremental) scheme builds a composite the way
some real materials are actually made: by adding inclusions a little
at a time, re-homogenizing after every increment. That construction
history turns out to matter as soon as there is more than one
inclusion phase.

## The differential scheme

[`DifferentialScheme`](@ref) integrates the **Norris differential
equation** [norris1985](@cite) over a fictitious loading parameter
``\tau \in [0, 1]``:

```math
\frac{d\mathbb{C}^*}{d\tau} = \frac{1}{1-f(\tau)}\,\big(\mathbb{C}_i-\mathbb{C}^*\big):\mathbb{A}^*\big(\mathbb{C}^*\big)\,\frac{df}{d\tau},
```

starting from ``\mathbb{C}^*(\tau=0) = \mathbb{C}_0`` (the matrix) and
adding inclusions at each step *into the current effective medium*
``\mathbb{C}^*`` rather than into the original matrix. Each increment is
itself a dilute estimate — so the scheme is, in effect, an infinite
sequence of infinitesimal dilute corrections, integrated along a
**trajectory** `nsteps` steps long.

## Trajectories for multi-phase RVEs

With a single inclusion phase, that trajectory is unambiguous: fraction
grows monotonically from 0 to its target value. With **two or more**
inclusion phases, the *order* in which they are added along the way
becomes a modeling choice, exposed through the `trajectory` keyword:

- `Proportional()` (default) — every phase grows in lock-step,
  proportionally to its target fraction, at each step.
- `Sequential([:A, :B])` — phase `:A` is added to completion first, then
  `:B` is added into the resulting effective medium.
- `CustomPath(Dict(:A => τ -> ..., ...))` — an arbitrary fraction-of-τ
  schedule per phase.

## Seeing the path dependence

```@example tutdiff
using MeanFieldHom
using TensND
using LinearAlgebra
using Plots
gr()  # headless backend; GKSwstype is set to "100" in make.jl

C_m = iso_stiffness(30.0, 10.0)
C_stiff = iso_stiffness(60.0, 20.0)
C_soft = iso_stiffness(5.0, 2.0)

function build(f_total)
    r = RVE(:M)
    add_matrix!(r, Ellipsoid(1.0), Dict(:C => C_m))
    add_phase!(r, :STIFF, Ellipsoid(1.0), Dict(:C => C_stiff); fraction = f_total / 2)
    add_phase!(r, :SOFT, Ellipsoid(1.0), Dict(:C => C_soft); fraction = f_total / 2)
    return r
end

trajectories = [
    (Proportional(), "Proportional", :black),
    (Sequential([:STIFF, :SOFT]), "Sequential(STIFF, SOFT)", :blue),
    (Sequential([:SOFT, :STIFF]), "Sequential(SOFT, STIFF)", :red),
]

fs = collect(range(0.005, 0.3; length = 25))
plt = plot(;
    xlabel = "total inclusion fraction f₁+f₂", ylabel = "k_eff",
    legend = :topleft, framestyle = :box, size = (760, 480),
)
for (traj, label, color) in trajectories
    ks = [k_mu(homogenize(build(f), DifferentialScheme(; trajectory = traj, nsteps = 100), :C))[1] for f in fs]
    plot!(plt, fs, ks; label = label, color = color, lw = 2)
end
plt
```

The three curves coincide in the dilute limit (small total fraction,
where increments barely interact regardless of order) and fan out as
the total fraction grows. Adding the stiff phase first leaves it
embedded in a matrix that later softens further when the soft phase is
added; adding the soft phase first has the opposite effect. Neither
order is "more correct" in the abstract — the differential scheme
encodes a **construction history**, and for a real composite the
mixing or loading order is itself a physical modeling choice, not a
numerical one. This closes the elastic-scheme arc of these tutorials;
the next pages turn to cracks, viscoelasticity, and derivatives.
