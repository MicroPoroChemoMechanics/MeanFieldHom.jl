# Nonlinear solvers for the self-consistent fixed point

The [self-consistent scheme](02_bounds_and_schemes.md) is not a one-shot
formula — it is a fixed point ``\mathbb C = \mathrm{step}(\mathbb C)``,
solved iteratively. Every tutorial so far has used the package's
built-in solver, a damped Picard iteration
([`AndersonDefault`](@ref), the default). `MeanFieldHom` also ships a
dependency-free Newton-Raphson solver ([`NewtonDefault`](@ref)), and —
the subject of this page — a weak extension,
`MeanFieldHomNonlinearSolveExt`, that hands the same fixed point to any
algorithm from [NonlinearSolve.jl](https://github.com/SciML/NonlinearSolve.jl)
(`NewtonRaphson`, `TrustRegion`, `LevenbergMarquardt`, …).

This matters for two reasons: some SciML algorithms converge faster or
more robustly than Picard on stiff, high-contrast problems, and — since
`MeanFieldHom`'s sensitivities (the [previous](08_sensitivities.md) two
tutorials) are built on `ForwardDiff` — differentiating *through* an
external nonlinear solve must not silently break or, worse, silently
give the wrong answer. This page checks both, and closes with the
strength-criterion example from the [capstone tutorial](09_strength_criteria.md)
computed through a SciML solver instead of Picard.

## The three solver families

```@example tutnls
using MeanFieldHom
using TensND
using ForwardDiff
using NonlinearSolve
using LinearAlgebra
using Printf

const k_m, μ_m = 30.0, 10.0
const k_i, μ_i = 60.0, 20.0

rve = RVE(:M)
add_matrix!(rve, Ellipsoid(1.0), Dict(:C => TensISO{3}(k_m, μ_m)))
add_phase!(rve, :I, Ellipsoid(1.0), Dict(:C => TensISO{3}(k_i, μ_i)); fraction = 0.3)

C_picard = homogenize(rve, SelfConsistent())                                  # built-in damped Picard (default)
C_newton = homogenize(rve, SelfConsistent(; algorithm = NewtonDefault()))     # built-in Newton, dependency-free
C_nr     = homogenize(rve, SelfConsistent(; algorithm = NewtonRaphson()))     # NonlinearSolve.jl
C_tr     = homogenize(rve, SelfConsistent(; algorithm = TrustRegion()))       # NonlinearSolve.jl
C_auto   = homogenize(rve, SelfConsistent(; algorithm = AutoNonlinear()))     # auto-resolving

k_mu.((C_picard, C_newton, C_nr, C_tr, C_auto))
```

All five agree on the same fixed point — as they must, since they solve
the same equation. [`AutoNonlinear`](@ref) resolves, at runtime, to a
globalized SciML algorithm (`TrustRegion`) when the extension is active,
and to the built-in [`NewtonDefault`](@ref) otherwise — a solver choice
that works whether or not `NonlinearSolve.jl` happens to be loaded.
`AsymmetricSelfConsistent` accepts `algorithm` the same way, for both
its stiffness- and compliance-form branches.

!!! note "Why `AndersonDefault` stays the default"
    A root-finder is not guaranteed to track the *physical* branch of
    the self-consistent equation through the porous-percolation
    bifurcation the way Picard's positive-definite guard and
    `select_best` do (see the [porous benchmark tutorial](04_porous_benchmark.md)).
    `AutoNonlinear` and explicit SciML algorithms are opt-in for that
    reason — reach for them on well-conditioned, high-contrast problems
    away from a bifurcation, where they can be markedly faster.

## Benchmark: time and memory

A single `homogenize` call is cheap enough that the interesting
comparison is many repeated solves — exactly the situation in a
parameter sweep or an optimization loop. The numbers below use plain
`@elapsed`/`@allocated` (minimum of a few samples after warm-up); they
are illustrative of the *relative* cost between solvers on this
particular problem, not absolute guarantees — see
`scripts/bench/bench_sc_solvers.jl` for the full, repeatable benchmark
this page summarizes.

```@example tutnls
function bench(f, label; warmups = 3, samples = 5)
    for _ in 1:warmups
        f()
    end
    t_min, bytes = Inf, 0
    for _ in 1:samples
        b = @allocated f()
        t = @elapsed f()
        if t < t_min
            t_min, bytes = t, b
        end
    end
    return (; label, t_ms = t_min * 1.0e3, mib = bytes / 2^20)
end

results = [
    bench(() -> homogenize(rve, SelfConsistent()), "Picard"),
    bench(() -> homogenize(rve, SelfConsistent(; algorithm = NewtonDefault())), "Newton (built-in)"),
    bench(() -> homogenize(rve, SelfConsistent(; algorithm = NewtonRaphson())), "NewtonRaphson (NLS)"),
    bench(() -> homogenize(rve, SelfConsistent(; algorithm = TrustRegion())), "TrustRegion (NLS)"),
]
for r in results
    @printf "%-22s  %8.3f ms   %8.4f MiB\n" r.label r.t_ms r.mib
end
```

On this well-conditioned, moderate-contrast problem, all four solvers
converge in a handful of iterations; differences in time/memory here
mostly reflect iteration count and per-iterate overhead (Picard has the
cheapest iterate but sometimes needs more of them; Newton-type solvers
converge quadratically but pay for a Jacobian each step). Which solver
wins depends on the contrast and proximity to a bifurcation — this is
exactly why the choice is exposed as a keyword rather than hard-coded.

## Differentiating through a SciML solve

This is the part that needs care. `ForwardDiff` computes
`derivative(rve, scheme, p)` by seeding a `Dual` number into the RVE and
re-running `homogenize` — so if `scheme` solves its fixed point via
`NonlinearSolve.jl`, that solver internally sees `Dual`-valued inputs
too. Naively handing those to a general-purpose nonlinear solver risks
**nested** `ForwardDiff.Dual`s: the solver's own Jacobian routine seeds
*its own* `Dual` on top of the caller's, which is fragile (tag ordering)
and wasteful.

`MeanFieldHomNonlinearSolveExt` avoids this with an
implicit-function-theorem (IFT) **lift**: it solves the *primal*
problem — every input stripped to `Float64` via `ForwardDiff.value`,
with an *explicit* finite-difference Jacobian so `NonlinearSolve` never
seeds a `Dual` of its own — and then recovers the caller's partials with
one linear-algebra correction,

```math
p^\star_{\text{dual}} = p^\star - \Big(\frac{\partial F}{\partial p}\Big)^{-1} F(p^\star; \theta),
```

where ``p^\star`` is the primal root and ``\theta`` the differentiated
parameter — exactly the implicit function theorem applied once, at the
root. Only the *caller's* `Dual` tag is ever present; nothing is nested.
The built-in [`NewtonDefault`](@ref) instead differentiates straight
through its own iterations (safe because it is entirely
hand-written), so both paths give the same answer by construction —
this section checks that they do:

```@example tutnls
idxC = C -> get_array(C)[1, 1, 1, 1]

d_picard = derivative(rve, SelfConsistent(), property(:I, :C, :bulk); indexer = idxC)
d_newton = derivative(rve, SelfConsistent(; algorithm = NewtonDefault()), property(:I, :C, :bulk); indexer = idxC)
d_nr     = derivative(rve, SelfConsistent(; algorithm = NewtonRaphson()), property(:I, :C, :bulk); indexer = idxC)
d_tr     = derivative(rve, SelfConsistent(; algorithm = TrustRegion()), property(:I, :C, :bulk); indexer = idxC)

(d_picard, d_newton, d_nr, d_tr)
```

```@example tutnls
# Central finite difference — a solver-independent ground truth.
function f_modulus(K_I)
    r = RVE(:M)
    add_matrix!(r, Ellipsoid(1.0), Dict(:C => TensISO{3}(k_m, μ_m)))
    add_phase!(r, :I, Ellipsoid(1.0), Dict(:C => TensISO{3}(K_I, μ_i)); fraction = 0.3)
    return idxC(homogenize(r, SelfConsistent()))
end
h = 1.0e-5
d_fd = (f_modulus(k_i + h) - f_modulus(k_i - h)) / (2h)
@printf "central FD = %.6f   |TrustRegion − FD| / |FD| = %.3e\n" d_fd abs(d_tr - d_fd) / abs(d_fd)
```

All four agree to within the solver tolerances — the IFT lift is exact
to first order, and cheaper than differentiating through iterations
(one linear solve at the root, versus propagating partials at every
Picard/Newton step). The same check holds differentiating with respect
to the *inclusion* fraction or the *matrix* shear modulus — parameters
that live on a phase other than the one the fixed point's initial
estimate is built from, which is exactly the case a naive
type-from-the-initial-guess dispatch can miss:

```@example tutnls
d_f_picard = derivative(rve, SelfConsistent(), amount(:I); indexer = idxC)
d_f_tr     = derivative(rve, SelfConsistent(; algorithm = TrustRegion()), amount(:I); indexer = idxC)
d_m_picard = derivative(rve, SelfConsistent(), property(:M, :C, :shear); indexer = idxC)
d_m_tr     = derivative(rve, SelfConsistent(; algorithm = TrustRegion()), property(:M, :C, :shear); indexer = idxC)

(; d_f_picard, d_f_tr, d_m_picard, d_m_tr)
```

## Strength criterion, revisited

The [capstone tutorial](09_strength_criteria.md) built a macroscopic
strength ellipse for a porous solid entirely from `ForwardDiff`
derivatives of `(k_{\mathrm{hom}}, \mu_{\mathrm{hom}})` with respect to
the solid's own shear modulus — no closed-form derivative written by
hand. Swapping the underlying SC solve for a SciML algorithm changes
nothing about that recipe, since `derivative` only ever sees
`homogenize`'s public interface:

```@example tutnls
const k_s, μs_value = 1.0e6, 1.0
const TINY = 1.0e-12
const ω_aspect = 0.1
const φ_value = 0.15

function C_hom_iso_2vec(μs::T, scheme) where {T}
    r = RVE(:SOLID; T = T)
    add_matrix!(r, Spheroid(ω_aspect), Dict(:C => iso_stiffness(convert(T, k_s), μs)); symmetrize = IsoSymmetrize())
    add_phase!(r, :PORE, Spheroid(ω_aspect), Dict(:C => iso_stiffness(convert(T, TINY), convert(T, TINY)));
               fraction = convert(T, φ_value), symmetrize = IsoSymmetrize())
    C = homogenize(r, scheme, :C)
    return [k_mu(best_fit_iso(C))...]
end

function ellipse_AB(scheme)
    k_hom, μ_hom = C_hom_iso_2vec(μs_value, scheme)
    dk_dμs, dμ_dμs = ForwardDiff.derivative(μ -> C_hom_iso_2vec(μ, scheme), μs_value)
    A = (μs_value / k_hom)^2 * dk_dμs
    B = (μs_value / μ_hom)^2 * dμ_dμs
    return A, B
end

for (scheme, label) in (
        (SelfConsistent(; abstol = 1.0e-10, maxiters = 300, select_best = true), "Picard"),
        (SelfConsistent(; algorithm = TrustRegion()), "TrustRegion"),
    )
    A, B = ellipse_AB(scheme)
    @printf "%-12s  A = %.6g   B = %.6g\n" label A B
end
```

The two solver families produce the same ellipse — the strength
criterion is a property of the *scheme*, not of how its fixed point
happens to be solved. `AutoNonlinear` (or an explicit SciML algorithm)
is a drop-in accelerator when a sweep over many porosities or contrasts
makes Picard's iteration count add up; `AndersonDefault` remains the
right default whenever the sweep crosses a percolation-like
bifurcation, as in the [porous benchmark](04_porous_benchmark.md).
