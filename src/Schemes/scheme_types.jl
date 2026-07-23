# =============================================================================
#  scheme_types.jl — concrete homogenization-scheme types and the
#  differential-trajectory hierarchy.
#
#  Each scheme is a subtype of `HomogenizationScheme`; bounds and one-shot
#  schemes are singletons (`Voigt`, `MoriTanaka`, …) while the iterative ones
#  (`SelfConsistent`, `AsymmetricSelfConsistent`, `DifferentialScheme`) carry
#  configuration. Solver algorithms for the self-consistent family are
#  themselves marker types (`AndersonDefault`, `NewtonDefault`, plus any
#  algorithm provided by the SciML `NonlinearSolve` weak extension).
# =============================================================================

"""
    HomogenizationScheme

Supertype for every mean-field homogenization scheme. Concrete subtypes:

- bounds — [`Voigt`](@ref), [`Reuss`](@ref) ;
- one-shot with matrix — [`Dilute`](@ref), [`DiluteDual`](@ref),
  [`MoriTanaka`](@ref), [`Maxwell`](@ref), [`PonteCastanedaWillis`](@ref) ;
- iterative — [`SelfConsistent`](@ref), [`AsymmetricSelfConsistent`](@ref) ;
- trajectory-based — [`DifferentialScheme`](@ref).
"""
abstract type HomogenizationScheme end

# ── Bounds ───────────────────────────────────────────────────────────────────

"""
    Voigt() <: HomogenizationScheme

Voigt (uniform-strain) upper bound: ``\\langle \\mathbb C \\rangle``.
"""
struct Voigt <: HomogenizationScheme end

"""
    Reuss() <: HomogenizationScheme

Reuss (uniform-stress) lower bound: ``\\langle \\mathbb S \\rangle^{-1}``.
"""
struct Reuss <: HomogenizationScheme end

# ── One-shot schemes (require a matrix phase) ────────────────────────────────

"""
    Dilute() <: HomogenizationScheme

Dilute scheme: ``\\mathbb C_{\\mathrm{eff}} = \\mathbb C_0 + \\sum_i f_i \\mathbb N_i``
where ``\\mathbb N_i = (\\mathbb C_i - \\mathbb C_0):\\mathbb A_{\\varepsilon\\varepsilon}^{(i)}``
is the size-independent stiffness contribution
([Eshelby 1957](@cite eshelby1957);
[Kachanov & Sevostianov 2018](@cite kachanov2018)).
"""
struct Dilute <: HomogenizationScheme end

"""
    DiluteDual() <: HomogenizationScheme

Dual dilute scheme on the compliance:
``\\mathbb S_{\\mathrm{eff}} = \\mathbb S_0 + \\sum_i f_i \\mathbb H_i``,
returning ``\\mathbb C_{\\mathrm{eff}} = \\mathbb S_{\\mathrm{eff}}^{-1}``.
"""
struct DiluteDual <: HomogenizationScheme end

"""
    MoriTanaka() <: HomogenizationScheme

Mori-Tanaka scheme ([Mori & Tanaka 1973](@cite mori1973);
[Christensen 1990](@cite christensen1990)).
"""
struct MoriTanaka <: HomogenizationScheme end

"""
    Maxwell() <: HomogenizationScheme

Maxwell homogenization, using the RVE's distribution shape as the
reference for the Hill polarisation tensor.
"""
struct Maxwell <: HomogenizationScheme end

"""
    PonteCastanedaWillis() <: HomogenizationScheme

Ponte-Castañeda & Willis 1995 scheme — distribution-shape-aware
generalization of Mori-Tanaka.
"""
struct PonteCastanedaWillis <: HomogenizationScheme end

# ── Self-consistent solvers (built-in markers) ───────────────────────────────

"""
    AndersonDefault()

Marker selecting the built-in Anderson-accelerated fixed-point solver
(default for [`SelfConsistent`](@ref)).  Pure Julia, Dual-safe.
"""
struct AndersonDefault end

"""
    NewtonDefault()

Marker selecting the built-in Newton-Raphson solver with ForwardDiff
Jacobian (alternative to [`AndersonDefault`](@ref)).
"""
struct NewtonDefault end

"""
    SelfConsistent(; algorithm = AndersonDefault(), kwargs...) <: HomogenizationScheme

Self-consistent scheme. The `algorithm` selects the non-linear solver;
default is the built-in Anderson acceleration. Any solver from the
SciML `NonlinearSolve.jl` package can be passed once `using NonlinearSolve`
activates the weak extension `MeanFieldHomNonlinearSolveExt`.

Standard kwargs forwarded to the solver: `abstol`, `reltol`,
`maxiters`, `damping`, `verbose`.
"""
struct SelfConsistent{A, K <: NamedTuple} <: HomogenizationScheme
    algorithm::A
    options::K
end
SelfConsistent(; algorithm = AndersonDefault(), kwargs...) =
    SelfConsistent(algorithm, NamedTuple(kwargs))

"""
    AsymmetricSelfConsistent(; algorithm = AndersonDefault(), kwargs...) <: HomogenizationScheme

Asymmetric self-consistent scheme: iterates in stiffness or compliance
space depending on the matrix-vs-Voigt-bound contrast, providing a
better behavior than [`SelfConsistent`](@ref) in matrix-stiff /
inclusion-soft regimes.
"""
struct AsymmetricSelfConsistent{A, K <: NamedTuple} <: HomogenizationScheme
    algorithm::A
    options::K
end
AsymmetricSelfConsistent(; algorithm = AndersonDefault(), kwargs...) =
    AsymmetricSelfConsistent(algorithm, NamedTuple(kwargs))

# ── Differential scheme + trajectories ───────────────────────────────────────

"""
    DifferentialTrajectory

Supertype describing the path through the multi-phase volume-fraction
space used by the [`DifferentialScheme`](@ref) scheme. Concrete subtypes:

- [`Proportional`](@ref) (default) — every phase grows linearly with
  the fictitious incorporation time `τ ∈ [0, 1]`, all phases reaching
  their target simultaneously.
- [`Sequential`](@ref) — phases are introduced one after the other in
  the user-supplied order, each occupying a contiguous slice of `τ`.
- [`CustomPath`](@ref) — explicit per-phase trajectory as a vector of
  monotone non-decreasing values; piecewise-linear interpolated along
  `τ ∈ [0, 1]`.
- [`Path`](@ref) — explicit per-phase trajectory as a callable
  `τ -> f(τ)` (auto-differentiated by `ForwardDiff`); the natural API
  for the multi-phase incorporation-sequence ODE
  ([Norris 1985](@cite norris1985); the user's hand-written DEM note).
"""
abstract type DifferentialTrajectory end

"""
    Proportional() <: DifferentialTrajectory

All phases grow proportionally during the differential integration.
"""
struct Proportional <: DifferentialTrajectory end

"""
    Sequential(order::Vector{Symbol}) <: DifferentialTrajectory

Introduce phases in the given `order`. The first phase ramps from 0 to
its target fraction over the steps it owns, then is frozen; the next
phase ramps over its allotted steps; and so on.
"""
struct Sequential <: DifferentialTrajectory
    order::Vector{Symbol}
end

"""
    CustomPath(path::Dict{Symbol, <:AbstractVector{<:Real}}) <: DifferentialTrajectory

Explicit per-phase trajectory. `path[:phase]` must be a length-`N`
monotone vector with `path[:phase][1] = 0` and `path[:phase][end] = 1`,
where `N` is the number of differential steps.
"""
struct CustomPath{D <: AbstractDict{Symbol, <:AbstractVector{<:Real}}} <: DifferentialTrajectory
    path::D
end

"""
    Path(path::Dict{Symbol, <:Function}) <: DifferentialTrajectory

Explicit per-phase trajectory as a callable.  `path[:phase]` is a
function of the fictitious incorporation time `τ ∈ [0, 1]` returning
the **effective volume fraction ratio** `f_α(τ) / f_α^∞ ∈ [0, 1]` for
solid phases, or the **density ratio** `ε_α(τ) / ε_α^∞` for crack
phases — `f(0) = 0`, `f(1) = 1`, monotone non-decreasing.

The derivative `df/dτ` is computed by `ForwardDiff.derivative` at
each ODE step.  This is the natural API for the multi-phase DEM
incorporation-sequence ODE :

```math
\\frac{\\mathrm d \\mathbb C^{hom}}{\\mathrm d \\tau}
  = \\sum_i \\frac{\\mathrm d \\varphi_i}{\\mathrm d \\tau}
            (\\mathbb C_i - \\mathbb C^{hom}):\\mathbb A_i^{dil}(\\mathbb C^{hom})
```

with the volumetric balance `dφ = (𝟙 - f ⊗ 𝐔)^{-1} · df` (Sherman-
Morrison) inverted at each `τ` to translate user-supplied `f_α(τ)`
into the increments `dφ_α / dτ`.

The single-phase case is degenerate (`f₁` itself serves as `τ`) and
does not require a `Path` — the default [`Proportional`](@ref) is
sufficient.
"""
struct Path{D <: AbstractDict{Symbol, <:Function}} <: DifferentialTrajectory
    path::D
end

"""
    DifferentialScheme(; trajectory = Proportional(), nsteps::Int = 100,
                         abstol::Real = 1e-8, reltol::Real = 1e-6,
                         alg = nothing, kwargs...)

Differential scheme : integrates the Norris ODE on the fictitious
incorporation time `τ ∈ [0, 1]` ([Norris 1985](@cite norris1985)) :

```math
\\frac{\\mathrm d \\mathbb C^{hom}}{\\mathrm d \\tau}
  = \\sum_i \\frac{\\mathrm d \\varphi_i}{\\mathrm d \\tau}\\,
            (\\mathbb C_i - \\mathbb C^{hom}):\\mathbb A_i^{dil}(\\mathbb C^{hom})
```

with the volume balance `df = (𝟙 - f ⊗ 𝐔)·dφ` inverted by Sherman-
Morrison so the user supplies effective volume fractions `f_α(τ)`
along the chosen `trajectory`.

# Keyword arguments

- `trajectory` — one of [`Proportional`](@ref), [`Sequential`](@ref),
  [`CustomPath`](@ref), [`Path`](@ref).  Default `Proportional()`.
- `nsteps` — density of save points along `τ` (passed as `saveat` to
  the SciML ODE solver).  The integration step is controlled by
  `abstol` / `reltol`, **not** by `nsteps`.
- `abstol`, `reltol` — ODE solver tolerances (forwarded to
  `OrdinaryDiffEq.solve`).
- `alg` — explicit ODE algorithm.  `nothing` selects `Tsit5()` (5th
  order adaptive Runge-Kutta).  Pass any `OrdinaryDiffEqAlgorithm`
  instance to override (e.g. `Vern9()` for higher accuracy, `Rosenbrock23()`
  for stiff cases).
- `kwargs...` — additional kwargs stored verbatim and forwarded to the
  scheme-specific dispatch.
"""
struct DifferentialScheme{P <: DifferentialTrajectory, K <: NamedTuple} <: HomogenizationScheme
    trajectory::P
    options::K
end
DifferentialScheme(;
    trajectory = Proportional(),
    nsteps::Int = 100,
    abstol::Real = 1.0e-8,
    reltol::Real = 1.0e-6,
    alg = nothing,
    kwargs...
) =
    DifferentialScheme(
    trajectory,
    (; nsteps, abstol, reltol, alg, kwargs...)
)
