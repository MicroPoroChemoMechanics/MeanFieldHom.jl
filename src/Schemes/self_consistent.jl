# =============================================================================
#  self_consistent.jl — SelfConsistent + AsymmetricSelfConsistent.
#
#  Iterates `C^{n+1} = step(C^n)` where `step` is a Mori-Tanaka-like
#  evaluation that uses C^n as the reference matrix. The dispatcher
#  `_solve_sc` picks the non-linear solver:
#
#   * `AndersonDefault`            — built-in damped Picard fixed point.
#                                    Pure Julia, Dual-safe. Default.
#   * `NewtonDefault`              — placeholder for the SciML weak
#                                    extension (`MeanFieldHomNonlinearSolveExt`).
#                                    Triggers an explicit error if invoked
#                                    without `using NonlinearSolve`.
#   * any algorithm from
#     `NonlinearSolve.jl`          — handled by the weak extension.
#
#  A future native Anderson with memory > 1 will replace the current
#  `AndersonDefault` (currently Picard with relaxation, equivalent to
#  Anderson with memory = 1).
# =============================================================================

# ── Public _evaluate ────────────────────────────────────────────────────────

"""
    _evaluate(rve, sc::SelfConsistent, ::Val{p}; kw...) -> AbstractTens

Self-consistent scheme for property `:p`
([McLaughlin 1977](@cite mclaughlin1977)). Iterates

```math
\\mathbb C^{(n+1)} = \\Big(\\sum_i f_i\\,\\mathbb C_i \\!:\\! \\mathbb A_\\mathrm{dil}^{(i)}(\\mathbb C^{(n)})\\Big)
                     :\\Big(\\sum_i f_i\\,\\mathbb A_\\mathrm{dil}^{(i)}(\\mathbb C^{(n)})\\Big)^{-1}
```

with the dilute concentration tensor evaluated against the current
estimate `C^{(n)}` itself (rather than the matrix property).

The solver algorithm is selected by `sc.algorithm`; convergence kwargs
in `sc.options` (`abstol`, `maxiters`, `damping`, `verbose`) override
their defaults. External algorithms from `NonlinearSolve.jl` are
supported via the weak extension `MeanFieldHomNonlinearSolveExt`.

Cracks are not natively handled by this stiffness-form SC (the strain
concentration tensor is singular). For mixed RVEs (solid + crack) use
[`AsymmetricSelfConsistent`](@ref) instead.
"""
function _evaluate(rve::RVE, sc::SelfConsistent, ::Val{p}; kw...) where {p}
    P_init = matrix_property(rve, p)
    solver_kw, step_kw = _split_sc_kwargs(kw)
    step = C -> _sc_step(rve, C, p; step_kw...)
    return _solve_sc(sc.algorithm, step, P_init; sc.options..., solver_kw...)
end

# Solver-only kwargs are intercepted at this level and never forwarded
# to the underlying `_sc_step` (which would otherwise leak them down to
# `hill_tensor` and trigger a `MethodError` on the unknown kwarg).
const _SC_SOLVER_KWARGS = (:abstol, :reltol, :maxiters, :damping,
                            :verbose, :select_best)

function _split_sc_kwargs(kw)
    solver = Dict{Symbol, Any}()
    step   = Dict{Symbol, Any}()
    for (k, v) in kw
        if k in _SC_SOLVER_KWARGS
            solver[k] = v
        else
            step[k] = v
        end
    end
    return pairs(NamedTuple(solver)), pairs(NamedTuple(step))
end

# ── SC step (one Mori-Tanaka-like iterate against current estimate) ─────────

function _sc_step(rve::RVE, C_n, prop::Symbol; kw...)
    return _sc_step_dispatch(rve, C_n, prop; kw...)
end

# 4th-order — symmetric (Hill 1965 / Budiansky 1965) self-consistent
# iteration : all phases (matrix included) carry a non-trivial dilute strain
# concentration A_α = inv(I + P(C_α - C_n)) computed in the iterating
# effective medium C_n. This gives the textbook SC fixed point with
# Hashin-Shtrikman lower-bound percolation behaviour for porous media.
function _sc_step_dispatch(rve::RVE, C_n::TensND.AbstractTens{4, 3}, prop::Symbol;
                           kw...)
    A_avg  = zero(C_n)
    CA_avg = zero(C_n)
    # Matrix and inclusion phases enter the SC sum on equal footing.
    for name in rve.phase_names
        if name === rve.matrix_name
            f = matrix_volume_fraction(rve)
        else
            a = rve.amounts[name]
            a isa VolumeFraction || continue   # cracks ignored (use ASC)
            f = amount_value(a)
        end
        P_α = phase_property(rve, name, prop)
        A_dil = _phase_dilute_concentration(rve, name, prop, C_n; kw...)
        A_avg  += f * A_dil
        CA_avg += f * (P_α ⊡ A_dil)
    end
    return CA_avg ⊡ inv(A_avg)
end

# 2nd-order — same symmetric SC structure for conductivity / diffusion.
function _sc_step_dispatch(rve::RVE, K_n::TensND.AbstractTens{2, 3}, prop::Symbol;
                           kw...)
    A_avg  = zero(K_n)
    KA_avg = zero(K_n)
    for name in rve.phase_names
        if name === rve.matrix_name
            f = matrix_volume_fraction(rve)
        else
            a = rve.amounts[name]
            a isa VolumeFraction || continue
            f = amount_value(a)
        end
        P_α = phase_property(rve, name, prop)
        A_dil = _phase_dilute_concentration(rve, name, prop, K_n; kw...)
        A_avg  += f * A_dil
        KA_avg += f * (P_α ⋅ A_dil)
    end
    return KA_avg ⋅ inv(A_avg)
end

# ── Built-in solvers ────────────────────────────────────────────────────────

"""
    _solve_sc(algo, step, x0; abstol, reltol, maxiters, damping, verbose,
              select_best, kw...) -> AbstractTens

Generic solver dispatcher for SC fixed points. Built-in:

- [`AndersonDefault`](@ref) — Picard with relaxation
  (`x_{n+1} = (1-damping)·step(x_n) + damping·x_n`). `damping = 0.0`
  default; raise to ≈ 0.5 for high-contrast iterations that overshoot.
  Convergence near a bifurcation (e.g. SC at the porous-percolation
  threshold) is intrinsically slow because the Picard Jacobian
  eigenvalue approaches 1 there; in that regime, set
  `select_best = true` to return the best iterate observed during the
  loop, or load `NonlinearSolve.jl` and switch to Newton/Anderson via
  the `algorithm` keyword.
- [`NewtonDefault`](@ref) — SciML Newton-Raphson, available only when
  the weak extension `MeanFieldHomNonlinearSolveExt` is loaded
  (`using NonlinearSolve`).

Other algorithms from `NonlinearSolve.jl` are supported through the
same weak extension.

Convergence is declared when `‖x_new − x_old‖ ≤ abstol + reltol · ‖x_old‖`
(absolute *and* relative tolerance, additive convention; pass
`abstol = 0` to require purely relative convergence). Default values:
`abstol = 1e-12`, `reltol = 1e-8`.

When `select_best = true`, the solver tracks the best iterate seen
during the loop (smallest residual on the value field) and returns it
at the end. Useful for high-contrast iterations where Picard
oscillates around the fixed point: the *last* iterate may be worse
than an earlier one. Default is `false` (return last iterate).

Non-convergence is reported via `@debug` (silent by default; set
`JULIA_DEBUG=MeanFieldHom` to surface it) rather than `@warn`. Near
bifurcation points the Picard step intrinsically slows down (the
linearised step has a Jacobian eigenvalue ≈ 1) and the residual stalls
above `tol_eff` while still being negligibly small compared to the
matrix-property scale; the returned iterate is informative even when
the strict tolerance is not reached, so a default warning would be
noise.
"""
function _solve_sc(::AndersonDefault, step, x0::TensND.AbstractTens;
                   abstol::Real = 1.0e-12, reltol::Real = 1.0e-8,
                   maxiters::Int = 100,
                   damping::Real = 0.0, verbose::Bool = false,
                   select_best::Bool = false,
                   kw...)
    x = x0
    last_resid = _sc_residual_zero(x0)
    x_best = x0
    resid_best_val = typemax(_value_eltype(x0))
    for k in 1:maxiters
        x_new = step(x)
        last_resid = _sc_residual_norm(x_new, x)
        norm_x  = _sc_residual_norm(x, zero(x))
        tol_eff = abstol + reltol * _scalar_value(norm_x)
        verbose && @info "SC iter $k : ‖Δ‖ = $last_resid   tol = $tol_eff"
        if select_best
            v = _scalar_value(last_resid)
            if v < resid_best_val
                resid_best_val = v
                x_best = x_new
            end
        end
        _sc_converged(last_resid, tol_eff) && return x_new
        x = (one(real(eltype(x))) - damping) * x_new + damping * x
    end
    # Non-convergence is reported as a `@debug` message rather than a
    # `@warn` so it stays out of normal output. Set
    # `JULIA_DEBUG=MeanFieldHom` (or pass `verbose = true`) to surface
    # the diagnostics. Near bifurcation points (porous SC at
    # percolation, …) the Picard step intrinsically slows down and
    # `last_resid` may stall above the requested tolerance while
    # remaining negligible compared to the matrix-property scale; the
    # returned iterate is still informative.
    @debug "SC (AndersonDefault/Picard) did not reach tolerance" maxiters last_resid abstol reltol
    return select_best ? x_best : x
end

function _solve_sc(::NewtonDefault, step, x0::TensND.AbstractTens; kw...)
    error("NewtonDefault requires NonlinearSolve.jl: load it with " *
          "`using NonlinearSolve` to activate the MeanFieldHomNonlinearSolveExt " *
          "extension. Default solver `AndersonDefault` works without any extra dependency.")
end

# Residual norm operating on stored components (works for any tensor type).
_sc_residual_norm(a::TensND.AbstractTens, b::TensND.AbstractTens) =
    sqrt(sum(abs2, get_array(a) .- get_array(b)))

# Initial zero residual with eltype matching `x0` (Dual-preserving).
_sc_residual_zero(x0::TensND.AbstractTens) = zero(real(eltype(x0)))

# Scalar Float64 view of a residual: extract the `.value` field for
# `ForwardDiff.Dual` (via duck-typing — no hard dependency on ForwardDiff)
# so that comparisons and best-iterate tracking work uniformly.
_scalar_value(r::Real) = float(r)
_scalar_value(r) = hasfield(typeof(r), :value) ? float(_scalar_value(getfield(r, :value))) :
                   throw(ArgumentError("cannot reduce residual of type $(typeof(r)) to a Float64"))

# Float64 type of `_scalar_value(zero_like_eltype(x0))`. Used to pick a
# safe `typemax` for `resid_best_val` regardless of whether x0 is real
# or Dual. ForwardDiff.Dual is `<: Real` so we can't dispatch on `T<:Real`
# alone — we must inspect the `:value` field first, only falling through
# to `Real` for plain numeric eltypes.
_value_eltype(::TensND.AbstractTens{<:Any, <:Any, T}) where {T} = _value_eltype(T)
function _value_eltype(::Type{T}) where {T}
    hasfield(T, :value) && return _value_eltype(fieldtype(T, :value))
    T <: Real && return float(T)
    return Float64
end

# Convergence criterion. For Real residuals it is the obvious `r < abstol`.
# For Dual residuals (ForwardDiff), require both the value AND every
# partial derivative to be below `abstol`: a fixed-point iteration where
# ‖value‖ has converged but ‖partials‖ has not gives a derivative that
# is numerically wrong (each Picard step propagates a contraction
# coefficient of order ‖∂step/∂x‖ — partials converge as fast as the
# value only if the contraction is well-behaved, which is not guaranteed
# for ill-conditioned SC iterations).
_sc_converged(r::Real, abstol::Real) = r < abstol
function _sc_converged(r, abstol::Real)
    hasfield(typeof(r), :value) || return float(r) < abstol
    abs(_scalar_value(getfield(r, :value))) < abstol || return false
    if hasfield(typeof(r), :partials)
        p = getfield(r, :partials)
        # `partials` is a `ForwardDiff.Partials{N,T}` wrapper; iterate values.
        for v in p
            _sc_converged(v, abstol) || return false
        end
    end
    return true
end

# =============================================================================
#  AsymmetricSelfConsistent — Mori-Tanaka-style iteration with self-reference.
#
#  At each iteration, the dilute concentration tensor is computed in the
#  current effective medium (stiffness `C^n` or compliance `S^n`), but
#  the resulting dilute correction is added to the matrix property
#  (`C_m` or `S_m`) — not to the iterating estimate. This is the
#  "asymmetric" formulation: the matrix retains its privileged role
#  even in the SC iteration. The fixed point coincides with the
#  Hill-symmetric SC fixed point, but the iteration dynamics differ
#  (different basins of attraction).
#
#  Stiffness form  : C^{n+1} = C_m + Σ_i f_i (C_i − C_m) A_{εε,i}^{(C^n)}
#  Compliance form : S^{n+1} = S_m + (⟨A_{εε}⟩ − S_m ⟨C_i A_{εε}⟩) S^n
#
#  The branch is picked by the matrix-vs-Voigt-bound contrast: when the
#  matrix is much stiffer than the Voigt upper bound (matrix-stiff /
#  inclusions-soft, e.g. a porous medium in elasticity), the compliance
#  form is contractive; otherwise stiffness form. Selection follows the
#  C++ reference's squared-Frobenius-norm criterion.
# =============================================================================

"""
    _evaluate(rve, asc::AsymmetricSelfConsistent, ::Val{p}; kw...) -> AbstractTens

Asymmetric self-consistent scheme. Iterates a Mori-Tanaka-like update
where the dilute concentration tensors use the *current effective
medium* as reference but the dilute correction is added to the
*matrix property*. The iteration dynamics differ from the Hill
symmetric SC even though the fixed point is the same — for porous and
crack RVEs the asymmetric form converges to a different physical
branch (the matrix-distinguished branch).

Compliance form is selected when the matrix property is "stiffer"
(in squared Frobenius norm) than the Voigt upper bound — the
contractivity argument otherwise breaks down for high-contrast
matrix-stiff systems.
"""
function _evaluate(rve::RVE, asc::AsymmetricSelfConsistent, ::Val{p}; kw...) where {p}
    if _asc_use_stiffness(rve, p; kw...)
        return _asc_iterate_stiffness(rve, asc, Val(p); kw...)
    else
        return _asc_iterate_compliance(rve, asc, Val(p); kw...)
    end
end

# Squared Frobenius norm — matches the C++ reference selection.
function _asc_use_stiffness(rve::RVE, prop::Symbol; kw...)
    P₀      = matrix_property(rve, prop)
    P_voigt = _evaluate(rve, Voigt(), Val(prop); kw...)
    return _frob_sq(P_voigt) ≥ _frob_sq(P₀)
end

_frob_sq(t::TensND.AbstractTens) = sum(abs2, get_array(t))

# ── Stiffness-form iteration ────────────────────────────────────────────────
#
# C^{n+1} = C_m + Σ_i f_i (C_i − C_m) A_{εε,i}^{(C^n)}
#
# Sum runs over INCLUSIONS only (matrix excluded, treated implicitly via
# C_m). This is the C++ reference's `evaluate_dilute(X)` body.

function _asc_iterate_stiffness(rve::RVE, asc::AsymmetricSelfConsistent,
                                 ::Val{p}; kw...) where {p}
    C_m  = matrix_property(rve, p)
    solver_kw, step_kw = _split_sc_kwargs(kw)
    step = C_n -> _asc_step_stiffness(rve, p, C_m, C_n; step_kw...)
    return _solve_sc(asc.algorithm, step, C_m; asc.options..., solver_kw...)
end

function _asc_step_stiffness(rve::RVE, prop::Symbol, C_m, C_n; kw...)
    return _asc_step_stiffness_dispatch(rve, prop, C_m, C_n; kw...)
end

function _asc_step_stiffness_dispatch(rve::RVE, prop::Symbol,
                                       C_m::TensND.AbstractTens{4, 3},
                                       C_n::TensND.AbstractTens{4, 3}; kw...)
    Δ = zero(C_n)
    for name in inclusion_phase_names(rve)
        a = rve.amounts[name]
        a isa VolumeFraction || continue   # cracks not supported in stiffness form
        f = amount_value(a)
        C_i = phase_property(rve, name, prop)
        A_dil = _phase_dilute_concentration(rve, name, prop, C_n; kw...)
        Δ += f * ((C_i - C_m) ⊡ A_dil)
    end
    return C_m + Δ
end

function _asc_step_stiffness_dispatch(rve::RVE, prop::Symbol,
                                       K_m::TensND.AbstractTens{2, 3},
                                       K_n::TensND.AbstractTens{2, 3}; kw...)
    Δ = zero(K_n)
    for name in inclusion_phase_names(rve)
        a = rve.amounts[name]
        a isa VolumeFraction || continue
        f = amount_value(a)
        K_i = phase_property(rve, name, prop)
        A_dil = _phase_dilute_concentration(rve, name, prop, K_n; kw...)
        Δ += f * ((K_i - K_m) ⋅ A_dil)
    end
    return K_m + Δ
end

# ── Compliance-form iteration ───────────────────────────────────────────────
#
# Iterating estimate is the *stiffness* `C^n` (we invert at the end);
# the formula reads
#
#   S^{n+1} = S_m + ⟨A_{εε}^{(C^n)}⟩ S^n − S_m ⟨C_i A_{εε}^{(C^n)}⟩ S^n
#
# where S_n = inv(C^n), S_m = inv(C_m), and the sum is over INCLUSIONS
# only. We then invert `S^{n+1}` to recover `C^{n+1}` so the iteration
# stays in stiffness space and the same `_solve_sc` driver is reused.

function _asc_iterate_compliance(rve::RVE, asc::AsymmetricSelfConsistent,
                                  ::Val{p}; kw...) where {p}
    C_m = matrix_property(rve, p)
    S_m = inv(C_m)
    solver_kw, step_kw = _split_sc_kwargs(kw)
    step = C_n -> _asc_step_compliance(rve, p, C_m, S_m, C_n; step_kw...)
    return _solve_sc(asc.algorithm, step, C_m; asc.options..., solver_kw...)
end

function _asc_step_compliance(rve::RVE, prop::Symbol, C_m, S_m, C_n; kw...)
    return _asc_step_compliance_dispatch(rve, prop, C_m, S_m, C_n; kw...)
end

function _asc_step_compliance_dispatch(rve::RVE, prop::Symbol,
                                        C_m::TensND.AbstractTens{4, 3},
                                        S_m::TensND.AbstractTens{4, 3},
                                        C_n::TensND.AbstractTens{4, 3}; kw...)
    S_n = inv(C_n)
    A_avg  = zero(C_n)
    CA_avg = zero(C_n)
    for name in inclusion_phase_names(rve)
        a = rve.amounts[name]
        # Volume-fraction inclusions: standard A_εε.
        # Crack densities: the strain concentration is singular; their
        # "compliance contribution" is the natural object instead and is
        # added directly to S^{n+1} via `compliance_contribution`.
        if a isa VolumeFraction
            f = amount_value(a)
            C_i = phase_property(rve, name, prop)
            A_dil = _phase_dilute_concentration(rve, name, prop, C_n; kw...)
            A_avg  += f * A_dil
            CA_avg += f * (C_i ⊡ A_dil)
        end
    end
    S_new = S_m + (A_avg - S_m ⊡ CA_avg) ⊡ S_n
    # Crack contributions add directly to the compliance.
    for name in inclusion_phase_names(rve)
        a = rve.amounts[name]
        a isa VolumeFraction && continue
        H = _phase_compliance_contribution(rve, name, prop, C_n; kw...)
        S_new += H
    end
    return inv(S_new)
end

function _asc_step_compliance_dispatch(rve::RVE, prop::Symbol,
                                        K_m::TensND.AbstractTens{2, 3},
                                        R_m::TensND.AbstractTens{2, 3},
                                        K_n::TensND.AbstractTens{2, 3}; kw...)
    R_n = inv(K_n)
    A_avg  = zero(K_n)
    KA_avg = zero(K_n)
    for name in inclusion_phase_names(rve)
        a = rve.amounts[name]
        a isa VolumeFraction || continue
        f = amount_value(a)
        K_i = phase_property(rve, name, prop)
        A_dil = _phase_dilute_concentration(rve, name, prop, K_n; kw...)
        A_avg  += f * A_dil
        KA_avg += f * (K_i ⋅ A_dil)
    end
    R_new = R_m + (A_avg - R_m ⋅ KA_avg) ⋅ R_n
    return inv(R_new)
end

# ── Legacy compliance-space RVE builder (kept for reference) ────────────────
function _rve_in_compliance_space(rve::RVE{T, S}, prop::Symbol) where {T, S}
    rve_S = RVE(rve.matrix_name; T = T,
                distribution_shape = rve.distribution_shape)
    m_phase = matrix_phase(rve)
    add_matrix!(rve_S, m_phase.geometry, Dict(:S => inv(m_phase.properties[prop])))
    for name in inclusion_phase_names(rve)
        ph = rve.phases[name]
        a  = rve.amounts[name]
        new_props = Dict(:S => inv(ph.properties[prop]))
        if a isa VolumeFraction
            add_phase!(rve_S, name, ph.geometry, new_props;
                       fraction = amount_value(a))
        else
            add_phase!(rve_S, name, ph.geometry, new_props;
                       density = amount_value(a))
        end
    end
    return rve_S
end
