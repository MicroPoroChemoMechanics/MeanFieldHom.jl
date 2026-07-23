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
    # Symmetric Hill / Budiansky SC iteration on the ECHOES `B · A^{-1}`
    # body : crack phases contribute their compliance contribution
    # `H_c(C_n)` to the denominator A_avg via `_phase_compliance_contribution`,
    # and their stiffness contribution `ΔC_c(C_n)` to the numerator
    # CA_avg (traction-free → no stress contribution from solid side).
    # The eigenvalue guard `_sc_pd_guard` mirrors ECHOES
    # `homogenization_scheme.h::evaluate` and prevents the iteration
    # from collapsing to the trivial percolated fixed point at moderate
    # density.
    P_init = matrix_property(rve, p)
    solver_kw, step_kw = _split_sc_kwargs(kw)
    step = C -> _sc_step(rve, C, p; step_kw...)
    return _solve_sc(sc.algorithm, step, P_init; sc.options..., solver_kw...)
end

# Solver-only kwargs are intercepted at this level and never forwarded
# to the underlying `_sc_step` (which would otherwise leak them down to
# `hill_tensor` and trigger a `MethodError` on the unknown kwarg).
const _SC_SOLVER_KWARGS = (
    :abstol, :reltol, :maxiters, :damping,
    :verbose, :select_best,
)

function _split_sc_kwargs(kw)
    solver = Dict{Symbol, Any}()
    step = Dict{Symbol, Any}()
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
function _sc_step_dispatch(
        rve::RVE, C_n::TensND.AbstractTens{4, 3}, prop::Symbol;
        kw...
    )
    # ECHOES `homogenization_scheme.h::evaluate` body :
    #   strain_Stress_α  = A_α(C_n) · S_n   (solid) [`compute_strain_Stress`]
    #   strain_Stress_c  = sym(H_c(C_n))    (void crack — NO trailing S_n!)
    #   stress_Stress_α  = C_α · strain_Stress_α
    #   stress_Stress_c  = 0                 (traction-free)
    # Accumulators :  A_E = Σ f_α·sym(strain_Stress_α)
    #                 B_E = Σ f_α·sym(stress_Stress_α)
    # Result  : C_eff = sym(B_E · A_E^{-vol}).
    # The trailing `S_n` factor cancels between A_E and B_E for solid
    # phases — but NOT for cracks, whose `strain_Stress` is the bare
    # compliance contribution `H_c`.  This breaks the cancellation and
    # gives a different fixed point than the textbook
    # `(Σ f·C·A)·(Σ f·A)^{-1}` SC body when cracks are present.
    A_avg = zero(C_n)   # = Σ_solids f·sym(A_α)
    CA_avg = zero(C_n)   # = Σ_solids f·sym(C_α·A_α)
    H_total = zero(C_n)   # = Σ_cracks ε·sym(H_c)
    has_cracks = false
    for name in rve.phase_names
        if name === rve.matrix_name
            f = matrix_volume_fraction(rve)
        else
            a = rve.amounts[name]
            a isa VolumeFraction || continue
            f = amount_value(a)
        end
        P_α = phase_property(rve, name, prop)
        A_dil = _phase_dilute_concentration(rve, name, prop, C_n; kw...)
        A_avg += f * A_dil
        CA_avg += f * (P_α ⊡ A_dil)
    end
    for name in inclusion_phase_names(rve)
        a = rve.amounts[name]
        a isa CrackDensity || continue
        H = _phase_compliance_contribution(rve, name, prop, C_n; kw...)
        H_total += H
        has_cracks = true
    end
    if has_cracks
        S_n = inv(C_n)
        A_E = (A_avg ⊡ S_n) + H_total
        B_E = CA_avg ⊡ S_n
        return B_E ⊡ inv(A_E)
    else
        return CA_avg ⊡ inv(A_avg)
    end
end

# 2nd-order — same symmetric SC structure for conductivity / diffusion.
function _sc_step_dispatch(
        rve::RVE, K_n::TensND.AbstractTens{2, 3}, prop::Symbol;
        kw...
    )
    # Conduction analogue of the 4th-order ECHOES SC body :
    # solids have `gradient_Flux = A · R_n` (R_n = inv(K_n) — resistivity),
    # cracks contribute the bare resistivity contribution `R_c`.
    A_avg = zero(K_n)
    KA_avg = zero(K_n)
    R_total = zero(K_n)
    has_cracks = false
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
        A_avg += f * A_dil
        KA_avg += f * (P_α ⋅ A_dil)
    end
    for name in inclusion_phase_names(rve)
        a = rve.amounts[name]
        a isa CrackDensity || continue
        R = _phase_compliance_contribution(rve, name, prop, K_n; kw...)
        R_total += R
        has_cracks = true
    end
    if has_cracks
        R_n = inv(K_n)
        A_E = (A_avg ⋅ R_n) + R_total
        B_E = KA_avg ⋅ R_n
        return B_E ⋅ inv(A_E)
    else
        return KA_avg ⋅ inv(A_avg)
    end
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
function _solve_sc(
        ::AndersonDefault, step, x0::TensND.AbstractTens;
        abstol::Real = 1.0e-12, reltol::Real = 1.0e-8,
        maxiters::Int = 100,
        damping::Real = 0.0, verbose::Bool = false,
        select_best::Bool = false,
        kw...
    )
    x = x0
    last_resid = _sc_residual_zero(x0)
    x_best = x0
    resid_best_val = typemax(_value_eltype(x0))
    ε_pos = _sc_pd_eps(x0)
    for k in 1:maxiters
        x = _sc_pd_guard_apply(x, ε_pos)
        x_new = step(x)
        last_resid = _sc_residual_norm(x_new, x)
        norm_x = _sc_residual_norm(x, zero(x))
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

"""
    _solve_sc(::NewtonDefault, step, x0::AbstractTens; …) -> AbstractTens

Built-in Newton-Raphson SC solver, parameterising the iterating
estimate by its symmetry-class **canonical components**
(`TensND.get_data` → `(α, β)` for iso, `(ℓ₁, …, ℓ₆)` for TI / Walpole,
9 components for ortho).  At each Newton step:

1. Build the residual `F(p) = canonical(step(rebuild(p))) − p`,
2. Compute the Jacobian `J = ∂F/∂p` via `ForwardDiff.jacobian`,
3. Take the Newton step `Δp = −J⁻¹·F(p)` with backtracking line
   search (Armijo with shrinking factor 1/2, minimum step 1e-6).
4. Fall back to a single Picard step when the line search fails.

Compared to the SciML weak-extension path, this is dependency-free and
specialised to the small parameter spaces of `MeanFieldHom` symmetry
classes (≤ 21 components for the most general aniso 4-tensor); the
Jacobian is computed once per iteration through the same `step`
function the AndersonDefault loop calls.
"""
function _solve_sc(
        ::NewtonDefault, step, x0::TensND.AbstractTens;
        abstol::Real = 1.0e-12, reltol::Real = 1.0e-8,
        maxiters::Int = 50,
        damping::Real = 0.0, verbose::Bool = false,
        select_best::Bool = false,
        kw...
    )
    p0 = _tens_to_param_vec(x0)
    L = length(p0)
    Tref = float(eltype(p0))
    rebuild = p -> _tens_from_param_vec(x0, p)
    ε_pos = _sc_pd_eps(x0)
    residual_vec = function (p)
        x_in = _sc_pd_guard_apply(rebuild(p), ε_pos)
        x_out = step(x_in)
        return _tens_to_param_vec(x_out) .- _tens_to_param_vec(x_in)
    end
    p = collect(Tref, p0)
    p_best = copy(p); resid_best = Inf
    for iter in 1:maxiters
        r = residual_vec(p)
        norm_r = sqrt(sum(abs2, r))
        norm_p = sqrt(sum(abs2, p))
        tol_eff = abstol + reltol * norm_p
        verbose && @info "SC-Newton iter $iter : ‖F‖ = $norm_r   tol = $tol_eff"
        if select_best && norm_r < resid_best
            resid_best = norm_r
            p_best .= p
        end
        if norm_r ≤ tol_eff
            return rebuild(p)
        end
        # Jacobian via ForwardDiff (strong dependency).
        J = ForwardDiff.jacobian(residual_vec, p)
        # Solve J·δ = -r with a regularising fallback if J is singular.
        δ = try
            J \ (-r)
        catch err
            @debug "SC-Newton: linear solve failed, applying tiny Tikhonov" err
            (J + 1.0e-10 * sqrt(sum(abs2, J)) * LinearAlgebra.I) \ (-r)
        end
        # Backtracking line search (Armijo).
        α_step = one(Tref)
        accepted = false
        for _ in 1:30
            p_new = p .+ α_step .* δ
            r_new = residual_vec(p_new)
            if sqrt(sum(abs2, r_new)) ≤ (1 - 1.0e-4 * α_step) * norm_r
                p .= p_new
                accepted = true
                break
            end
            α_step /= 2
            α_step < 1.0e-6 && break
        end
        if !accepted
            # Fall back to a damped Picard step.
            verbose && @info "SC-Newton: line search failed, taking Picard step"
            p .= residual_vec(p) .+ p   # one Picard step
        end
    end
    @debug "SC-Newton: maxiters reached without convergence" maxiters
    return rebuild(select_best ? p_best : p)
end

# ── Positive-definite guard for the SC running estimate ───────────────────
#
# At high contrast (porous SC near the percolation threshold, cracks at
# moderate density), the SC iteration map can have a stable fixed point
# at the trivial null tensor (`C = 0`).  A Picard iteration starting
# from `C_M` may drift into this percolation fixed point even when a
# physically meaningful finite fixed point exists nearby.  The ECHOES
# `homogenization_scheme.h::evaluate` mitigates this by detecting a
# negative-definite running estimate and resetting it to a tiny
# positive identity (`mX = EPSILON · I`) before each step — this
# prevents the iteration from collapsing to zero and lets it find the
# physical finite-modulus branch when it exists.  We mirror the same
# guard here: when any canonical eigenvalue of the running estimate
# falls below a relative threshold, we reset it to `ε · α_M_init`
# (matrix scale).

function _sc_pd_guard(x, x0)
    return _sc_pd_guard_apply(x, _sc_pd_eps(x0))
end

# `x0` is fixed for the entire SC iteration (Picard or Newton) — computing
# this once and passing `ε_pos` into `_sc_pd_guard_apply` directly avoids
# recomputing `_max_canonical_value(x0)` (which itself does a `try`/`catch`,
# blocking inlining) on every single guard call within the loop.
function _sc_pd_eps(x0)
    α0_max = _max_canonical_value(x0)
    return max(α0_max * sqrt(eps(real(_value_eltype(x0)))), 1.0e-12)
end

# Iso 4-tensor: check (α, β); reset to ε·𝕁 + ε·𝕂 if either component
# is non-positive.
function _sc_pd_guard_apply(x::TensND.TensISO{3}, ε_pos)
    α, β = TensND.get_data(x)
    if real(α) ≤ ε_pos || real(β) ≤ ε_pos
        return TensND.TensISO{3}(max(real(α), ε_pos), max(real(β), ε_pos))
    end
    return x
end

# Default: try to detect the smallest canonical component value.  If
# any canonical component is non-positive, reset to a tiny positive
# baseline of the same shape.  Falls back to passthrough for tensors
# whose `get_data` does not give a meaningful "smallest eigenvalue"
# notion.
function _sc_pd_guard_apply(x::TensND.AbstractTens, ε_pos)
    try
        d = TensND.get_data(x)
        for v in d
            real(v) ≤ ε_pos && return _rebuild_min_eps(x, ε_pos)
        end
    catch
    end
    return x
end

function _rebuild_min_eps(x::TensND.AbstractTens, ε_pos)
    d = TensND.get_data(x)
    new_d = Tuple(real(v) ≤ ε_pos ? ε_pos : v for v in d)
    return TensND._rebuild(x, new_d)
end

# Used to scale the eigenvalue tolerance.
_max_canonical_value(x::TensND.AbstractTens) = begin
    try
        d = TensND.get_data(x)
        return maximum(abs(real(v)) for v in d)
    catch
        return 1.0
    end
end

# ── Tens ↔ canonical parameter vector helpers ──────────────────────────────
#
# `TensND.get_data(t)` already returns the canonical tuple of an iso /
# TI / ortho / Walpole tensor (`(α, β)` for iso, `(ℓ₁, …, ℓ₆)` for
# Walpole TI, etc.).  We `collect` it into a `Vector` for the
# ForwardDiff-friendly residual function, and rebuild the same
# concrete type via the canonical constructor inferred from the
# prototype.
_tens_to_param_vec(t::TensND.AbstractTens) = collect(TensND.get_data(t))
_tens_from_param_vec(prototype::TensND.AbstractTens, p::AbstractVector) =
    _rebuild_from_data(prototype, p)

# Default rebuild: `_rebuild` is TensND-internal; for known types we use
# the public constructor.
_rebuild_from_data(::TensND.TensISO{3}, p) = TensND.TensISO{3}(p[1], p[2])
_rebuild_from_data(t::TensND.AbstractTens, p) =
    TensND._rebuild(t, ntuple(i -> p[i], length(p)))

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
# Crack phases (`CrackDensity`) carry no volume but still soften the
# composite ; the Voigt heuristic above ignores them and may wrongly
# select the stiffness form (which skips cracks → no degradation).
# Force the compliance form whenever an RVE has at least one crack.
function _asc_use_stiffness(rve::RVE, prop::Symbol; kw...)
    if any(a isa CrackDensity for a in values(rve.amounts))
        return false
    end
    P₀ = matrix_property(rve, prop)
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

function _asc_iterate_stiffness(
        rve::RVE, asc::AsymmetricSelfConsistent,
        ::Val{p}; kw...
    ) where {p}
    C_m = matrix_property(rve, p)
    solver_kw, step_kw = _split_sc_kwargs(kw)
    step = C_n -> _asc_step_stiffness(rve, p, C_m, C_n; step_kw...)
    return _solve_sc(asc.algorithm, step, C_m; asc.options..., solver_kw...)
end

function _asc_step_stiffness(rve::RVE, prop::Symbol, C_m, C_n; kw...)
    return _asc_step_stiffness_dispatch(rve, prop, C_m, C_n; kw...)
end

function _asc_step_stiffness_dispatch(
        rve::RVE, prop::Symbol,
        C_m::TensND.AbstractTens{4, 3},
        C_n::TensND.AbstractTens{4, 3}; kw...
    )
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

function _asc_step_stiffness_dispatch(
        rve::RVE, prop::Symbol,
        K_m::TensND.AbstractTens{2, 3},
        K_n::TensND.AbstractTens{2, 3}; kw...
    )
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

function _asc_iterate_compliance(
        rve::RVE, asc::AsymmetricSelfConsistent,
        ::Val{p}; kw...
    ) where {p}
    C_m = matrix_property(rve, p)
    S_m = inv(C_m)
    solver_kw, step_kw = _split_sc_kwargs(kw)
    step = C_n -> _asc_step_compliance(rve, p, C_m, S_m, C_n; step_kw...)
    return _solve_sc(asc.algorithm, step, C_m; asc.options..., solver_kw...)
end

function _asc_step_compliance(rve::RVE, prop::Symbol, C_m, S_m, C_n; kw...)
    return _asc_step_compliance_dispatch(rve, prop, C_m, S_m, C_n; kw...)
end

function _asc_step_compliance_dispatch(
        rve::RVE, prop::Symbol,
        C_m::TensND.AbstractTens{4, 3},
        S_m::TensND.AbstractTens{4, 3},
        C_n::TensND.AbstractTens{4, 3}; kw...
    )
    S_n = inv(C_n)
    A_avg = zero(C_n)
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
            A_avg += f * A_dil
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

function _asc_step_compliance_dispatch(
        rve::RVE, prop::Symbol,
        K_m::TensND.AbstractTens{2, 3},
        R_m::TensND.AbstractTens{2, 3},
        K_n::TensND.AbstractTens{2, 3}; kw...
    )
    R_n = inv(K_n)
    A_avg = zero(K_n)
    KA_avg = zero(K_n)
    for name in inclusion_phase_names(rve)
        a = rve.amounts[name]
        a isa VolumeFraction || continue
        f = amount_value(a)
        K_i = phase_property(rve, name, prop)
        A_dil = _phase_dilute_concentration(rve, name, prop, K_n; kw...)
        A_avg += f * A_dil
        KA_avg += f * (K_i ⋅ A_dil)
    end
    R_new = R_m + (A_avg - R_m ⋅ KA_avg) ⋅ R_n
    return inv(R_new)
end

# ── Legacy compliance-space RVE builder (kept for reference) ────────────────
function _rve_in_compliance_space(rve::RVE{T, S}, prop::Symbol) where {T, S}
    rve_S = RVE(
        rve.matrix_name; T = T,
        distribution_shape = rve.distribution_shape
    )
    m_phase = matrix_phase(rve)
    add_matrix!(rve_S, m_phase.geometry, Dict(:S => inv(m_phase.properties[prop])))
    for name in inclusion_phase_names(rve)
        ph = rve.phases[name]
        a = rve.amounts[name]
        new_props = Dict(:S => inv(ph.properties[prop]))
        if a isa VolumeFraction
            add_phase!(
                rve_S, name, ph.geometry, new_props;
                fraction = amount_value(a)
            )
        else
            add_phase!(
                rve_S, name, ph.geometry, new_props;
                density = amount_value(a)
            )
        end
    end
    return rve_S
end
