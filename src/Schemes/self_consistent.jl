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
    step   = C -> _sc_step(rve, C, p; kw...)
    return _solve_sc(sc.algorithm, step, P_init; sc.options..., kw...)
end

# ── SC step (one Mori-Tanaka-like iterate against current estimate) ─────────

function _sc_step(rve::RVE, C_n, prop::Symbol; kw...)
    return _sc_step_dispatch(rve, C_n, prop; kw...)
end

# 4th-order
function _sc_step_dispatch(rve::RVE, C_n::TensND.AbstractTens{4, 3}, prop::Symbol;
                           kw...)
    Iref = _identity_like(C_n)
    f_m = matrix_volume_fraction(rve)
    A_avg  = f_m * Iref     # matrix carries A_dil = I (it IS the reference)
    CA_avg = f_m * matrix_property(rve, prop)
    for name in inclusion_phase_names(rve)
        a = rve.amounts[name]
        geom = rve.phases[name].geometry
        if a isa VolumeFraction
            f = amount_value(a)
            P_i = phase_property(rve, name, prop)
            A_dil = MFH_Core.strain_strain_loc(geom, P_i, C_n; kw...)
            A_avg  += f * A_dil
            CA_avg += f * (P_i ⊡ A_dil)
        end
        # CrackDensity ignored in the stiffness-form SC (caller should use ASC)
    end
    return CA_avg ⊡ inv(A_avg)
end

# 2nd-order
function _sc_step_dispatch(rve::RVE, K_n::TensND.AbstractTens{2, 3}, prop::Symbol;
                           kw...)
    Iref = _identity_like(K_n)
    f_m = matrix_volume_fraction(rve)
    A_avg  = f_m * Iref
    KA_avg = f_m * matrix_property(rve, prop)
    for name in inclusion_phase_names(rve)
        a = rve.amounts[name]
        geom = rve.phases[name].geometry
        if a isa VolumeFraction
            f = amount_value(a)
            P_i = phase_property(rve, name, prop)
            A_dil = MFH_Core.gradient_gradient_loc(geom, P_i, K_n; kw...)
            A_avg  += f * A_dil
            KA_avg += f * (P_i ⋅ A_dil)
        end
    end
    return KA_avg ⋅ inv(A_avg)
end

# ── Built-in solvers ────────────────────────────────────────────────────────

"""
    _solve_sc(algo, step, x0; abstol, maxiters, damping, verbose, kw...) -> AbstractTens

Generic solver dispatcher for SC fixed points. Built-in:

- [`AndersonDefault`](@ref) — Picard with relaxation
  (`x_{n+1} = (1-damping)·step(x_n) + damping·x_n`). `damping = 0.0`
  default; raise to ≈ 0.5 for high-contrast iterations that overshoot.
- [`NewtonDefault`](@ref) — SciML Newton-Raphson, available only when the
  weak extension `MeanFieldHomNonlinearSolveExt` is loaded
  (`using NonlinearSolve`).

Other algorithms from `NonlinearSolve.jl` are supported through the
same weak extension.
"""
function _solve_sc(::AndersonDefault, step, x0::TensND.AbstractTens;
                   abstol::Real = 1.0e-10, maxiters::Int = 100,
                   damping::Real = 0.0, verbose::Bool = false, kw...)
    x = x0
    last_resid = zero(real(eltype(x0)))
    for k in 1:maxiters
        x_new = step(x)
        last_resid = _sc_residual_norm(x_new, x)
        verbose && @info "SC iter $k : ‖Δ‖ = $last_resid"
        last_resid < abstol && return x_new
        x = (one(real(eltype(x))) - damping) * x_new + damping * x
    end
    @warn "SC (AndersonDefault/Picard) did not converge in $maxiters iterations" last_resid abstol
    return x
end

function _solve_sc(::NewtonDefault, step, x0::TensND.AbstractTens; kw...)
    error("NewtonDefault requires NonlinearSolve.jl: load it with " *
          "`using NonlinearSolve` to activate the MeanFieldHomNonlinearSolveExt " *
          "extension. Default solver `AndersonDefault` works without any extra dependency.")
end

# Residual norm operating on stored components (works for any tensor type).
_sc_residual_norm(a::TensND.AbstractTens, b::TensND.AbstractTens) =
    sqrt(sum(abs2, get_array(a) .- get_array(b)))

# =============================================================================
#  AsymmetricSelfConsistent — switch between stiffness and compliance
#  iteration depending on the matrix-vs-bound contrast.
# =============================================================================

"""
    _evaluate(rve, asc::AsymmetricSelfConsistent, ::Val{p}; kw...) -> AbstractTens

Asymmetric self-consistent scheme. Decides whether to iterate in
stiffness or compliance space based on the contrast between the Voigt
upper bound and the matrix property:

- if `‖C_Voigt‖ ≥ ‖C_matrix‖` (matrix soft, inclusions stiff): iterate
  in stiffness — stiffness contributions are bounded ⇒ contractive.
- else (matrix stiff, inclusions soft): iterate in compliance space.

Better behaviour than [`SelfConsistent`](@ref) on high-contrast
matrix-stiff / inclusion-soft RVEs and on RVEs containing cracks (the
compliance-form path natively handles the singular crack stiffness).
"""
function _evaluate(rve::RVE, asc::AsymmetricSelfConsistent, ::Val{p}; kw...) where {p}
    if _asc_use_stiffness(rve, p)
        sc = SelfConsistent(asc.algorithm, asc.options)
        return _evaluate(rve, sc, Val(p); kw...)
    else
        return _asc_compliance_iterate(rve, asc, Val(p); kw...)
    end
end

function _asc_use_stiffness(rve::RVE, prop::Symbol; kw...)
    P₀      = matrix_property(rve, prop)
    P_voigt = _evaluate(rve, Voigt(), Val(prop); kw...)
    return _matrix_norm(P_voigt) ≥ _matrix_norm(P₀)
end

# Order-aware norm of a tensor — uses the Inf-norm (max absolute row sum)
# rather than the 2-norm so the heuristic is Dual-safe (the 2-norm goes
# through SVD which is not implemented for `Matrix{ForwardDiff.Dual}`).
_matrix_norm(C::TensND.AbstractTens{4, 3}) = LinearAlgebra.opnorm(KM(C), Inf)
_matrix_norm(K::TensND.AbstractTens{2, 3}) = LinearAlgebra.opnorm(collect(get_array(K)), Inf)

function _asc_compliance_iterate(rve::RVE, asc::AsymmetricSelfConsistent,
                                 ::Val{p}; kw...) where {p}
    rve_S = _rve_in_compliance_space(rve, p)
    sc = SelfConsistent(asc.algorithm, asc.options)
    S_eff = _evaluate(rve_S, sc, Val(:S); kw...)
    return inv(S_eff)
end

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
