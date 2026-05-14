# =============================================================================
#  differential.jl — Differential homogenisation scheme as a SciML ODE.
#
#  Integrates the multi-phase Norris ODE on the fictitious incorporation
#  time `τ ∈ [0, 1]` (cf. user's hand-written DEM note ; @norris1985) :
#
#      dC^hom / dτ = Σ_α dφ_α/dτ · (C_α - C^hom) ⊡ A_α^dil(C^hom)
#                  + Σ_c dε_c/dτ · ΔC^crack_c(C^hom)
#
#  with the volumetric balance `df = (𝟙 − f ⊗ 𝐔) · dφ` inverted by
#  Sherman-Morrison so the user supplies effective volume fractions
#  `f_α(τ)` along the chosen `trajectory`.  The ODE is solved by
#  `OrdinaryDiffEq.solve` with an adaptive RK (default `Tsit5`) and
#  the result returned at `τ = 1`.
#
#  Cracks have negligible volume — their density `ε_c(τ)` enters the
#  RHS directly without going through Sherman-Morrison.
# =============================================================================

"""
    _evaluate(rve, scheme::DifferentialScheme, ::Val{p}; kw...) -> AbstractTens

Differential scheme for property `:p` ([Norris 1985](@cite norris1985)).
Integrates the multi-phase incorporation-sequence ODE on `τ ∈ [0, 1]`
with the SciML `OrdinaryDiffEq.solve` driver (default `Tsit5`).
"""
function _evaluate(rve::RVE, scheme::DifferentialScheme, ::Val{p}; kw...) where {p}
    nsteps = get(scheme.options, :nsteps, 100)
    abstol = get(scheme.options, :abstol, 1.0e-8)
    reltol = get(scheme.options, :reltol, 1.0e-6)
    alg = get(scheme.options, :alg, nothing)
    paths = _resolve_paths(scheme.trajectory, rve, nsteps)
    P_init = matrix_property(rve, p)
    return _diff_integrate_ode(
        rve, paths, p, P_init;
        nsteps, abstol, reltol, alg, kw...
    )
end

# ── ODE integrator ──────────────────────────────────────────────────────────

function _diff_integrate_ode(
        rve::RVE,
        paths::AbstractDict{Symbol},
        prop::Symbol,
        P_init::TensND.AbstractTens;
        nsteps::Int = 100,
        abstol::Real = 1.0e-8,
        reltol::Real = 1.0e-6,
        alg = nothing,
        kw...
    )
    # Split inclusion phases between solids and cracks.  Targets are the
    # final values reached at `τ = 1` (volume fractions for solids,
    # densities for cracks).  The target eltype is preserved (it can be
    # `ForwardDiff.Dual` when the user differentiates `homogenize`
    # through its scalar inputs).
    solid_names = Symbol[]
    crack_names = Symbol[]
    incl = inclusion_phase_names(rve)
    T_target = isempty(incl) ? Float64 :
        promote_type((typeof(amount_value(rve.amounts[name])) for name in incl)...)
    targets = Dict{Symbol, T_target}()
    for name in incl
        a = rve.amounts[name]
        if a isa VolumeFraction
            push!(solid_names, name)
        else  # CrackDensity
            push!(crack_names, name)
        end
        targets[name] = T_target(amount_value(a))
    end

    # State : canonical components of `P_init` (length 2 / 5 / 9 for
    # iso/TI/ortho 4-tensors, 1 / 3 / 6 for the corresponding 2-tensors,
    # 36 / 9 for the fully-anisotropic Mandel fallback).  When any
    # phase carries a contribution that may break the matrix's
    # symmetry class (an inclusion / crack with `symmetrize=:none` is
    # the typical case), we fall back to the full Mandel
    # representation so the ODE state can absorb it.  For ForwardDiff
    # sensitivity we promote the state eltype to whichever type
    # accommodates both the matrix property and the per-phase targets.
    sym_tag = _diff_state_tag(rve, P_init)
    x0 = _get_state(sym_tag, P_init)
    T_state = promote_type(eltype(x0), T_target)
    x0 = T_state.(x0)
    ode_kw = (
        rve = rve,
        prop = prop,
        paths = paths,
        solid_names = solid_names,
        crack_names = crack_names,
        targets = targets,
        sym_tag = sym_tag,
        proto = P_init,
        kw = kw,
    )
    rhs! = (du, u, p, τ) -> _diff_ode_rhs!(du, u, p, τ)
    prob = ODEProblem(rhs!, x0, (0.0, 1.0), ode_kw)
    sol = solve(
        prob,
        alg === nothing ? Tsit5() : alg;
        abstol, reltol,
        saveat = range(0.0, 1.0; length = max(nsteps, 1) + 1),
        dense = false
    )
    return _reconstruct_tens(sym_tag, P_init, sol.u[end])
end

# ── RHS ─────────────────────────────────────────────────────────────────────

function _diff_ode_rhs!(du, u, p, τ)
    Cτ = _reconstruct_tens(p.sym_tag, p.proto, u)
    # Sherman-Morrison inversion of dφ = (𝟙 − f ⊗ 𝐔)^{-1} · df.
    n_solid = length(p.solid_names)
    f = Vector{eltype(u)}(undef, n_solid)
    df = Vector{eltype(u)}(undef, n_solid)
    @inbounds for (i, name) in enumerate(p.solid_names)
        nt = p.paths[name]
        f[i] = nt.f(τ) * p.targets[name]
        df[i] = nt.df(τ) * p.targets[name]
    end
    f0 = one(eltype(u)) - sum(f; init = zero(eltype(u)))
    sum_df = sum(df; init = zero(eltype(u)))
    Δ = zero(Cτ)
    @inbounds for (i, name) in enumerate(p.solid_names)
        # dφ_i / dτ = df_i + (f_i / f_0) · sum(df)   (Sherman-Morrison)
        dφᵢ = df[i] + (f[i] / f0) * sum_df
        iszero(dφᵢ) && continue
        Δ += dφᵢ * _diff_dilute_correction(p.rve, name, p.prop, Cτ; p.kw...)
    end
    @inbounds for name in p.crack_names
        nt = p.paths[name]
        dε = nt.df(τ) * p.targets[name]
        iszero(dε) && continue
        Δ += dε * _diff_crack_density_kernel(p.rve, name, p.prop, Cτ; p.kw...)
    end
    _set_state!(du, p.sym_tag, Δ)
    return nothing
end

# ── Symmetry tag + reconstruction helpers ──────────────────────────────────

# We use a small sum type stored in `p.sym_tag` to dispatch the
# reconstruction (and the initial flatten) without a dynamic
# constructor lookup at every RHS step.
# Decide which state representation the ODE solver should use, given
# the matrix property and the per-phase symmetrize keywords.  An
# isotropic running estimate may evolve out of the iso class when a
# crack carries `NoSymmetrize` (the crack's TI contribution in its
# canonical frame leaks anisotropy into the running matrix at every
# RHS step).  Solid VolumeFraction phases remain iso-compatible when
# both their property and the matrix are iso, regardless of the
# `symmetrize` flag.  In the anisotropy-leaking case we fall back to
# the full Mandel 6×6 representation (`:full_4`).
function _diff_state_tag(rve::RVE, P_init::TensND.AbstractTens)
    base = _symmetry_tag(P_init)
    base === Val(:iso_4) || return base
    for name in inclusion_phase_names(rve)
        a = rve.amounts[name]
        a isa CrackDensity || continue
        sym = phase_symmetrize(rve, name)
        sym isa NoSymmetrize && return Val(:full_4)
    end
    return base
end

_symmetry_tag(::TensND.TensISO{4}) = Val(:iso_4)
_symmetry_tag(::TensND.TensISO{2}) = Val(:iso_2)
_symmetry_tag(::TensND.TensTI) = Val(:ti)
_symmetry_tag(::TensND.TensOrtho) = Val(:ortho)
# Fully-anisotropic 4-tensor (`TensCanonical`, `Tens`, …) — flatten
# the full 6×6 Mandel matrix as the ODE state (36 components ; not
# minimal but generic).
_symmetry_tag(::TensND.AbstractTens{4, 3}) = Val(:full_4)
_symmetry_tag(::TensND.AbstractTens{2, 3}) = Val(:full_2)

# Flatten / reconstruct helpers : the iso / TI / ortho cases use the
# canonical components ; the full-aniso fallback uses the Mandel
# array.

# Initial state vector for the ODE solver.
_get_state(::Val{:iso_4}, t) = collect(TensND.get_data(t))
_get_state(::Val{:iso_2}, t) = collect(TensND.get_data(t))
_get_state(::Val{:ti}, t) = collect(TensND.get_data(t))
_get_state(::Val{:ortho}, t) = collect(TensND.get_data(t))
_get_state(::Val{:full_4}, t) = vec(collect(TensND.KM(t)))
_get_state(::Val{:full_2}, t) = vec(collect(TensND.KM(t)))

# Push a tensor into a flat vector for `du`.
_set_state!(du, ::Val{:iso_4}, Δ) = (du .= TensND.get_data(Δ); nothing)
_set_state!(du, ::Val{:iso_2}, Δ) = (du .= TensND.get_data(Δ); nothing)
_set_state!(du, ::Val{:ti}, Δ) = (du .= TensND.get_data(Δ); nothing)
_set_state!(du, ::Val{:ortho}, Δ) = (du .= TensND.get_data(Δ); nothing)
_set_state!(du, ::Val{:full_4}, Δ) = (du .= vec(TensND.KM(Δ)); nothing)
_set_state!(du, ::Val{:full_2}, Δ) = (du .= vec(TensND.KM(Δ)); nothing)

# Reconstruct the tensor from a state vector.
_reconstruct_tens(::Val{:iso_4}, ::TensND.AbstractTens, u) =
    TensND.TensISO{3}(u[1], u[2])
_reconstruct_tens(::Val{:iso_2}, ::TensND.AbstractTens, u) =
    TensND.TensISO{3}(u[1])
_reconstruct_tens(::Val{:ti}, proto::TensND.AbstractTens, u) =
    TensND._rebuild(proto, ntuple(i -> u[i], length(u)))
_reconstruct_tens(::Val{:ortho}, proto::TensND.AbstractTens, u) =
    TensND._rebuild(proto, ntuple(i -> u[i], length(u)))
# Fully-anisotropic fallback : rebuild from the Mandel 6×6 form.
_reconstruct_tens(::Val{:full_4}, ::TensND.AbstractTens, u) =
    TensND.inv_KM(reshape(u, 6, 6))
_reconstruct_tens(::Val{:full_2}, ::TensND.AbstractTens, u) =
    TensND.inv_KM(reshape(u, 3, 3))

# ── Per-phase contribution helpers ──────────────────────────────────────────

# Dilute correction `(C_i − C) ⊡ A_dil(C)` for a solid inclusion phase
# (symmetrize honoured through `_phase_dilute_concentration`).
function _diff_dilute_correction(
        rve::RVE, name::Symbol, prop::Symbol,
        P_curr::TensND.AbstractTens{4, 3}; kw...
    )
    P_i = phase_property(rve, name, prop)
    A = _phase_dilute_concentration(rve, name, prop, P_curr; kw...)
    return (P_i - P_curr) ⊡ A
end

function _diff_dilute_correction(
        rve::RVE, name::Symbol, prop::Symbol,
        P_curr::TensND.AbstractTens{2, 3}; kw...
    )
    P_i = phase_property(rve, name, prop)
    A = _phase_dilute_concentration(rve, name, prop, P_curr; kw...)
    return (P_i - P_curr) ⋅ A
end

# Crack contribution kernel **per unit density** : returns
# `delta_stiffness(geom, N, 1.0)` (or `delta_conductivity` for 2-tensor)
# so the RHS multiplies by `dε/dτ` directly.
function _diff_crack_density_kernel(
        rve::RVE, name::Symbol, prop::Symbol,
        P_curr::TensND.AbstractTens{4, 3}; kw...
    )
    geom = rve.phases[name].geometry
    sym = phase_symmetrize(rve, name)
    P₀_proj = _project_matrix(P_curr, sym)
    K_int = _crack_interface_K4(rve, name)
    N = MFH_Core.stiffness_contribution(
        geom, P₀_proj;
        K_interface = K_int, kw...
    )
    return _apply_symmetrize(MFH_Core.delta_stiffness(geom, N, 1.0), sym)
end

function _diff_crack_density_kernel(
        rve::RVE, name::Symbol, prop::Symbol,
        P_curr::TensND.AbstractTens{2, 3}; kw...
    )
    geom = rve.phases[name].geometry
    sym = phase_symmetrize(rve, name)
    P₀_proj = _project_matrix(P_curr, sym)
    α_int = _crack_interface_α(rve, name)
    N = MFH_Core.conductivity_contribution(
        geom, P₀_proj;
        α_interface = α_int, kw...
    )
    return _apply_symmetrize(MFH_Core.delta_conductivity(geom, N, 1.0), sym)
end
