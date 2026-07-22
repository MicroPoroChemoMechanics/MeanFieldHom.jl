# =============================================================================
#  schemes_alv_extra.jl — Ponte-Castañeda & Willis (PCW),
#  Asymmetric Self-Consistent (ASC) and Differential (DIFF) schemes
#  in ageing linear viscoelasticity.
#
#  All operate on the discrete `(6n × 6n)` block matrices produced by
#  `trapezoidal_matrix` (or its `_trapezoidal_relaxation` wrapper for
#  `:creep`-mode laws).  The Volterra product is the regular matrix
#  product (`*`); the Volterra inverse is [`volterra_inverse`](@ref).
#
#  References:
#    * PCW : Ponte-Castañeda & Willis 1995 — coincides with Maxwell in
#      the single-shape case.
#    * ASC : C++ ECHOES `homogenization_asc.h` ; same fixed point as
#      `self_consistent_alv` but the iteration update is anchored on
#      the matrix property `C_M` rather than on the running estimate.
#    * DIFF : Norris 1985 ; user's hand-written DEM N-component note.
#      Solved as a SciML ODE on the fictitious incorporation time
#      `τ ∈ [0, 1]` (`Tsit5` default ; user-overridable via the
#      `alg = …` kwarg of [`DifferentialScheme`](@ref)).
# =============================================================================

# ── PCW : single-shape distribution ⇒ algebraically identical to Maxwell ───

"""
    pcw_alv(C_0, contribs, fractions; H_dist) -> Matrix

Ponte-Castañeda & Willis (1995) viscoelastic homogenization in
single-distribution-shape form.  The formula is algebraically
identical to [`maxwell_alv`](@ref) ; the only difference is that the
Hill kernel `H_dist` is computed against the **distribution shape**
stored in `rve.distribution_shape`, not against any individual phase.

Use [`homogenize_alv`](@ref) with `scheme = PonteCastanedaWillis()` —
the dispatcher reads the RVE's distribution shape and forwards here.
"""
function pcw_alv(
        C_0::AbstractMatrix,
        contribs::AbstractVector{<:AbstractMatrix},
        fractions::AbstractVector;
        H_dist::AbstractMatrix
    )
    return maxwell_alv(C_0, contribs, fractions; H_0 = H_dist)
end

# ── Asymmetric Self-Consistent (ASC) — stiffness form ──────────────────────

"""
    asymmetric_self_consistent_alv(rve::RVE, prop::Symbol; times,
                                    abstol = 1e-10, reltol = 1e-8,
                                    maxiters = 200, damping = 0.0,
                                    verbose = false, select_best = false)
        -> Matrix{T}

Asymmetric self-consistent viscoelastic homogenization.  The
iteration update reads

    `C^{n+1} = C_M + Σ_i f_i (C_i − C_M) ∘ A^{dil,i}(C^n)`,

mirroring the C++ ECHOES `homogenization_asc.h` form.  Returns the
`(6n × 6n)` effective relaxation matrix once the residual
`‖C^{n+1} − C^n‖_F` falls below `abstol + reltol · ‖C^n‖_F` (or after
`maxiters` Picard steps).
"""
function asymmetric_self_consistent_alv(
        rve::RVE, prop::Symbol;
        times::AbstractVector{<:Real},
        abstol::Real = 1.0e-10,
        reltol::Real = 1.0e-8,
        maxiters::Int = 200,
        damping::Real = 0.0,
        verbose::Bool = false,
        select_best::Bool = false
    )
    C_M_law = matrix_property(rve, prop)
    C_M_law isa ViscoLaw ||
        throw(ArgumentError("asymmetric_self_consistent_alv: matrix property is not a ViscoLaw"))
    C_M = _trapezoidal_relaxation(C_M_law, times, 6)
    incl_names = inclusion_phase_names(rve)
    # Eltype-generic containers (Dual-safe — see schemes_alv_sc.jl).
    C_phases = Matrix[]
    geometries = Any[]
    fractions = typeof(matrix_volume_fraction(rve))[]   # eltype of the RVE amounts
    symmetrizes = AbstractSymmetrize[]
    crack_data = Tuple{
        Any, Any, AbstractSymmetrize,
        Union{Nothing, AbstractMatrix},
        Union{Nothing, AbstractMatrix},
    }[]
    for name in incl_names
        ph = rve.phases[name]
        a = rve.amounts[name]
        if a isa CrackDensity
            Rn_mat = haskey(ph.properties, :Rn) ?
                _trapezoidal_relaxation_scalar(ph.properties[:Rn], times) : nothing
            Rt_mat = haskey(ph.properties, :Rt) ?
                _trapezoidal_relaxation_scalar(ph.properties[:Rt], times) : nothing
            push!(
                crack_data, (
                    ph.geometry, a.value,
                    phase_symmetrize(rve, name), Rn_mat, Rt_mat,
                )
            )
            continue
        end
        C_r_law = phase_property(rve, name, prop)
        C_r_law isa ViscoLaw ||
            throw(ArgumentError("asymmetric_self_consistent_alv: phase $name property is not a ViscoLaw"))
        push!(C_phases, _trapezoidal_relaxation(C_r_law, times, 6))
        push!(geometries, ph.geometry)
        push!(fractions, _amount_value(rve, name))
        push!(symmetrizes, phase_symmetrize(rve, name))
    end

    U_M_phases = Matrix[_tens_to_mandel66(tens_UA(g)) for g in geometries]
    V_M_phases = Matrix[_tens_to_mandel66(tens_VA(g)) for g in geometries]

    Tp = _alv_promoted_eltype(
        vcat(Matrix[C_M], C_phases), fractions, U_M_phases, crack_data
    )
    n = length(times)
    Id = _identity_alv(n, Tp)
    C_n = Tp.(C_M)
    best_resid = Inf
    C_best = C_n

    for iter in 1:maxiters
        C_n_new = _asc_alv_step(
            C_M, C_n, C_phases, U_M_phases, V_M_phases,
            fractions, symmetrizes, n, Id
        )
        # Crack contribution (Budiansky-O'Connell SC):
        # `ΔJ̃_cracks(C_n)` against the running estimate, added to the
        # compliance side of the ASC solid update.
        if !isempty(crack_data)
            ΔJ = zeros(eltype(C_n), size(C_n)...)
            J_n = volterra_inverse(C_n; block_size = 6)
            @inbounds for (geom, ε, sym, Rn_mat, Rt_mat) in crack_data
                Ñ = stiffness_contribution_alv_at(
                    geom, C_n;
                    Rn_mat = Rn_mat,
                    Rt_mat = Rt_mat
                )
                ΔC = delta_stiffness_alv(geom, Ñ, ε)
                ΔJ_block = -(J_n * ΔC * J_n)
                ΔJ_block = _maybe_symmetrize_alv(ΔJ_block, sym)
                ΔJ .+= ΔJ_block
            end
            J_solid_new = volterra_inverse(C_n_new; block_size = 6)
            C_n_new = volterra_inverse(J_solid_new .+ ΔJ; block_size = 6)
        end
        Δ = norm(C_n_new - C_n)
        norm_C = norm(C_n)
        tol_eff = abstol + reltol * norm_C
        verbose && @info "ASC-ALV iter $iter : ‖Δ‖ = $(Δ)   tol = $tol_eff"
        if select_best && Δ < best_resid
            best_resid = Δ
            C_best = C_n_new
        end
        if Δ ≤ tol_eff
            return C_n_new
        end
        C_n = (1 - damping) .* C_n_new .+ damping .* C_n
    end

    @debug "asymmetric_self_consistent_alv: maxiters=$(maxiters) reached without convergence" abstol reltol
    return select_best ? C_best : C_n
end

# Single ASC step — increment is anchored on C_M.
function _asc_alv_step(
        C_M::AbstractMatrix,
        C_n::AbstractMatrix,
        C_phases::AbstractVector{<:AbstractMatrix},
        U_M_phases::AbstractVector{<:AbstractMatrix},
        V_M_phases::AbstractVector{<:AbstractMatrix},
        fractions::AbstractVector{<:Real},
        symmetrizes::AbstractVector{<:AbstractSymmetrize},
        n::Int, Id::AbstractMatrix
    )
    sz = size(C_n, 1)
    T = eltype(C_n)
    Δ = zeros(T, sz, sz)
    α_n, β_n = iso_params_from_blocks(C_n)
    M_long = @. (α_n + 2 * β_n) / 3
    M_shear = β_n ./ 2
    J_long = volterra_inverse(M_long; block_size = 1)
    J_shear = volterra_inverse(M_shear; block_size = 1)

    @inbounds for r in eachindex(C_phases)
        U_M = U_M_phases[r]
        V_M = V_M_phases[r]
        D_M = V_M .- U_M
        P_r = zeros(T, sz, sz)
        for i in 1:n, j in 1:i
            block = J_long[i, j] .* U_M .+ J_shear[i, j] .* D_M
            rows = (6 * (i - 1) + 1):(6 * i)
            cols = (6 * (j - 1) + 1):(6 * j)
            P_r[rows, cols] = block
        end
        ΔCr = C_phases[r] - C_n
        A_dil = volterra_inverse(Id + P_r * ΔCr; block_size = 6)
        contrib = (C_phases[r] - C_M) * A_dil
        sym = symmetrizes[r]
        contrib = _maybe_symmetrize_alv(contrib, sym)
        f = fractions[r]
        @. Δ += f * contrib
    end
    return C_M .+ Δ
end

# =============================================================================
#  Differential ALV — SciML ODE on the fictitious incorporation time τ.
#
#  State : `vec(C̃::Matrix{Float64})` of length `(6n)²`.
#  RHS    : reshape `u → C̃`, build the per-phase Hill kernel against
#           `C̃`, assemble Norris terms via Sherman-Morrison and crack
#           density terms, return the flattened tensor derivative.
#  Solver : `OrdinaryDiffEq.solve` with `Tsit5` default.
# =============================================================================

"""
    differential_alv(rve::RVE, prop::Symbol; times,
                      nsteps = 100, trajectory = nothing,
                      abstol = 1e-8, reltol = 1e-6, alg = nothing) -> Matrix{T}

Differential homogenization in ageing linear viscoelasticity, solved
as a SciML ODE on the fictitious incorporation time `τ ∈ [0, 1]`
([Norris 1985](@cite norris1985); user's hand-written DEM note) :

```math
\\frac{\\mathrm d \\tilde{\\mathbb C}}{\\mathrm d \\tau}
  = \\sum_\\alpha \\frac{\\mathrm d \\varphi_\\alpha}{\\mathrm d \\tau}
                  (\\tilde{\\mathbb C}_\\alpha - \\tilde{\\mathbb C})
                  \\circ \\tilde{\\mathbb A}_\\alpha^{dil}(\\tilde{\\mathbb C})
   + \\sum_c \\frac{\\mathrm d \\varepsilon_c}{\\mathrm d \\tau}
              \\Delta\\tilde{\\mathbb C}^{crack}_c(\\tilde{\\mathbb C})
```

with the volume balance `df = (𝟙 − f ⊗ 𝐔)·dφ` inverted by Sherman-
Morrison for solid phases (cracks contribute their density derivative
directly).  `nsteps` is the density of save points along τ ; the
integration step is controlled by `abstol` / `reltol`.
"""
function differential_alv(
        rve::RVE, prop::Symbol;
        times::AbstractVector{<:Real},
        nsteps::Int = 100,
        trajectory = nothing,
        abstol::Real = 1.0e-8,
        reltol::Real = 1.0e-6,
        alg = nothing
    )
    C_M_law = matrix_property(rve, prop)
    C_M_law isa ViscoLaw ||
        throw(ArgumentError("differential_alv: matrix property is not a ViscoLaw"))
    C_M_full = _trapezoidal_relaxation(C_M_law, times, 6)
    n = length(times)

    # Per-phase data : split solids vs cracks.
    solid_data = NamedTuple[]
    crack_data = NamedTuple[]
    for name in inclusion_phase_names(rve)
        ph = rve.phases[name]
        amt = rve.amounts[name]
        if amt isa Schemes.VolumeFraction
            C_r_law = phase_property(rve, name, prop)
            C_r_law isa ViscoLaw ||
                throw(ArgumentError("differential_alv: phase $name property is not a ViscoLaw"))
            push!(
                solid_data, (
                    name = name,
                    C_r = _trapezoidal_relaxation(C_r_law, times, 6),
                    geom = ph.geometry,
                    target = amt.value,
                    sym = phase_symmetrize(rve, name),
                    U_M = _tens_to_mandel66(tens_UA(ph.geometry)),
                    V_M = _tens_to_mandel66(tens_VA(ph.geometry)),
                )
            )
        else  # CrackDensity
            ph.geometry isa MFH_Core.AbstractCrack ||
                throw(ArgumentError("differential_alv: phase $name has CrackDensity but geometry $(typeof(ph.geometry)) is not a crack"))
            Rn_mat = haskey(ph.properties, :Rn) ?
                _trapezoidal_relaxation_scalar(ph.properties[:Rn], times) : nothing
            Rt_mat = haskey(ph.properties, :Rt) ?
                _trapezoidal_relaxation_scalar(ph.properties[:Rt], times) : nothing
            push!(
                crack_data, (
                    name = name,
                    geom = ph.geometry,
                    target = amt.value,
                    sym = phase_symmetrize(rve, name),
                    Rn_mat = Rn_mat,
                    Rt_mat = Rt_mat,
                )
            )
        end
    end

    # Trajectory : default proportional path.
    paths = trajectory === nothing ?
        _resolve_paths_alv(Schemes.Proportional(), rve, nsteps) :
        _resolve_paths_alv(trajectory, rve, nsteps)

    # ODE state and parameters — the state carries the promoted eltype of
    # every input (Dual-safe, cf. `_alv_promoted_eltype`).
    Tp = eltype(C_M_full)
    for sd in solid_data
        Tp = promote_type(Tp, eltype(sd.C_r), typeof(sd.target), eltype(sd.U_M))
    end
    for cd in crack_data
        Tp = promote_type(Tp, typeof(cd.target))
        cd.Rn_mat === nothing || (Tp = promote_type(Tp, eltype(cd.Rn_mat)))
        cd.Rt_mat === nothing || (Tp = promote_type(Tp, eltype(cd.Rt_mat)))
    end
    sz = 6 * n
    x0 = vec(Tp.(C_M_full))
    ode_p = (
        n = n, sz = sz,
        solid_data = solid_data,
        crack_data = crack_data,
        paths = paths,
    )
    rhs! = (du, u, p, τ) -> _diff_alv_ode_rhs!(du, u, p, τ)
    prob = ODEProblem(rhs!, x0, (0.0, 1.0), ode_p)
    sol = solve(
        prob,
        alg === nothing ? Tsit5() : alg;
        abstol, reltol,
        saveat = range(0.0, 1.0; length = max(nsteps, 1) + 1),
        dense = false
    )
    return reshape(sol.u[end], sz, sz)
end

# ── ALV ODE RHS ─────────────────────────────────────────────────────────────

function _diff_alv_ode_rhs!(du, u, p, τ)
    n, sz = p.n, p.sz
    C_curr = reshape(u, sz, sz)
    Δ = zeros(eltype(u), sz, sz)
    Id = _identity_alv(n, eltype(u))

    n_solid = length(p.solid_data)
    if n_solid > 0
        # Sherman-Morrison : dφ_α/dτ = df_α/dτ + (f_α / f_0) · sum(df).
        f = Vector{eltype(u)}(undef, n_solid)
        df = Vector{eltype(u)}(undef, n_solid)
        @inbounds for (i, r) in enumerate(p.solid_data)
            nt = p.paths[r.name]
            f[i] = nt.f(τ) * r.target
            df[i] = nt.df(τ) * r.target
        end
        f0 = one(eltype(u)) - sum(f; init = zero(eltype(u)))
        sum_df = sum(df; init = zero(eltype(u)))

        # Pre-compute Volterra inverses against C_curr (shared across solid phases).
        α_c, β_c = iso_params_from_blocks(C_curr)
        M_long = @. (α_c + 2 * β_c) / 3
        M_shear = β_c ./ 2
        J_long = volterra_inverse(M_long; block_size = 1)
        J_shear = volterra_inverse(M_shear; block_size = 1)

        @inbounds for (i, r) in enumerate(p.solid_data)
            dφᵢ = df[i] + (f[i] / f0) * sum_df
            iszero(dφᵢ) && continue
            # Per-phase Hill kernel against C_curr.
            D_M = r.V_M .- r.U_M
            P_r = zeros(eltype(u), sz, sz)
            for ii in 1:n, jj in 1:ii
                block = J_long[ii, jj] .* r.U_M .+ J_shear[ii, jj] .* D_M
                rows = (6 * (ii - 1) + 1):(6 * ii)
                cols = (6 * (jj - 1) + 1):(6 * jj)
                P_r[rows, cols] = block
            end
            ΔC = r.C_r - C_curr
            A_dil = volterra_inverse(Id + P_r * ΔC; block_size = 6)
            contrib = ΔC * A_dil
            contrib = _maybe_symmetrize_alv(contrib, r.sym)
            @. Δ += dφᵢ * contrib
        end
    end

    # Crack contributions : dε_α/dτ · ΔC̃^crack_α(C_curr) where the
    # crack ΔC̃^crack is the dilute stiffness contribution evaluated
    # against the running matrix `C_curr` (with optional Sevostianov
    # interface stiffness correction).
    @inbounds for r in p.crack_data
        nt = p.paths[r.name]
        dε = nt.df(τ) * r.target
        iszero(dε) && continue
        Ñ = stiffness_contribution_alv_at(
            r.geom, C_curr;
            Rn_mat = r.Rn_mat, Rt_mat = r.Rt_mat
        )
        ΔC = delta_stiffness_alv(r.geom, Ñ, 1.0)
        ΔC = _maybe_symmetrize_alv(ΔC, r.sym)
        @. Δ += dε * ΔC
    end

    du .= vec(Δ)
    return nothing
end

# ── Trajectory path resolution for ALV — share the elastic
#    `_resolve_paths` (which now returns callables, not vectors) ──────────────
function _resolve_paths_alv(
        trajectory::Schemes.DifferentialTrajectory,
        rve::RVE, nsteps::Int
    )
    return Schemes._resolve_paths(trajectory, rve, nsteps)
end
