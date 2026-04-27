# =============================================================================
#  schemes_alv_extra.jl — Ponte-Castañeda & Willis (PCW),
#  Asymmetric Self-Consistent (ASC) and Differential (DIFF) schemes
#  in ageing linear viscoelasticity.
#
#  All three operate on the discrete `(6n × 6n)` block matrices produced
#  by `trapezoidal_matrix` (or its `_trapezoidal_relaxation` wrapper for
#  `:creep`-mode laws).  The Volterra product is the regular matrix
#  product (`*`); the Volterra inverse is [`volterra_inverse`](@ref).
#
#  References:
#    * PCW : Ponte-Castañeda & Willis 1995 — coincides with Maxwell in
#      the single-shape case (see `Schemes/pcw.jl`).
#    * ASC : C++ ECHOES `homogenization_asc.h` ; same fixed point as
#      `self_consistent_alv` but the iteration update is anchored on the
#      matrix property `C_M` rather than on the running estimate.
#    * DIFF : Norris 1985 ; ECHOES `homogenization_differential::compute_property`.
# =============================================================================

# ── PCW : single-shape distribution ⇒ algebraically identical to Maxwell ───

"""
    pcw_alv(C_0, contribs, fractions; H_dist) -> Matrix

Ponte-Castañeda & Willis (1995) viscoelastic homogenisation in
single-distribution-shape form.  The formula is algebraically
identical to [`maxwell_alv`](@ref) ; the only difference is that the
Hill kernel `H_dist` is computed against the **distribution shape**
stored in `rve.distribution_shape`, not against any individual phase.

Use [`homogenize_alv`](@ref) with `scheme = PonteCastanedaWillis()` —
the dispatcher reads the RVE's distribution shape and forwards here.
"""
function pcw_alv(C_0::AbstractMatrix,
                  contribs::AbstractVector{<:AbstractMatrix},
                  fractions::AbstractVector;
                  H_dist::AbstractMatrix)
    return maxwell_alv(C_0, contribs, fractions; H_0 = H_dist)
end

# ── Asymmetric Self-Consistent (ASC) — stiffness form ──────────────────────
#
# C^{n+1} = C_M + Σ_i f_i (C_i − C_M) ∘ A^{dil,i}(C^n)
#
# where the dilute concentration `A^{dil,i}(C^n) = (𝟙 + P̃_i(C^n) ∘ (C_i − C^n))^{-vol}`
# is built against the running estimate `C^n` but the increment is
# applied to the matrix property `C_M`.  Same fixed point as the
# Hill-symmetric SC (`self_consistent_alv`) but a different basin of
# attraction and convergence rate.

"""
    asymmetric_self_consistent_alv(rve::RVE, prop::Symbol; times,
                                    abstol = 1e-10, reltol = 1e-8,
                                    maxiters = 200, damping = 0.0,
                                    verbose = false, select_best = false)
        -> Matrix{T}

Asymmetric self-consistent viscoelastic homogenisation.  The
iteration update reads

    `C^{n+1} = C_M + Σ_i f_i (C_i − C_M) ∘ A^{dil,i}(C^n)`,

mirroring the C++ ECHOES `homogenization_asc.h` form.  Returns the
`(6n × 6n)` effective relaxation matrix once the residual
`‖C^{n+1} − C^n‖_F` falls below `abstol + reltol · ‖C^n‖_F` (or after
`maxiters` Picard steps).
"""
function asymmetric_self_consistent_alv(rve::RVE, prop::Symbol;
                                          times::AbstractVector{<:Real},
                                          abstol::Real = 1.0e-10,
                                          reltol::Real = 1.0e-8,
                                          maxiters::Int = 200,
                                          damping::Real = 0.0,
                                          verbose::Bool = false,
                                          select_best::Bool = false)
    C_M_law = matrix_property(rve, prop)
    C_M_law isa ViscoLaw ||
        throw(ArgumentError("asymmetric_self_consistent_alv: matrix property is not a ViscoLaw"))
    C_M = _trapezoidal_relaxation(C_M_law, times, 6)
    incl_names = inclusion_phase_names(rve)
    C_phases = Matrix{eltype(C_M)}[]
    geometries = Any[]
    fractions = Float64[]
    symmetrizes = AbstractSymmetrize[]
    crack_data = Tuple{Any, Float64, AbstractSymmetrize}[]
    for name in incl_names
        ph = rve.phases[name]
        a = rve.amounts[name]
        if a isa CrackDensity
            push!(crack_data, (ph.geometry, Float64(a.value),
                                phase_symmetrize(rve, name)))
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

    U_M_phases = Matrix{Float64}[_tens_to_mandel66(tens_UA(g)) for g in geometries]
    V_M_phases = Matrix{Float64}[_tens_to_mandel66(tens_VA(g)) for g in geometries]

    n = length(times)
    Id = _identity_alv(n, eltype(C_M))
    C_n = copy(C_M)
    best_resid = Inf
    C_best = C_n

    for iter in 1:maxiters
        C_n_new = _asc_alv_step(C_M, C_n, C_phases, U_M_phases, V_M_phases,
                                 fractions, symmetrizes, n, Id)
        # Crack contribution (Budiansky-O'Connell SC):
        # `ΔJ̃_cracks(C_n)` against the running estimate, added to the
        # compliance side of the ASC solid update.
        if !isempty(crack_data)
            ΔJ = zeros(eltype(C_n), size(C_n)...)
            J_n = volterra_inverse(C_n; block_size = 6)
            @inbounds for (geom, ε, sym) in crack_data
                Ñ = stiffness_contribution_alv_at(geom, C_n)
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
function _asc_alv_step(C_M::AbstractMatrix,
                        C_n::AbstractMatrix,
                        C_phases::AbstractVector{<:AbstractMatrix},
                        U_M_phases::AbstractVector{<:AbstractMatrix},
                        V_M_phases::AbstractVector{<:AbstractMatrix},
                        fractions::AbstractVector{<:Real},
                        symmetrizes::AbstractVector{<:AbstractSymmetrize},
                        n::Int, Id::AbstractMatrix)
    sz = size(C_n, 1)
    T = eltype(C_n)
    Δ = zeros(T, sz, sz)
    α_n, β_n = iso_params_from_blocks(C_n)
    M_long  = @. (α_n + 2 * β_n) / 3
    M_shear = β_n ./ 2
    J_long  = volterra_inverse(M_long;  block_size = 1)
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

# ── Differential scheme — explicit Euler integration of the Norris ODE ──
#
# At each step k:
#   1. Build the volume-balance system  M Φ = ΔF  (same as elastic).
#   2. For every solid phase i, add  Φ_i · (C_i − C_curr) ∘ A^{dil,i}(C_curr)
#      to the running matrix C_curr (Volterra products on (6n)×(6n)).
#   3. For every crack phase, add a (1 / nsteps) fractional dilute
#      contribution evaluated at C_curr.

"""
    differential_alv(rve::RVE, prop::Symbol; times,
                      nsteps = 100, trajectory = nothing) -> Matrix{T}

Differential homogenisation in ageing linear viscoelasticity.
Integrates the Norris ODE
    `dC̃ / df_i = (C̃_i − C̃) ∘ A^{dil,i}(C̃)`
along the trajectory defined by `trajectory` over `nsteps` explicit
Euler steps, using the **same volume-balance recurrence** as the
elastic [`DifferentialScheme`](@ref) of `Schemes/differential.jl`.

When no trajectory is provided, the proportional path is used by
default — i.e. every phase fraction is filled at the same fractional
rate, mirroring the C++ ECHOES default.
"""
function differential_alv(rve::RVE, prop::Symbol;
                            times::AbstractVector{<:Real},
                            nsteps::Int = 100,
                            trajectory = nothing)
    C_M_law = matrix_property(rve, prop)
    C_M_law isa ViscoLaw ||
        throw(ArgumentError("differential_alv: matrix property is not a ViscoLaw"))
    C_curr = _trapezoidal_relaxation(C_M_law, times, 6)
    n = length(times)
    Id = _identity_alv(n, eltype(C_curr))

    # Inclusion data (separate solid vs crack phases).
    incl_names = inclusion_phase_names(rve)
    solid_data  = NamedTuple[]
    crack_data  = NamedTuple[]
    for name in incl_names
        ph = rve.phases[name]
        C_r_law = phase_property(rve, name, prop)
        C_r_law isa ViscoLaw ||
            throw(ArgumentError("differential_alv: phase $name property is not a ViscoLaw"))
        amt = rve.amounts[name]
        record = (
            name      = name,
            C_r       = _trapezoidal_relaxation(C_r_law, times, 6),
            geom      = ph.geometry,
            target    = Float64(amt.value),
            sym       = phase_symmetrize(rve, name),
            U_M       = _tens_to_mandel66(tens_UA(ph.geometry)),
            V_M       = _tens_to_mandel66(tens_VA(ph.geometry)),
        )
        if amt isa Schemes.VolumeFraction
            push!(solid_data, record)
        else
            push!(crack_data, record)
        end
    end

    # Trajectory : default proportional path (linear from 0 to 1 in nsteps).
    paths = if trajectory === nothing
        Dict{Symbol, Vector{Float64}}(
            r.name => collect(range(0.0, 1.0; length = nsteps + 1))
            for r in solid_data)
    else
        _resolve_paths_alv(trajectory, solid_data, nsteps)
    end

    n_solid = length(solid_data)
    for k in 1:nsteps
        if n_solid > 0
            f_prev = [paths[r.name][k]     * r.target for r in solid_data]
            f_curr = [paths[r.name][k + 1] * r.target for r in solid_data]
            ΔF = f_curr .- f_prev
            Mlin = Matrix{Float64}(LinearAlgebra.I, n_solid, n_solid)
            for i in 1:n_solid
                Mlin[i, :] .-= f_prev[i]
            end
            Φ = Mlin \ ΔF

            # Pre-compute Volterra inverses against C_curr for the Hill kernel
            # (shared across phases at this step).
            α_c, β_c = iso_params_from_blocks(C_curr)
            M_long  = @. (α_c + 2 * β_c) / 3
            M_shear = β_c ./ 2
            J_long  = volterra_inverse(M_long;  block_size = 1)
            J_shear = volterra_inverse(M_shear; block_size = 1)

            sz = size(C_curr, 1)
            T = eltype(C_curr)
            for (i, r) in enumerate(solid_data)
                Φi = Φ[i]
                iszero(Φi) && continue
                # Phase Hill kernel against current matrix.
                D_M = r.V_M .- r.U_M
                P_r = zeros(T, sz, sz)
                @inbounds for ii in 1:n, jj in 1:ii
                    block = J_long[ii, jj] .* r.U_M .+ J_shear[ii, jj] .* D_M
                    rows = (6 * (ii - 1) + 1):(6 * ii)
                    cols = (6 * (jj - 1) + 1):(6 * jj)
                    P_r[rows, cols] = block
                end
                ΔC = r.C_r - C_curr
                A_dil = volterra_inverse(Id + P_r * ΔC; block_size = 6)
                contrib = ΔC * A_dil
                contrib = _maybe_symmetrize_alv(contrib, r.sym)
                @. C_curr += Φi * contrib
            end
        end
        # Cracks not yet routed through the differential ALV pipeline —
        # their compliance contribution is geometry-only and would require
        # `compliance_contribution_alv` integration per step.  Deferred.
    end
    return C_curr
end

# Resolve the path map for differential ALV given a trajectory descriptor.
# Currently supports `Schemes.Proportional` only — extend with `Sequential`
# / `CustomPath` analogues when needed.
function _resolve_paths_alv(::Schemes.Proportional, solid_data,
                              nsteps::Int)
    return Dict{Symbol, Vector{Float64}}(
        r.name => collect(range(0.0, 1.0; length = nsteps + 1))
        for r in solid_data)
end
