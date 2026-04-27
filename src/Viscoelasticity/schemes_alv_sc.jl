# =============================================================================
#  schemes_alv_sc.jl — time-domain self-consistent (SC) homogenisation.
#
#  Iterates the symmetric SC fixed point on the discrete `(6n × 6n)`
#  effective relaxation matrix:
#
#    C̃_{m+1} = (Σ_α f_α C̃_α ∘ Ã_α^dil(C̃_m)) ∘ (Σ_α f_α Ã_α^dil(C̃_m))^{-vol}
#
#  where the sum runs over **all** phases (the matrix included), and the
#  dilute concentration `Ã_α^dil(C̃_m)` is computed against the running
#  estimate `C̃_m` itself.  Picard iteration with optional damping;
#  convergence on the Frobenius norm of the residual `C̃_{m+1} − C̃_m`.
#
#  Reference: Sanahuja IJSS 2013 §3.2 ; Barthélémy et al. IJES 2019 §4 ;
#  ECHOES manual `viscoelasticity_time.qmd` § "SC ALV scheme".
# =============================================================================

"""
    self_consistent_alv(rve, prop; times,
                        abstol = 1e-10, reltol = 1e-8, maxiters = 200,
                        damping = 0.0, verbose = false,
                        select_best = false) -> Matrix

Self-consistent ALV homogenisation.  Iterates the symmetric Picard
fixed point on the `(6n × 6n)` block matrix until convergence.

The initial estimate is the discretised matrix kernel `C̃^0`. Each
iteration rebuilds the per-phase Hill kernels using the current
estimate's iso parameters, computes the dilute concentration tensors,
and forms `C̃_{m+1}`.

Returns the converged effective relaxation matrix.

# Keyword arguments

- `abstol`     — absolute Frobenius tolerance on `‖C̃_{m+1} − C̃_m‖`.
- `reltol`     — additive relative tolerance (multiplied by `‖C̃_m‖`).
- `maxiters`   — hard iteration cap.
- `damping`    — Picard relaxation `0 ≤ damping < 1` (0 = no damping).
- `verbose`    — print residual norms each iteration.
- `select_best`— return the best iterate seen (rather than the last)
  when convergence stalls.
"""
function self_consistent_alv(rve::RVE, prop::Symbol;
                             times::AbstractVector{<:Real},
                             abstol::Real = 1.0e-10,
                             reltol::Real = 1.0e-8,
                             maxiters::Int = 200,
                             damping::Real = 0.0,
                             verbose::Bool = false,
                             select_best::Bool = false)
    # 1. Discretise every phase's kernel once.
    C_M_law = matrix_property(rve, prop)
    C_M_law isa ViscoLaw ||
        throw(ArgumentError("self_consistent_alv: matrix property is not a ViscoLaw"))
    C_0 = _trapezoidal_relaxation(C_M_law, times, 6)
    f_M = matrix_volume_fraction(rve)
    incl_names = inclusion_phase_names(rve)
    C_phases = Matrix{eltype(C_0)}[C_0]
    geometries = Any[rve.phases[rve.matrix_name].geometry]
    fractions = Float64[f_M]
    symmetrizes = AbstractSymmetrize[NoSymmetrize()]
    crack_data = Tuple{Any, Float64, AbstractSymmetrize}[]
    for name in incl_names
        ph = rve.phases[name]
        a = rve.amounts[name]
        if a isa CrackDensity
            ph.geometry isa MFH_Core.AbstractCrack ||
                throw(ArgumentError("self_consistent_alv: phase $name has CrackDensity but geometry is not a crack"))
            push!(crack_data, (ph.geometry, Float64(a.value),
                                phase_symmetrize(rve, name)))
            continue
        end
        C_r_law = phase_property(rve, name, prop)
        C_r_law isa ViscoLaw ||
            throw(ArgumentError("self_consistent_alv: phase $name property is not a ViscoLaw"))
        push!(C_phases, _trapezoidal_relaxation(C_r_law, times, 6))
        push!(geometries, ph.geometry)
        push!(fractions, _amount_value(rve, name))
        push!(symmetrizes, phase_symmetrize(rve, name))
    end

    # 2. Pre-compute the Mandel forms of U^A, V^A for each phase.
    U_M_phases = Matrix{Float64}[_tens_to_mandel66(tens_UA(g)) for g in geometries]
    V_M_phases = Matrix{Float64}[_tens_to_mandel66(tens_VA(g)) for g in geometries]

    # 3. Iterate.
    n = length(times)
    sz = 6 * n
    Id = _identity_alv(n, eltype(C_0))
    C_m = copy(C_0)
    best_resid = Inf
    C_best = C_m

    for iter in 1:maxiters
        C_m_new = _sc_alv_step(C_m, C_phases, U_M_phases, V_M_phases,
                               fractions, n, Id, symmetrizes)
        # Crack contribution (Budiansky-O'Connell SC):
        # ΔJ̃_cracks evaluated against the **running estimate C_m**,
        # then added to the compliance side of the solid-SC update.
        # `C_n+1 = inv(inv(C_solid_SC(C_m)) + ΔJ_cracks(C_m))`
        if !isempty(crack_data)
            ΔJ = zeros(eltype(C_m), size(C_m)...)
            J_m = volterra_inverse(C_m; block_size = 6)
            @inbounds for (geom, ε, sym) in crack_data
                Ñ = stiffness_contribution_alv_at(geom, C_m)
                ΔC = delta_stiffness_alv(geom, Ñ, ε)
                ΔJ_block = -(J_m * ΔC * J_m)
                ΔJ_block = _maybe_symmetrize_alv(ΔJ_block, sym)
                ΔJ .+= ΔJ_block
            end
            J_solid_new = volterra_inverse(C_m_new; block_size = 6)
            C_m_new = volterra_inverse(J_solid_new .+ ΔJ; block_size = 6)
        end
        Δ = norm(C_m_new - C_m)
        norm_C = norm(C_m)
        tol_eff = abstol + reltol * norm_C
        verbose && @info "SC-ALV iter $iter : ‖Δ‖ = $(Δ)   tol = $tol_eff"
        if select_best && Δ < best_resid
            best_resid = Δ
            C_best = C_m_new
        end
        if Δ ≤ tol_eff
            return C_m_new
        end
        # Picard with relaxation.
        C_m = (1 - damping) .* C_m_new .+ damping .* C_m
    end

    @debug "self_consistent_alv: maxiters=$(maxiters) reached without convergence" abstol reltol
    return select_best ? C_best : C_m
end

# Single SC step.
function _sc_alv_step(C_m::AbstractMatrix,
                      C_phases::AbstractVector{<:AbstractMatrix},
                      U_M_phases::AbstractVector{<:AbstractMatrix},
                      V_M_phases::AbstractVector{<:AbstractMatrix},
                      fractions::AbstractVector{<:Real},
                      n::Int, Id::AbstractMatrix,
                      symmetrizes::AbstractVector{<:AbstractSymmetrize})
    sz = size(C_m, 1)
    T = eltype(C_m)
    A_avg = zeros(T, sz, sz)
    CA_avg = zeros(T, sz, sz)
    # Iso parameters of the running estimate → scalar Volterra inverses
    # for the Hill-kernel time-space decoupling.
    α_m, β_m = iso_params_from_blocks(C_m)
    M_long  = @. (α_m + 2 * β_m) / 3
    M_shear = β_m ./ 2
    J_long  = volterra_inverse(M_long;  block_size = 1)
    J_shear = volterra_inverse(M_shear; block_size = 1)

    @inbounds for α in eachindex(C_phases)
        # Phase Hill kernel against current estimate C_m.
        U_M = U_M_phases[α]
        V_M = V_M_phases[α]
        D_M = V_M .- U_M
        P_α = zeros(T, sz, sz)
        for i in 1:n, j in 1:i
            block = J_long[i, j] .* U_M .+ J_shear[i, j] .* D_M
            rows = (6 * (i - 1) + 1):(6 * i)
            cols = (6 * (j - 1) + 1):(6 * j)
            P_α[rows, cols] = block
        end
        # Dilute concentration & scaled contribution.
        ΔC = C_phases[α] - C_m
        A_dil = volterra_inverse(Id + P_α * ΔC; block_size = 6)
        CA   = C_phases[α] * A_dil
        # Apply orientation-averaging projection (`symmetrize=[ISO]` for
        # ECHOES) to both the dilute concentration and its scaled
        # contribution.  This matters per SC iteration so that the
        # running C_m converges to the iso-symmetrized fixed point.
        sym = symmetrizes[α]
        A_dil = _maybe_symmetrize_alv(A_dil, sym)
        CA    = _maybe_symmetrize_alv(CA,    sym)
        f = fractions[α]
        @. A_avg += f * A_dil
        @. CA_avg += f * CA
    end
    return CA_avg * volterra_inverse(A_avg; block_size = 6)
end
