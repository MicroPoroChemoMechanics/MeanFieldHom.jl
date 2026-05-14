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
function self_consistent_alv(
        rve::RVE, prop::Symbol;
        times::AbstractVector{<:Real},
        abstol::Real = 1.0e-10,
        reltol::Real = 1.0e-8,
        maxiters::Int = 200,
        damping::Real = 0.0,
        verbose::Bool = false,
        select_best::Bool = false
    )
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
    # crack_data tuple: (geom, density, sym, Rn_mat::Union{Nothing, Matrix},
    #                    Rt_mat::Union{Nothing, Matrix}).  The interface
    # matrices are pre-discretised once on the times grid and reused by
    # every SC iteration when computing the crack stiffness contribution
    # against the running estimate C_m.
    crack_data = Tuple{
        Any, Float64, AbstractSymmetrize,
        Union{Nothing, Matrix{Float64}},
        Union{Nothing, Matrix{Float64}},
    }[]
    for name in incl_names
        ph = rve.phases[name]
        a = rve.amounts[name]
        if a isa CrackDensity
            ph.geometry isa MFH_Core.AbstractCrack ||
                throw(ArgumentError("self_consistent_alv: phase $name has CrackDensity but geometry is not a crack"))
            Rn_mat = haskey(ph.properties, :Rn) ?
                _trapezoidal_relaxation_scalar(ph.properties[:Rn], times) : nothing
            Rt_mat = haskey(ph.properties, :Rt) ?
                _trapezoidal_relaxation_scalar(ph.properties[:Rt], times) : nothing
            push!(
                crack_data, (
                    ph.geometry, Float64(a.value),
                    phase_symmetrize(rve, name), Rn_mat, Rt_mat,
                )
            )
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
        # ECHOES SC body (cf. `homogenization_scheme.h::evaluate` and
        # `inclusion_ellipsoid::compute_strain_Stress`,
        # `inclusion_crack::compute_strain_Stress_void_crack`):
        #   strain_Stress_α  = A_α(C_m) · J_m   (solid, J_m = inv(C_m))
        #   strain_Stress_c  = sym(H_c(C_m))    (void crack — no J_m)
        #   stress_Stress_α  = C_α · strain_Stress_α
        #   stress_Stress_c  = 0                (traction-free)
        # Accumulators :  A_E = Σ f_α·sym(strain_Stress_α)
        #                 B_E = Σ f_α·sym(stress_Stress_α)
        # Result : C_eff = B_E · A_E^{-vol}.
        # The trailing `J_m` cancels for solid-only RVEs but NOT when
        # cracks are present (the crack term has no `J_m` factor),
        # giving the ECHOES SC fixed point that doesn't match the
        # textbook `(Σ f·C·A)·(Σ f·A)^{-1}` form.
        if isempty(crack_data)
            C_m_new = _sc_alv_step(
                C_m, C_phases, U_M_phases, V_M_phases,
                fractions, n, Id, symmetrizes
            )
        else
            C_m_new = _sc_alv_step_echoes_form(
                C_m, C_phases,
                U_M_phases, V_M_phases,
                fractions, n, Id, symmetrizes,
                crack_data
            )
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

# ── ECHOES SC body for ALV with cracks ─────────────────────────────────────
#
# Mirrors the elastic ECHOES SC body (cf.
# `homogenization_scheme.h::evaluate` and the C++
# `inclusion_*::compute_strain_Stress` family) :
#   strain_Stress_α  = A_α(C_m) · J_m   (solid, J_m = volterra-inv(C_m))
#   strain_Stress_c  = sym(H_c(C_m))    (void crack — NO trailing J_m)
#   stress_Stress_α  = C_α · strain_Stress_α
#   stress_Stress_c  = 0                (traction-free)
# Accumulators :  A_E = Σ f_α·sym(strain_Stress_α)
#                 B_E = Σ f_α·sym(stress_Stress_α)
# Result   : C_eff = B_E · A_E^{-vol}.  The trailing `J_m` factor
# cancels for solid-only RVEs (recovering `(Σ f·CA)·(Σ f·A)^{-vol}`)
# but not for cracks, whose `strain_Stress` is the bare compliance
# contribution `H_c` without `J_m` factor.  This is the ECHOES SC
# fixed point.
function _sc_alv_step_echoes_form(
        C_m::AbstractMatrix,
        C_phases::AbstractVector{<:AbstractMatrix},
        U_M_phases::AbstractVector{<:AbstractMatrix},
        V_M_phases::AbstractVector{<:AbstractMatrix},
        fractions::AbstractVector{<:Real},
        n::Int, Id::AbstractMatrix,
        symmetrizes::AbstractVector{<:AbstractSymmetrize},
        crack_data
    )
    sz = size(C_m, 1)
    T = eltype(C_m)
    A_avg = zeros(T, sz, sz)   # = Σ f·sym(A_α)            (no J_m yet)
    CA_avg = zeros(T, sz, sz)   # = Σ f·sym(C_α·A_α)        (no J_m yet)
    α_m, β_m = iso_params_from_blocks(C_m)
    M_long = @. (α_m + 2 * β_m) / 3
    M_shear = β_m ./ 2
    J_long = volterra_inverse(M_long; block_size = 1)
    J_shear = volterra_inverse(M_shear; block_size = 1)
    @inbounds for α in eachindex(C_phases)
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
        ΔC = C_phases[α] - C_m
        A_dil = volterra_inverse(Id + P_α * ΔC; block_size = 6)
        sym = symmetrizes[α]
        A_dil_sym = _maybe_symmetrize_alv(A_dil, sym)
        CA_sym = _maybe_symmetrize_alv(C_phases[α] * A_dil, sym)
        f = T(fractions[α])
        @. A_avg += f * A_dil_sym
        @. CA_avg += f * CA_sym
    end
    # Crack compliance contributions (without J_m factor — that's the
    # essential ECHOES form difference compared to solids).
    H_total = _build_sc_crack_extra_J(C_m, crack_data)
    # Apply the ECHOES `B · A^{-vol}` formula with explicit `J_m`.
    J_m = volterra_inverse(C_m; block_size = 6)
    A_E = (A_avg * J_m) .+ H_total
    B_E = CA_avg * J_m
    return B_E * volterra_inverse(A_E; block_size = 6)
end

# ── Crack compliance contribution against the running estimate ─────────────
#
# Computes `Σ_c ε·sym(H̃_c(C_m))` where `H̃_c(C_m)` is the (Sevostianov-
# corrected) crack compliance contribution against the running estimate
# `C_m`.  Used in two ways:
#   (a) `_build_sc_crack_extra_J(C_m, crack_data)`  — appended to the
#       compliance `J̃ = inv(C_m_solid_SC)` in the Budiansky-O'Connell
#       branch of the SC iteration (the default, robust path).
#   (b) `_build_sc_crack_extra_A(C_m, crack_data)`  — alias used by the
#       experimental ECHOES-form `B · A^{-vol}` MT body (kept available
#       for the Newton-Raphson SC solver, which uses the same algebra
#       but solves `F(C) = 0` instead of iterating Picard).
function _build_sc_crack_extra_J(C_m::AbstractMatrix, crack_data)
    sz = size(C_m, 1)
    T = eltype(C_m)
    extra = zeros(T, sz, sz)
    isempty(crack_data) && return extra
    _is_iso_block(C_m) ||
        error("self_consistent_alv with cracks: only iso running estimate is supported")
    α_c, β_c = _iso_pair(C_m)
    α_p_2β = α_c .+ 2β_c
    α_p_βh = α_c .+ β_c ./ 2
    α_p_β = α_c .+ β_c
    βα1 = β_c * α_p_βh
    βα2 = β_c * α_p_β
    B_n_base = (8 / (3π)) .* volterra_left_divide(βα1, α_p_2β)
    B_t_base = (32 / (9π)) .* volterra_left_divide(βα2, α_p_2β)
    Iₙ = Matrix{T}(LinearAlgebra.I, size(α_c, 1), size(α_c, 1))
    @inbounds for (geom, ε, sym, Rn_mat, Rt_mat) in crack_data
        B_n = B_n_base
        B_t = B_t_base
        if Rn_mat !== nothing || Rt_mat !== nothing
            b = semi_minor(geom)
            if Rn_mat !== nothing
                KB = Rn_mat * B_n; @. KB *= b; @. KB += Iₙ
                B_n = B_n * volterra_inverse(KB; block_size = 1)
            end
            if Rt_mat !== nothing
                KB = Rt_mat * B_t; @. KB *= b; @. KB += Iₙ
                B_t = B_t * volterra_inverse(KB; block_size = 1)
            end
        end
        Z = zeros(T, size(α_c))
        ℓ₁ = (3 / 4) .* B_n
        ℓ₆ = (3 / 8) .* B_t
        H_TI = ti_blocks_from_params((ℓ₁, copy(Z), copy(Z), copy(Z), copy(Z), ℓ₆))
        H_full = _maybe_symmetrize_alv(delta_compliance_alv(geom, H_TI, ε), sym)
        @. extra += H_full
    end
    return extra
end

# Alias used by the ECHOES-form MT body (numerically identical : both
# represent the iso-symmetrized crack compliance contribution scaled
# by the Budiansky concentration factor `(4π/3)·ε`).
@inline _build_sc_crack_extra_A(C_m, crack_data) = _build_sc_crack_extra_J(C_m, crack_data)

# Single SC step using the MT body with `m_ref_is_matrix = true`
# (matrix special-cased as `f_M·𝟙` / `f_M·C_m`, solids iterated
# against C_m, cracks added to A).  This is exactly the ECHOES MT body
# applied iteratively with the running estimate as reference — a
# strict implicit-MT view of SC.
function _sc_alv_mt_body_against(
        C_m::AbstractMatrix,
        fractions::AbstractVector{<:Real},
        f_M::Real,
        C_phases::AbstractVector{<:AbstractMatrix},
        U_M_phases::AbstractVector{<:AbstractMatrix},
        V_M_phases::AbstractVector{<:AbstractMatrix},
        n::Int, Id::AbstractMatrix,
        symmetrizes::AbstractVector{<:AbstractSymmetrize},
        crack_data
    )
    sz = size(C_m, 1)
    T = eltype(C_m)
    A = T(f_M) .* Id
    B = T(f_M) .* C_m
    # Solid inclusions iterated with the Hill kernel computed against C_m.
    α_m, β_m = iso_params_from_blocks(C_m)
    M_long = @. (α_m + 2 * β_m) / 3
    M_shear = β_m ./ 2
    J_long = volterra_inverse(M_long; block_size = 1)
    J_shear = volterra_inverse(M_shear; block_size = 1)
    @inbounds for s in 2:length(C_phases)
        U_M = U_M_phases[s]
        V_M = V_M_phases[s]
        D_M = V_M .- U_M
        P_s = zeros(T, sz, sz)
        for i in 1:n, j in 1:i
            block = J_long[i, j] .* U_M .+ J_shear[i, j] .* D_M
            rows = (6 * (i - 1) + 1):(6 * i)
            cols = (6 * (j - 1) + 1):(6 * j)
            P_s[rows, cols] = block
        end
        ΔC = C_phases[s] - C_m
        A_dil = volterra_inverse(Id + P_s * ΔC; block_size = 6)
        sym = symmetrizes[s]
        AC = A_dil * C_m
        AC = _maybe_symmetrize_alv(AC, sym)
        CAC = C_phases[s] * AC
        CAC = _maybe_symmetrize_alv(CAC, sym)
        f = T(fractions[s])
        @. A += f * AC
        @. B += f * CAC
    end
    # Cracks : recomputed against C_m at every iteration.
    @inbounds for (geom, ε, sym, Rn_mat, Rt_mat) in crack_data
        _is_iso_block(C_m) ||
            error("self_consistent_alv with cracks: only iso running estimate is supported")
        α_c, β_c = _iso_pair(C_m)
        α_p_2β = α_c .+ 2β_c
        α_p_βh = α_c .+ β_c ./ 2
        α_p_β = α_c .+ β_c
        βα1 = β_c * α_p_βh
        βα2 = β_c * α_p_β
        B_n = (8 / (3π)) .* volterra_left_divide(βα1, α_p_2β)
        B_t = (32 / (9π)) .* volterra_left_divide(βα2, α_p_2β)
        if Rn_mat !== nothing || Rt_mat !== nothing
            Iₙ = Matrix{T}(LinearAlgebra.I, size(α_c, 1), size(α_c, 1))
            b = semi_minor(geom)
            if Rn_mat !== nothing
                KB = Rn_mat * B_n; @. KB *= b; @. KB += Iₙ
                B_n = B_n * volterra_inverse(KB; block_size = 1)
            end
            if Rt_mat !== nothing
                KB = Rt_mat * B_t; @. KB *= b; @. KB += Iₙ
                B_t = B_t * volterra_inverse(KB; block_size = 1)
            end
        end
        Z = zeros(T, size(α_c))
        ℓ₁ = (3 / 4) .* B_n
        ℓ₆ = (3 / 8) .* B_t
        H_TI = ti_blocks_from_params((ℓ₁, copy(Z), copy(Z), copy(Z), copy(Z), ℓ₆))
        H_full = _maybe_symmetrize_alv(delta_compliance_alv(geom, H_TI, ε), sym)
        AC = H_full * C_m
        @. A += AC
        # Traction-free: B contribution = 0.
    end
    return B * volterra_inverse(A; block_size = 6)
end

# Single SC step (legacy MFH form, retained for external callers /
# internal sub-steps that still rely on the standard
# `(Σ f A) ↔ (Σ f C A)` accumulators).
function _sc_alv_step(
        C_m::AbstractMatrix,
        C_phases::AbstractVector{<:AbstractMatrix},
        U_M_phases::AbstractVector{<:AbstractMatrix},
        V_M_phases::AbstractVector{<:AbstractMatrix},
        fractions::AbstractVector{<:Real},
        n::Int, Id::AbstractMatrix,
        symmetrizes::AbstractVector{<:AbstractSymmetrize};
        extra_A::Union{Nothing, AbstractMatrix} = nothing
    )
    sz = size(C_m, 1)
    T = eltype(C_m)
    A_avg = zeros(T, sz, sz)
    CA_avg = zeros(T, sz, sz)
    if extra_A !== nothing
        @. A_avg += extra_A
    end
    # Iso parameters of the running estimate → scalar Volterra inverses
    # for the Hill-kernel time-space decoupling.
    α_m, β_m = iso_params_from_blocks(C_m)
    M_long = @. (α_m + 2 * β_m) / 3
    M_shear = β_m ./ 2
    J_long = volterra_inverse(M_long; block_size = 1)
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
        CA = C_phases[α] * A_dil
        # Apply orientation-averaging projection (`symmetrize=[ISO]` for
        # ECHOES) to both the dilute concentration and its scaled
        # contribution.  This matters per SC iteration so that the
        # running C_m converges to the iso-symmetrized fixed point.
        sym = symmetrizes[α]
        A_dil = _maybe_symmetrize_alv(A_dil, sym)
        CA = _maybe_symmetrize_alv(CA, sym)
        f = fractions[α]
        @. A_avg += f * A_dil
        @. CA_avg += f * CA
    end
    return CA_avg * volterra_inverse(A_avg; block_size = 6)
end

# ECHOES-form SC step: matches `homogenization_scheme.h::evaluate`
# (lines 866-901) verbatim.
#
#   A_E = f_M · 𝟙 + Σ_solids f_s · strain_Strain_s + Σ_cracks ε_c · strain_Strain_c
#   B_E = f_M · X + Σ_solids f_s · stress_Strain_s + Σ_cracks 0 (TF)
#   C_m+1 = B_E · A_E^{-vol}
#
# where X = C_m (running estimate) and `strain_Strain = A_εε · C_m`,
# `stress_Strain = C_inc · A_εε · C_m`.  Cracks contribute the
# `(4π/3)·ε·H̃_iso(C_m)·C_m` term to `A_E`.
function _sc_alv_step_echoes(
        C_m::AbstractMatrix,
        C_phases::AbstractVector{<:AbstractMatrix},
        U_M_phases::AbstractVector{<:AbstractMatrix},
        V_M_phases::AbstractVector{<:AbstractMatrix},
        fractions::AbstractVector{<:Real},
        f_M::Real,
        n::Int, Id::AbstractMatrix,
        symmetrizes::AbstractVector{<:AbstractSymmetrize},
        crack_data
    )
    sz = size(C_m, 1)
    T = eltype(C_m)
    # Matrix part.
    A_E = T(f_M) .* Id
    B_E = T(f_M) .* C_m
    # Iso parameters of the running estimate for the Hill-kernel
    # time-space decoupling.
    α_m, β_m = iso_params_from_blocks(C_m)
    M_long = @. (α_m + 2 * β_m) / 3
    M_shear = β_m ./ 2
    J_long = volterra_inverse(M_long; block_size = 1)
    J_shear = volterra_inverse(M_shear; block_size = 1)
    # Solid inclusions (loop indices 2..end correspond to solids; the
    # matrix lives at C_phases[1] but is handled separately above).
    @inbounds for α in 2:length(C_phases)
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
        ΔC = C_phases[α] - C_m
        A_dil = volterra_inverse(Id + P_α * ΔC; block_size = 6)
        AC = A_dil * C_m                           # strain_Strain_s
        CAC = C_phases[α] * AC                     # stress_Strain_s
        sym = symmetrizes[α]
        AC = _maybe_symmetrize_alv(AC, sym)
        CAC = _maybe_symmetrize_alv(CAC, sym)
        f = T(fractions[α])
        @. A_E += f * AC
        @. B_E += f * CAC
    end
    # Cracks (traction-free or interface-stiffness Sevostianov).
    @inbounds for (geom, ε, sym, Rn_mat, Rt_mat) in crack_data
        _is_iso_block(C_m) ||
            error("self_consistent_alv with cracks: only iso running estimate is supported")
        α_c, β_c = _iso_pair(C_m)
        α_p_2β = α_c .+ 2β_c
        α_p_βh = α_c .+ β_c ./ 2
        α_p_β = α_c .+ β_c
        βα1 = β_c * α_p_βh
        βα2 = β_c * α_p_β
        B_n = (8 / (3π)) .* volterra_left_divide(βα1, α_p_2β)
        B_t = (32 / (9π)) .* volterra_left_divide(βα2, α_p_2β)
        if Rn_mat !== nothing || Rt_mat !== nothing
            Iₙ = Matrix{T}(LinearAlgebra.I, size(α_c, 1), size(α_c, 1))
            b = semi_minor(geom)
            if Rn_mat !== nothing
                KB = Rn_mat * B_n; @. KB *= b; @. KB += Iₙ
                B_n = B_n * volterra_inverse(KB; block_size = 1)
            end
            if Rt_mat !== nothing
                KB = Rt_mat * B_t; @. KB *= b; @. KB += Iₙ
                B_t = B_t * volterra_inverse(KB; block_size = 1)
            end
        end
        Z = zeros(T, size(α_c))
        ℓ₁ = (3 / 4) .* B_n
        ℓ₆ = (3 / 8) .* B_t
        H_TI = ti_blocks_from_params((ℓ₁, copy(Z), copy(Z), copy(Z), copy(Z), ℓ₆))
        H_full = _maybe_symmetrize_alv(delta_compliance_alv(geom, H_TI, ε), sym)
        AC = H_full * C_m                          # strain_Strain_crack
        @. A_E += AC
        # B_E: traction-free crack → C_inc = 0 → no contribution.
    end
    return B_E * volterra_inverse(A_E; block_size = 6)
end
