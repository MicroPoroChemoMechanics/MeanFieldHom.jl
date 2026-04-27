# =============================================================================
#  homogenize_alv.jl — public dispatcher routing ViscoLaw properties
#  through the ALV homogenisation pipeline.
#
#  Usage : `homogenize(rve, scheme, :C; times = T)` where `T` is a
#  `Vector{<:Real}` of monotonically increasing time points.  The
#  function detects whether `phase_property(rve, _, :C) isa ViscoLaw`
#  and, if so, calls `homogenize_alv` instead of the elastic pipeline.
# =============================================================================

"""
    has_visco_property(rve, prop::Symbol = :C) -> Bool

Return `true` if any phase of `rve` carries a [`ViscoLaw`](@ref) under
the key `prop`.
"""
function has_visco_property(rve::RVE, prop::Symbol = :C)
    try
        m = matrix_property(rve, prop)
        m isa ViscoLaw && return true
    catch
    end
    for name in inclusion_phase_names(rve)
        try
            v = phase_property(rve, name, prop)
            v isa ViscoLaw && return true
        catch
        end
    end
    return false
end

# ── Iso projection of a 6n×6n block matrix (ECHOES `symmetrize=[ISO]`) ─────
#
# ECHOES applies an orientation-averaging projection to each phase
# 4-tensor whose RVE definition carries `symmetrize=[ISO]`.  For each
# 6×6 Mandel block, the iso projection is :
#   α = (1/3) Σᵢⱼ T_iijj  =  (M[1,1] + M[2,2] + M[3,3] + 2(M[1,2] + M[1,3] + M[2,3])) / 3
#   β = (Σᵢⱼ T_ijij - α) / 5  =  (trace(M) - α) / 5
#   M_iso = α 𝕁 + β 𝕂  (rebuilt via `iso_blocks_from_params` on a 1×1
#                       parameter pair).
#
# We apply this block-by-block to the (6n)×(6n) `Ñ` (and `Ã`) of the
# inclusions when their phase's `phase_symmetrize` is `IsoSymmetrize`.
# Iso averaging of a TI block over a uniform orientation distribution
# matches ECHOES `symmetrize=[ISO]` exactly.

@inline function _iso_project_mandel66(M::AbstractMatrix)
    @assert size(M) == (6, 6)
    α = (M[1, 1] + M[2, 2] + M[3, 3] +
         2 * (M[1, 2] + M[1, 3] + M[2, 3])) / 3
    tr = M[1, 1] + M[2, 2] + M[3, 3] + M[4, 4] + M[5, 5] + M[6, 6]
    β = (tr - α) / 5
    return α, β
end

"""
    _iso_project_blocks(M::AbstractMatrix) -> Matrix

Project every 6×6 Mandel block of a `(6n × 6n)` ALV matrix to its iso
component (Reynolds average over the orthogonal group), returning a
new `(6n × 6n)` block matrix whose every block is iso.  Equivalent to
the ECHOES `symmetrize=[ISO]` orientation-averaging projection
applied to each `(t_i, t_j)` block independently.
"""
function _iso_project_blocks(M::AbstractMatrix)
    sz = size(M, 1)
    sz == size(M, 2) ||
        throw(ArgumentError("_iso_project_blocks: matrix must be square"))
    sz % 6 == 0 ||
        throw(ArgumentError("_iso_project_blocks: size $(sz) not divisible by 6"))
    n = sz ÷ 6
    T = eltype(M)
    α = zeros(T, n, n)
    β = zeros(T, n, n)
    @inbounds for i in 1:n, j in 1:n
        rows = (6 * (i - 1) + 1):(6 * i)
        cols = (6 * (j - 1) + 1):(6 * j)
        a, b = _iso_project_mandel66(view(M, rows, cols))
        α[i, j] = a
        β[i, j] = b
    end
    return iso_blocks_from_params(α, β)
end

"""
    _maybe_symmetrize_alv(M, sym) -> Matrix

Apply the orientation-averaging projection corresponding to `sym` to a
`(6n × 6n)` ALV block matrix.  Currently supports `NoSymmetrize`
(passthrough) and `IsoSymmetrize` (block-by-block iso projection).
TI projection on the ALV side is reserved for a follow-up.
"""
@inline _maybe_symmetrize_alv(M::AbstractMatrix, ::NoSymmetrize) = M
@inline _maybe_symmetrize_alv(M::AbstractMatrix, ::IsoSymmetrize) =
    _iso_project_blocks(M)

"""
    _trapezoidal_relaxation(law::ViscoLaw, times, B) -> Matrix

Build the discrete relaxation block matrix from a `ViscoLaw`, regardless
of whether the law is in `:relaxation` or `:creep` mode.  When the law
is `:creep`, the trapezoidal compliance matrix is inverted (block forward
substitution at `block_size = B`) to obtain the corresponding relaxation
matrix — the convention every ALV scheme assumes internally.

`B` is the block size (`6` for order-4 4-tensor / Mandel, `3` for
order-2 vector-tensor).
"""
function _trapezoidal_relaxation(law::ViscoLaw,
                                  times::AbstractVector{<:Real}, B::Int)
    M = trapezoidal_matrix(law, times)
    if visco_mode(law) === :creep
        return volterra_inverse(M; block_size = B)
    end
    return M
end

"""
    _alv_property_order(law::ViscoLaw, t) -> Int

Inspect the sample returned by `visco_eval(law, t, t)` and report the
tensor order (`2` for vector-tensor / 3×3, `4` for 4-tensor / 6×6
Mandel).  Used by [`homogenize_alv`](@ref) to dispatch between the
order-4 (stiffness / relaxation) and order-2 (conductivity / creep
admittance) pipelines.
"""
function _alv_property_order(law::ViscoLaw, t::Real)
    sample = visco_eval(law, t, t)
    if sample isa TensND.AbstractTens{2, 3}
        return 2
    elseif sample isa TensND.AbstractTens{4, 3}
        return 4
    elseif sample isa AbstractMatrix
        if size(sample) == (3, 3)
            return 2
        elseif size(sample) == (6, 6)
            return 4
        end
    end
    throw(ArgumentError("homogenize_alv: cannot infer ALV order from sample of type $(typeof(sample))"))
end

"""
    homogenize_alv(rve, scheme, prop::Symbol; times) -> Matrix

ALV pipeline: build the discrete block matrices of every phase, compute
the ALV Hill kernel and dilute concentration tensors, and dispatch on
`scheme` to the corresponding `_alv` scheme function.

The function dispatches on the **order of the matrix property** (read
once from the matrix `ViscoLaw` sample type):
  * order-4 (4-tensor / 6×6 Mandel kernel) → returns `(6n × 6n)`
    relaxation matrix following the standard Hill-kernel +
    dilute-concentration pipeline.
  * order-2 (2-tensor / 3×3 kernel; conductivity / diffusion /
    permittivity) → returns `(3n × 3n)` matrix via the order-2 ALV
    pipeline (`hill_kernel_order2`, time-space decoupling).

Supports two inclusion-geometry families (order-4 only):
  * Single-shape ellipsoidal (`Ellipsoid`, `Spheroid`): the standard
    Hill-kernel + dilute-concentration pipeline.
  * `LayeredSphere`: bulk + shear ALV recurrences (see
    [`bulk_localization_alv`](@ref) and
    [`shear_localization_alv`](@ref)) feed
    [`stiffness_contribution_alv`](@ref) and
    [`strain_strain_loc_alv`](@ref) directly — no Hill kernel is
    needed.  In this case `phase_property(rve, name, :C)` is ignored
    (the per-layer moduli stored in the geometry are used instead);
    pass any `ViscoLaw` (e.g. `heaviside_law(C_0)`) as a placeholder.
"""
function homogenize_alv(rve::RVE, scheme::HomogenizationScheme,
                        prop::Symbol; times::AbstractVector{<:Real}, kw...)
    # 1. Matrix kernel.
    C_M_law = matrix_property(rve, prop)
    C_M_law isa ViscoLaw ||
        throw(ArgumentError("homogenize_alv: matrix property $prop is not a ViscoLaw"))

    # Dispatch on the property order (2 vs 4) inferred from the sample.
    order = _alv_property_order(C_M_law, first(times))
    if order == 2
        return _homogenize_alv_order2(rve, scheme, prop;
                                       times = times, kw...)
    end

    C_0 = _trapezoidal_relaxation(C_M_law, times, 6)
    f_M = matrix_volume_fraction(rve)

    # 2. Loop on inclusions, separating SOLIDS (`VolumeFraction`) from
    #    CRACKS (`CrackDensity`).  Cracks contribute ΔC̃_crack to the
    #    numerator of the schemes (no volume → no denominator effect).
    incl_names = inclusion_phase_names(rve)
    fractions = Float64[]
    contribs = Matrix{eltype(C_0)}[]
    A_duts = Matrix{eltype(C_0)}[]
    C_phases = Matrix{eltype(C_0)}[C_0]
    H_phases = Matrix{eltype(C_0)}[]   # per-phase Hill kernels (for Maxwell distribution)
    crack_data = Tuple{Any, Float64, AbstractSymmetrize}[]   # (geom, density, sym)
    ΔC_cracks_M = zeros(eltype(C_0), size(C_0)...)   # cracks-against-C_M sum
    ΔJ_cracks_M = zeros(eltype(C_0), size(C_0)...)   # for Reuss/DiluteDual

    for name in incl_names
        ph = rve.phases[name]
        a = rve.amounts[name]
        sym = phase_symmetrize(rve, name)
        if a isa CrackDensity
            geom = ph.geometry
            geom isa MFH_Core.AbstractCrack ||
                throw(ArgumentError("homogenize_alv: phase $name has CrackDensity but geometry $(typeof(geom)) is not a crack"))
            ε = Float64(a.value)
            push!(crack_data, (geom, ε, sym))
            # Stiffness contribution of the crack against C̃_M.
            Ñ = stiffness_contribution_alv_at(geom, C_0)
            ΔC = delta_stiffness_alv(geom, Ñ, ε)
            ΔC = _maybe_symmetrize_alv(ΔC, sym)
            ΔC_cracks_M .+= ΔC
            # Compliance contribution (for Reuss / DiluteDual).
            H̃ = compliance_contribution_alv(geom, C_M_law, times)
            ΔJ = delta_compliance_alv(geom, H̃, ε)
            ΔJ = _maybe_symmetrize_alv(ΔJ, sym)
            ΔJ_cracks_M .+= ΔJ
        else
            C_r_law = phase_property(rve, name, prop)
            C_r, A_dut, N_dut, P_r = _inclusion_alv_quantities(
                ph.geometry, C_r_law, C_M_law, C_0, times)
            A_dut = _maybe_symmetrize_alv(A_dut, sym)
            N_dut = _maybe_symmetrize_alv(N_dut, sym)
            push!(C_phases, C_r)
            push!(A_duts, A_dut)
            push!(contribs, N_dut)
            push!(H_phases, P_r)
            push!(fractions, _amount_value(rve, name))
        end
    end

    return _homogenize_alv_dispatch(rve, scheme, prop, times,
                                    C_0, C_phases, A_duts, contribs,
                                    H_phases, fractions, f_M;
                                    crack_data = crack_data,
                                    ΔC_cracks_M = ΔC_cracks_M,
                                    ΔJ_cracks_M = ΔJ_cracks_M,
                                    kw...)
end

# ── Per-geometry inclusion quantities ───────────────────────────────────────

"""
    _inclusion_alv_quantities(geom, C_r_law, C_M_law, C_0, times)
        -> (C_r, A_dut, N_dut, P_r)

Compute the four `(6n × 6n)` matrices needed by the ALV scheme
dispatch for a single inclusion of geometry `geom`.  Default method
covers ellipsoidal geometries (Hill kernel + dilute formulas);
specialisations for `LayeredSphere` use the layered-sphere recurrences.
"""
function _inclusion_alv_quantities(geom, C_r_law,
                                    C_M_law::ViscoLaw,
                                    C_0::AbstractMatrix,
                                    times::AbstractVector{<:Real})
    C_r_law isa ViscoLaw ||
        throw(ArgumentError("homogenize_alv: phase property is not a ViscoLaw"))
    C_r = _trapezoidal_relaxation(C_r_law, times, 6)
    P_r = hill_kernel(geom, C_M_law, times)
    # Iso fast path : if every input matrix is iso (typical for
    # spherical inclusions in iso ALV matrix), compute the dilute
    # concentration and contribution as TWO scalar n×n Volterra
    # problems and lift back to (6n × 6n) at the end.  Avoids the
    # generic block-LU on the 6n×6n `(𝟙 + P̃ ∘ ΔC̃)` system.
    if _is_iso_block(C_r) && _is_iso_block(C_0) && _is_iso_block(P_r)
        αβ_E = _iso_pair(C_r)
        αβ_0 = _iso_pair(C_0)
        αβ_P = _iso_pair(P_r)
        αβ_A_dut = dilute_concentration_alv_iso(αβ_E, αβ_0, αβ_P)
        αβ_N_dut = dilute_contribution_alv_iso(αβ_E, αβ_0, αβ_P)
        A_dut = _iso_blocks(αβ_A_dut)
        N_dut = _iso_blocks(αβ_N_dut)
    elseif _is_ti_block(C_r) && _is_ti_block(C_0) && _is_ti_block(P_r)
        # TI fast path: shared canonical axis e_3 across phase, matrix
        # and Hill kernel.  Reduces the dilute concentration to a
        # (2n)×(2n) block-Volterra inverse + 2 scalar Volterra inverses.
        ℓ_E = _ti_pair(C_r)
        ℓ_0 = _ti_pair(C_0)
        ℓ_P = _ti_pair(P_r)
        ℓ_A_dut = dilute_concentration_alv_ti(ℓ_E, ℓ_0, ℓ_P)
        ℓ_N_dut = dilute_contribution_alv_ti(ℓ_E, ℓ_0, ℓ_P)
        A_dut = _ti_blocks(ℓ_A_dut)
        N_dut = _ti_blocks(ℓ_N_dut)
    elseif _is_ortho_block(C_r) && _is_ortho_block(C_0) && _is_ortho_block(P_r)
        # Ortho fast path: shared canonical material frame across phase,
        # matrix and Hill kernel.  Reduces the dilute concentration to a
        # (3n)×(3n) block-Volterra inverse + 3 scalar Volterra inverses.
        o_E = _ortho_pair(C_r)
        o_0 = _ortho_pair(C_0)
        o_P = _ortho_pair(P_r)
        o_A_dut = dilute_concentration_alv_ortho(o_E, o_0, o_P)
        o_N_dut = dilute_contribution_alv_ortho(o_E, o_0, o_P)
        A_dut = _ortho_blocks(o_A_dut)
        N_dut = _ortho_blocks(o_N_dut)
    else
        A_dut = dilute_concentration_alv(C_r, C_0, P_r)
        N_dut = dilute_contribution_alv(C_r, C_0, P_r)
    end
    return (C_r, A_dut, N_dut, P_r)
end

function _inclusion_alv_quantities(sphere::LayeredSphere, _C_r_law,
                                    C_M_law::ViscoLaw,
                                    C_0::AbstractMatrix,
                                    times::AbstractVector{<:Real})
    # Per-layer moduli are stored inside the LayeredSphere geometry; the
    # phase-level C_r_law is ignored (we accept any placeholder so that
    # the existing `add_phase!` API still works).
    A_dut = strain_strain_loc_alv(sphere, C_M_law, times)
    N_dut = stiffness_contribution_alv(sphere, C_M_law, times)
    # No single C_r is well-defined for a layered sphere ; expose the
    # dilute-effective stiffness (C_0 + N_dut) as a representative
    # monolithic estimate so the Voigt / Reuss code paths still type-check
    # (they remain non-physical for a layered inclusion).
    C_r = C_0 .+ N_dut
    # Placeholder Hill kernel (not used by Dilute / MT / Maxwell paths).
    P_r = zeros(eltype(C_0), size(C_0)...)
    return (C_r, A_dut, N_dut, P_r)
end

# Convenience: turn `volume_fraction(rve, name)` into a `Float64` even when
# the amount is wrapped in `VolumeFraction`/`CrackDensity`.
function _amount_value(rve::RVE, name::Symbol)
    amount = rve.amounts[name]
    return Float64(amount.value)
end

# ── Iso-symmetry detection for the scheme fast path ─────────────────────────
#
# When every (6n × 6n) block matrix supplied to the scheme step is in
# iso form, the scheme algebra reduces to two independent scalar n × n
# Volterra problems on (α, β).  This is ~108× cheaper for matrix-matrix
# products and ~18× for inversion compared to the generic block-LU
# `(6n × 6n)` path.

"""
    _try_iso_pairs(matrices) -> Vector{Tuple} or nothing

If every matrix in `matrices` passes the iso-form check
(`_is_iso_block`), return a `Vector` of `(α, β)` `n×n` parameter
tuples extracted from each.  Otherwise return `nothing`.

Used by the scheme fast paths to opt into the iso pipeline only when
all phases really are iso 4-tensors.
"""
function _try_iso_pairs(matrices::AbstractVector{<:AbstractMatrix})
    isempty(matrices) && return Tuple{Matrix{Float64}, Matrix{Float64}}[]
    out = Vector{Tuple{Matrix{eltype(matrices[1])}, Matrix{eltype(matrices[1])}}}()
    for M in matrices
        _is_iso_block(M) || return nothing
        push!(out, _iso_pair(M))
    end
    return out
end

"""
    _try_ti_tuples(matrices) -> Vector{NTuple{6, Matrix}} or nothing

If every matrix passes the TI-form check (`_is_ti_block`), return a
`Vector` of 6-tuples of `n×n` Walpole parameter matrices extracted
from each.  Otherwise return `nothing`.

Iso block matrices automatically satisfy the TI test (iso ⊂ TI), so a
mixed iso/TI phase setup with the canonical axis works.
"""
function _try_ti_tuples(matrices::AbstractVector{<:AbstractMatrix})
    T = isempty(matrices) ? Float64 : eltype(matrices[1])
    isempty(matrices) && return NTuple{6, Matrix{T}}[]
    out = Vector{NTuple{6, Matrix{T}}}()
    for M in matrices
        _is_ti_block(M) || return nothing
        push!(out, _ti_pair(M))
    end
    return out
end

"""
    _try_ortho_tuples(matrices) -> Vector{NTuple{12, Matrix}} or nothing

If every matrix passes the ortho-form check (`_is_ortho_block`), return
a `Vector` of 12-tuples of `n×n` ortho parameter matrices extracted
from each.  Otherwise return `nothing`.

Iso and TI (axis = e₃) block matrices automatically satisfy the ortho
test, so a mixed iso/TI/ortho phase setup with the canonical material
frame works.
"""
function _try_ortho_tuples(matrices::AbstractVector{<:AbstractMatrix})
    T = isempty(matrices) ? Float64 : eltype(matrices[1])
    isempty(matrices) && return NTuple{12, Matrix{T}}[]
    out = Vector{NTuple{12, Matrix{T}}}()
    for M in matrices
        _is_ortho_block(M) || return nothing
        push!(out, _ortho_pair(M))
    end
    return out
end

# ── Dispatch table on scheme types ──────────────────────────────────────────
#
# Crack handling.  Every dispatcher honours the optional kwargs
#    crack_data    :: Vector{Tuple{geom, density, sym}}
#    ΔC_cracks_M   :: pre-aggregated stiffness contribution against C̃_M
#    ΔJ_cracks_M   :: pre-aggregated compliance contribution against C̃_M
# computed once in `homogenize_alv` (see the matrix-reference cracks
# loop above).  When `isempty(crack_data)`, the iso/TI fast paths are
# attempted; otherwise the crack-aware generic path is used.

@inline _has_cracks(kw) = haskey(kw, :crack_data) && !isempty(kw[:crack_data])

function _homogenize_alv_dispatch(::RVE, ::Voigt, ::Symbol, ::AbstractVector,
                                  C_0, C_phases, A_duts, contribs,
                                  H_phases, fractions, f_M; kw...)
    # Voigt ignores cracks (zero-volume convention, mirroring elastic
    # `Schemes/voigt.jl`).  Result depends only on solid volume fractions.
    iso = _try_iso_pairs(C_phases)
    if iso !== nothing
        αβ_eff = voigt_alv_iso(iso, [f_M; fractions])
        return _iso_blocks(αβ_eff)
    end
    ti = _try_ti_tuples(C_phases)
    if ti !== nothing
        ℓ_eff = voigt_alv_ti(ti, [f_M; fractions])
        return _ti_blocks(ℓ_eff)
    end
    ortho = _try_ortho_tuples(C_phases)
    if ortho !== nothing
        o_eff = voigt_alv_ortho(ortho, [f_M; fractions])
        return _ortho_blocks(o_eff)
    end
    return voigt_alv(C_phases, [f_M; fractions])
end

function _homogenize_alv_dispatch(::RVE, ::Reuss, ::Symbol, ::AbstractVector,
                                  C_0, C_phases, A_duts, contribs,
                                  H_phases, fractions, f_M; kw...)
    # Reuss ignores cracks (same convention as elastic `Schemes/reuss.jl`).
    iso = _try_iso_pairs(C_phases)
    if iso !== nothing
        αβ_eff = reuss_alv_iso(iso, [f_M; fractions])
        return _iso_blocks(αβ_eff)
    end
    ti = _try_ti_tuples(C_phases)
    if ti !== nothing
        ℓ_eff = reuss_alv_ti(ti, [f_M; fractions])
        return _ti_blocks(ℓ_eff)
    end
    ortho = _try_ortho_tuples(C_phases)
    if ortho !== nothing
        o_eff = reuss_alv_ortho(ortho, [f_M; fractions])
        return _ortho_blocks(o_eff)
    end
    return reuss_alv(C_phases, [f_M; fractions])
end

function _homogenize_alv_dispatch(::RVE, ::Dilute, ::Symbol, ::AbstractVector,
                                  C_0, C_phases, A_duts, contribs,
                                  H_phases, fractions, f_M; kw...)
    if !_has_cracks(kw)
        iso_contribs = _try_iso_pairs(contribs)
        if iso_contribs !== nothing && _is_iso_block(C_0)
            αβ_0 = _iso_pair(C_0)
            αβ_eff = dilute_alv_iso(αβ_0, iso_contribs, fractions)
            return _iso_blocks(αβ_eff)
        end
        ti_contribs = _try_ti_tuples(contribs)
        if ti_contribs !== nothing && _is_ti_block(C_0)
            ℓ_0 = _ti_pair(C_0)
            ℓ_eff = dilute_alv_ti(ℓ_0, ti_contribs, fractions)
            return _ti_blocks(ℓ_eff)
        end
        ortho_contribs = _try_ortho_tuples(contribs)
        if ortho_contribs !== nothing && _is_ortho_block(C_0)
            o_0 = _ortho_pair(C_0)
            o_eff = dilute_alv_ortho(o_0, ortho_contribs, fractions)
            return _ortho_blocks(o_eff)
        end
        return dilute_alv(C_0, contribs, fractions)
    end
    # Cracks: additive — `C̃_dilute + ΔC̃_cracks_M`.
    return dilute_alv(C_0, contribs, fractions) .+ kw[:ΔC_cracks_M]
end

function _homogenize_alv_dispatch(::RVE, ::DiluteDual, ::Symbol, ::AbstractVector,
                                  C_0, C_phases, A_duts, contribs,
                                  H_phases, fractions, f_M; kw...)
    if _has_cracks(kw)
        # DiluteDual : invert C → J, add per-phase compliance contribs +
        # crack ΔJ̃, invert back.  We rebuild the per-phase compliance
        # contribs from the relaxation contribs : N̄_J = -J_M ∘ N̄_C ∘ J_M.
        J_M = volterra_inverse(C_0; block_size = 6)
        J_eff = copy(J_M)
        @inbounds for (f, N̄) in zip(fractions, contribs)
            term = -(J_M * N̄ * J_M)
            @. J_eff += f * term
        end
        J_eff .+= kw[:ΔJ_cracks_M]
        return volterra_inverse(J_eff; block_size = 6)
    end
    iso_contribs = _try_iso_pairs(contribs)
    if iso_contribs !== nothing && _is_iso_block(C_0)
        αβ_0 = _iso_pair(C_0)
        αβ_eff = dilute_dual_alv_iso(αβ_0, iso_contribs, fractions)
        return _iso_blocks(αβ_eff)
    end
    ti_contribs = _try_ti_tuples(contribs)
    if ti_contribs !== nothing && _is_ti_block(C_0)
        ℓ_0 = _ti_pair(C_0)
        ℓ_eff = dilute_dual_alv_ti(ℓ_0, ti_contribs, fractions)
        return _ti_blocks(ℓ_eff)
    end
    ortho_contribs = _try_ortho_tuples(contribs)
    if ortho_contribs !== nothing && _is_ortho_block(C_0)
        o_0 = _ortho_pair(C_0)
        o_eff = dilute_dual_alv_ortho(o_0, ortho_contribs, fractions)
        return _ortho_blocks(o_eff)
    end
    return dilute_dual_alv(C_0, contribs, fractions)
end

function _homogenize_alv_dispatch(::RVE, ::MoriTanaka, ::Symbol, ::AbstractVector,
                                  C_0, C_phases, A_duts, contribs,
                                  H_phases, fractions, f_M; kw...)
    if !_has_cracks(kw)
        iso_contribs = _try_iso_pairs(contribs)
        iso_A = _try_iso_pairs(A_duts)
        if iso_contribs !== nothing && iso_A !== nothing && _is_iso_block(C_0)
            αβ_0 = _iso_pair(C_0)
            αβ_eff = mori_tanaka_alv_iso(αβ_0, iso_A, iso_contribs, fractions, f_M)
            return _iso_blocks(αβ_eff)
        end
        ti_contribs = _try_ti_tuples(contribs)
        ti_A = _try_ti_tuples(A_duts)
        if ti_contribs !== nothing && ti_A !== nothing && _is_ti_block(C_0)
            ℓ_0 = _ti_pair(C_0)
            ℓ_eff = mori_tanaka_alv_ti(ℓ_0, ti_A, ti_contribs, fractions, f_M)
            return _ti_blocks(ℓ_eff)
        end
        ortho_contribs = _try_ortho_tuples(contribs)
        ortho_A = _try_ortho_tuples(A_duts)
        if ortho_contribs !== nothing && ortho_A !== nothing && _is_ortho_block(C_0)
            o_0 = _ortho_pair(C_0)
            o_eff = mori_tanaka_alv_ortho(o_0, ortho_A, ortho_contribs, fractions, f_M)
            return _ortho_blocks(o_eff)
        end
        return mori_tanaka_alv(C_0, A_duts, contribs, fractions, f_M)
    end
    # Crack-aware MT : append a virtual phase with N̄ = ΔC̃_cracks_M,
    # Ã = 0 (cracks have no volume in the denominator), f = 1.  This
    # injects ΔC̃ into the numerator without polluting the denominator.
    sz = size(C_0, 1)
    T = eltype(C_0)
    extra_A = zeros(T, sz, sz)
    contribs_aug = vcat(contribs, [kw[:ΔC_cracks_M]])
    A_duts_aug = vcat(A_duts, [extra_A])
    fractions_aug = vcat(fractions, [1.0])
    return mori_tanaka_alv(C_0, A_duts_aug, contribs_aug, fractions_aug, f_M)
end

function _homogenize_alv_dispatch(rve::RVE, ::Maxwell, ::Symbol,
                                  times::AbstractVector,
                                  C_0, C_phases, A_duts, contribs,
                                  H_phases, fractions, f_M; kw...)
    # Default distribution shape: spherical (matches the elastic Maxwell
    # default in `Schemes.maxwell`).
    C_M_law = matrix_property(rve, :C)
    H_0 = hill_kernel(Spheroid(1.0), C_M_law, times)
    if !_has_cracks(kw)
        iso_contribs = _try_iso_pairs(contribs)
        if iso_contribs !== nothing && _is_iso_block(C_0) && _is_iso_block(H_0)
            αβ_0 = _iso_pair(C_0)
            αβ_H_0 = _iso_pair(H_0)
            αβ_eff = maxwell_alv_iso(αβ_0, iso_contribs, fractions, αβ_H_0)
            return _iso_blocks(αβ_eff)
        end
        ti_contribs = _try_ti_tuples(contribs)
        if ti_contribs !== nothing && _is_ti_block(C_0) && _is_ti_block(H_0)
            ℓ_0 = _ti_pair(C_0)
            ℓ_H_0 = _ti_pair(H_0)
            ℓ_eff = maxwell_alv_ti(ℓ_0, ti_contribs, fractions, ℓ_H_0)
            return _ti_blocks(ℓ_eff)
        end
        ortho_contribs = _try_ortho_tuples(contribs)
        if ortho_contribs !== nothing && _is_ortho_block(C_0) && _is_ortho_block(H_0)
            o_0 = _ortho_pair(C_0)
            o_H_0 = _ortho_pair(H_0)
            o_eff = maxwell_alv_ortho(o_0, ortho_contribs, fractions, o_H_0)
            return _ortho_blocks(o_eff)
        end
        return maxwell_alv(C_0, contribs, fractions; H_0 = H_0)
    end
    # Crack-aware Maxwell : append cracks to the contribution sum (Σ).
    contribs_aug = vcat(contribs, [kw[:ΔC_cracks_M]])
    fractions_aug = vcat(fractions, [1.0])
    return maxwell_alv(C_0, contribs_aug, fractions_aug; H_0 = H_0)
end

# Self-Consistent ALV: re-routes to `self_consistent_alv` (different
# computation flow, since each iteration recomputes the per-phase Hill
# kernels against the running estimate).
function _homogenize_alv_dispatch(rve::RVE, sc::SelfConsistent, prop::Symbol,
                                  times::AbstractVector,
                                  C_0, C_phases, A_duts, contribs,
                                  H_phases, fractions, f_M; kw...)
    # SC reads cracks directly from the RVE; strip the pre-aggregated
    # crack kwargs that were meant for the simpler scheme dispatchers.
    kw_filt = Iterators.filter(p -> !(p[1] in
                                       (:crack_data, :ΔC_cracks_M, :ΔJ_cracks_M)),
                                kw)
    return self_consistent_alv(rve, prop; times = times,
                               sc.options..., kw_filt...)
end

# Asymmetric Self-Consistent ALV.  Same ingredients as `SelfConsistent`
# but the iteration update is anchored on the matrix property C_M
# rather than the running estimate (cf. `schemes_alv_extra.jl`).
function _homogenize_alv_dispatch(rve::RVE, asc::AsymmetricSelfConsistent, prop::Symbol,
                                  times::AbstractVector,
                                  C_0, C_phases, A_duts, contribs,
                                  H_phases, fractions, f_M; kw...)
    kw_filt = Iterators.filter(p -> !(p[1] in
                                       (:crack_data, :ΔC_cracks_M, :ΔJ_cracks_M)),
                                kw)
    return asymmetric_self_consistent_alv(rve, prop; times = times,
                                            asc.options..., kw_filt...)
end

# Ponte-Castañeda & Willis ALV.  Algebraically identical to Maxwell in
# the single-shape case, but uses the `rve.distribution_shape` for the
# Hill kernel instead of a fixed sphere.
function _homogenize_alv_dispatch(rve::RVE, ::PonteCastanedaWillis, ::Symbol,
                                  times::AbstractVector,
                                  C_0, C_phases, A_duts, contribs,
                                  H_phases, fractions, f_M; kw...)
    C_M_law = matrix_property(rve, :C)
    dist = rve.distribution_shape
    dist isa UniformDistribution ||
        throw(ArgumentError("PCW-ALV: only UniformDistribution is currently supported"))
    H_d = hill_kernel(dist.shape, C_M_law, times)
    if !_has_cracks(kw)
        iso_contribs = _try_iso_pairs(contribs)
        if iso_contribs !== nothing && _is_iso_block(C_0) && _is_iso_block(H_d)
            αβ_0 = _iso_pair(C_0)
            αβ_H = _iso_pair(H_d)
            αβ_eff = maxwell_alv_iso(αβ_0, iso_contribs, fractions, αβ_H)
            return _iso_blocks(αβ_eff)
        end
        ti_contribs = _try_ti_tuples(contribs)
        if ti_contribs !== nothing && _is_ti_block(C_0) && _is_ti_block(H_d)
            ℓ_0 = _ti_pair(C_0)
            ℓ_H = _ti_pair(H_d)
            ℓ_eff = maxwell_alv_ti(ℓ_0, ti_contribs, fractions, ℓ_H)
            return _ti_blocks(ℓ_eff)
        end
        ortho_contribs = _try_ortho_tuples(contribs)
        if ortho_contribs !== nothing && _is_ortho_block(C_0) && _is_ortho_block(H_d)
            o_0 = _ortho_pair(C_0)
            o_H = _ortho_pair(H_d)
            o_eff = maxwell_alv_ortho(o_0, ortho_contribs, fractions, o_H)
            return _ortho_blocks(o_eff)
        end
        return pcw_alv(C_0, contribs, fractions; H_dist = H_d)
    end
    # Crack-aware PCW : same as Maxwell with rve distribution shape.
    contribs_aug = vcat(contribs, [kw[:ΔC_cracks_M]])
    fractions_aug = vcat(fractions, [1.0])
    return pcw_alv(C_0, contribs_aug, fractions_aug; H_dist = H_d)
end

# Differential ALV.  Multi-step Euler integration of the Norris ODE.
function _homogenize_alv_dispatch(rve::RVE, sch::DifferentialScheme, ::Symbol,
                                  times::AbstractVector,
                                  C_0, C_phases, A_duts, contribs,
                                  H_phases, fractions, f_M; kw...)
    nsteps = get(sch.options, :nsteps, 100)
    return differential_alv(rve, :C; times = times,
                              nsteps = nsteps,
                              trajectory = sch.trajectory)
end
