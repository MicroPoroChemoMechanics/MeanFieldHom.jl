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

"""
    homogenize_alv(rve, scheme, prop::Symbol; times) -> Matrix

ALV pipeline: build the discrete `(6n × 6n)` block matrices of every
phase, compute the ALV Hill kernel and dilute concentration tensors,
and dispatch on `scheme` to the corresponding `_alv` scheme function.

Returns the effective relaxation matrix `C̃_eff` of size `(6n × 6n)`,
where `n = length(times)`.

Supports two inclusion-geometry families:
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
    C_0 = trapezoidal_matrix(C_M_law, times)
    f_M = matrix_volume_fraction(rve)

    # 2. Loop on inclusions.
    incl_names = inclusion_phase_names(rve)
    fractions = Float64[]
    contribs = Matrix{eltype(C_0)}[]
    A_duts = Matrix{eltype(C_0)}[]
    C_phases = Matrix{eltype(C_0)}[C_0]
    H_phases = Matrix{eltype(C_0)}[]   # per-phase Hill kernels (for Maxwell distribution)
    for name in incl_names
        ph = rve.phases[name]
        C_r_law = phase_property(rve, name, prop)
        C_r, A_dut, N_dut, P_r = _inclusion_alv_quantities(
            ph.geometry, C_r_law, C_M_law, C_0, times)
        push!(C_phases, C_r)
        push!(A_duts, A_dut)
        push!(contribs, N_dut)
        push!(H_phases, P_r)
        push!(fractions, _amount_value(rve, name))
    end

    return _homogenize_alv_dispatch(rve, scheme, prop, times,
                                    C_0, C_phases, A_duts, contribs,
                                    H_phases, fractions, f_M; kw...)
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
    C_r = trapezoidal_matrix(C_r_law, times)
    P_r = hill_kernel(geom, C_M_law, times)
    A_dut = dilute_concentration_alv(C_r, C_0, P_r)
    N_dut = dilute_contribution_alv(C_r, C_0, P_r)
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

# ── Dispatch table on scheme types ──────────────────────────────────────────

function _homogenize_alv_dispatch(::RVE, ::Voigt, ::Symbol, ::AbstractVector,
                                  C_0, C_phases, A_duts, contribs,
                                  H_phases, fractions, f_M; kw...)
    return voigt_alv(C_phases, [f_M; fractions])
end

function _homogenize_alv_dispatch(::RVE, ::Reuss, ::Symbol, ::AbstractVector,
                                  C_0, C_phases, A_duts, contribs,
                                  H_phases, fractions, f_M; kw...)
    return reuss_alv(C_phases, [f_M; fractions])
end

function _homogenize_alv_dispatch(::RVE, ::Dilute, ::Symbol, ::AbstractVector,
                                  C_0, C_phases, A_duts, contribs,
                                  H_phases, fractions, f_M; kw...)
    return dilute_alv(C_0, contribs, fractions)
end

function _homogenize_alv_dispatch(::RVE, ::DiluteDual, ::Symbol, ::AbstractVector,
                                  C_0, C_phases, A_duts, contribs,
                                  H_phases, fractions, f_M; kw...)
    return dilute_dual_alv(C_0, contribs, fractions)
end

function _homogenize_alv_dispatch(::RVE, ::MoriTanaka, ::Symbol, ::AbstractVector,
                                  C_0, C_phases, A_duts, contribs,
                                  H_phases, fractions, f_M; kw...)
    return mori_tanaka_alv(C_0, A_duts, contribs, fractions, f_M)
end

function _homogenize_alv_dispatch(rve::RVE, ::Maxwell, ::Symbol,
                                  times::AbstractVector,
                                  C_0, C_phases, A_duts, contribs,
                                  H_phases, fractions, f_M; kw...)
    # Default distribution shape: spherical (matches the elastic Maxwell
    # default in `Schemes.maxwell`).
    C_M_law = matrix_property(rve, :C)
    H_0 = hill_kernel(Spheroid(1.0), C_M_law, times)
    return maxwell_alv(C_0, contribs, fractions; H_0 = H_0)
end

# Self-Consistent ALV: re-routes to `self_consistent_alv` (different
# computation flow, since each iteration recomputes the per-phase Hill
# kernels against the running estimate).
function _homogenize_alv_dispatch(rve::RVE, sc::SelfConsistent, prop::Symbol,
                                  times::AbstractVector,
                                  C_0, C_phases, A_duts, contribs,
                                  H_phases, fractions, f_M; kw...)
    return self_consistent_alv(rve, prop; times = times,
                               sc.options..., kw...)
end
