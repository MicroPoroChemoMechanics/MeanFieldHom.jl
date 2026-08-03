# =============================================================================
#  interfaces_laminate.jl — imperfect interfaces of a periodic laminate.
#
#  The five interface types of `LayeredSpheres` are reused UNCHANGED. A planar
#  interface is the degenerate, curvature-free case of the spherical one, and
#  the algebra collapses to two additive terms — which is what makes the
#  laminate the cleanest possible test of the package's interface conventions.
#
#  ── Primal (spring, Kapitza): a field jump ────────────────────────────────
#  `[u] = 𝕂·(σ·n)`, the traction staying continuous. This is the limit of a
#  layer of vanishing thickness `h → 0` whose out-of-plane compliance is
#  `Ŝ/h` and whose in-plane stiffness is zero. In the laminate formula,
#
#      f·ℙ            → ℙ_int / L        (finite)
#      f·ℂ_IO ℂ_OO⁻¹  → 0                (no in-plane/out-of-plane coupling)
#      f·Schur_IP(ℂ)  → 0                (no in-plane stiffness)
#
#  so a primal interface adds to `⟨ℙ⟩` and to nothing else.
#
#  ── Dual (membrane, surface-conductive): a surface stiffness ──────────────
#  The interfaces being PLANAR, `divₛ σˢ = 0`: there is no traction jump at
#  all. The surface stress is driven by the in-plane strain, which is
#  continuous and equal to the macroscopic `E`, so it adds straight to the
#  macroscopic stress:  `ℂ_hom ← ℂ_hom + Σ_j ℂˢ_j / L`, in the in-plane block.
#
#  Both terms carry the weight `1/L` — an interface *density*. Doubling every
#  thickness at fixed fractions therefore halves the interface correction, and
#  `L → ∞` recovers the perfect interface. That size effect is the physical
#  content of storing thicknesses rather than fractions.
# =============================================================================

# ── Elasticity, primal: out-of-plane compliance block ───────────────────────

"""
    _interface_P(itf, ::Type{T}) -> SMatrix{6,6,T}

Contribution of one interface to `⟨ℙ⟩` (before the `1/L` weight), in the
layer frame and in Kelvin-Mandel form. Non-zero for the primal (field-jump)
types only.

For [`SpringInterface`](@ref)`(kn, kt)` — whose fields are *compliances* —
the jump `[u] = 𝕂·(σ·n)` with `𝕂 = kn n⊗n + kt (δ − n⊗n)` contributes the
added strain `(𝕂·(σ·n)) ⊗ˢ n`, i.e. the out-of-plane block
`diag(kn, kt/2, kt/2)` in Mandel slots. The halving of the tangential term is
the symmetrised product, and it is produced by
`Core.compliance_op_block` — the very helper that turns `𝐊⁻¹` into `ℙ`, so
the "interface = zero-thickness layer" statement is literal in the code.
"""
_interface_P(::PerfectInterface, ::Type{T}) where {T} = zero(SMatrix{6, 6, T})

function _interface_P(itf::SpringInterface, ::Type{T}) where {T}
    kn = convert(T, itf.kn)
    kt = convert(T, itf.kt)
    z = zero(T)
    𝕂 = SMatrix{3, 3, T}(kt, z, z, z, kt, z, z, z, kn)   # (ℓ, m, n) frame
    return MFH_Core._op_embed(MFH_Core.compliance_op_block(𝕂))
end

# Dual types contribute nothing to ⟨ℙ⟩ …
_interface_P(::MembraneInterface, ::Type{T}) where {T} = zero(SMatrix{6, 6, T})
# … and the transport types are not elastic at all.
_interface_P(::KapitzaInterface, ::Type{T}) where {T} = zero(SMatrix{6, 6, T})
_interface_P(::SurfaceConductiveInterface, ::Type{T}) where {T} = zero(SMatrix{6, 6, T})

# ── Elasticity, dual: in-plane surface stiffness block ──────────────────────

"""
    _interface_Cs(itf, ::Type{T}) -> SMatrix{6,6,T}

Contribution of one interface to `ℂ_hom` (before the `1/L` weight), in the
layer frame and in Kelvin-Mandel form. Non-zero for the dual (surface
stiffness) types only.

For [`MembraneInterface`](@ref)`(κs, μs)` — Gurtin-Murdoch surface elasticity
with `κs = λs + μs` the surface dilatation modulus, matching the convention
of `LayeredSpheres` and of Echoes' `DUALDISC` — the 2-D surface law
`σˢ = λs tr(εˢ) p + 2μs εˢ` gives `C^s_1111 = κs + μs`, `C^s_1122 = κs − μs`,
`C^s_1212 = μs`, hence the in-plane Mandel block
`[κs+μs κs−μs 0; κs−μs κs+μs 0; 0 0 2μs]`.
"""
_interface_Cs(::PerfectInterface, ::Type{T}) where {T} = zero(SMatrix{6, 6, T})
_interface_Cs(::SpringInterface, ::Type{T}) where {T} = zero(SMatrix{6, 6, T})

function _interface_Cs(itf::MembraneInterface, ::Type{T}) where {T}
    κ = convert(T, itf.κs)
    μ = convert(T, itf.μs)
    z = zero(T)
    B = SMatrix{3, 3, T}(κ + μ, κ - μ, z, κ - μ, κ + μ, z, z, z, 2μ)
    return MFH_Core._ip_embed(B)
end

_interface_Cs(::KapitzaInterface, ::Type{T}) where {T} = zero(SMatrix{6, 6, T})
_interface_Cs(::SurfaceConductiveInterface, ::Type{T}) where {T} = zero(SMatrix{6, 6, T})

# ── Conduction, primal / dual ───────────────────────────────────────────────

"""
    _interface_P2(itf, ::Type{T}) -> SMatrix{3,3,T}

Order-2 analogue of [`_interface_P`](@ref): the contribution of one interface
to the out-of-plane "compliance" average of a transport problem.
[`KapitzaInterface`](@ref)`(ρ)` imposes `[T] = ρ q_n`, hence `ρ · n⊗n`.
"""
_interface_P2(::PerfectInterface, ::Type{T}) where {T} = zero(SMatrix{3, 3, T})

function _interface_P2(itf::KapitzaInterface, ::Type{T}) where {T}
    z = zero(T)
    ρ = convert(T, itf.resistance)
    return SMatrix{3, 3, T}(z, z, z, z, z, z, z, z, ρ)
end

_interface_P2(::SurfaceConductiveInterface, ::Type{T}) where {T} = zero(SMatrix{3, 3, T})
_interface_P2(::SpringInterface, ::Type{T}) where {T} = zero(SMatrix{3, 3, T})
_interface_P2(::MembraneInterface, ::Type{T}) where {T} = zero(SMatrix{3, 3, T})

"""
    _interface_Ks(itf, ::Type{T}) -> SMatrix{3,3,T}

Order-2 analogue of [`_interface_Cs`](@ref): a highly conductive 2-D surface
layer ([`SurfaceConductiveInterface`](@ref)`(ks)`) carries a surface flux
driven by the in-plane gradient, adding `ks (δ − n⊗n)` to the effective
conductivity.
"""
_interface_Ks(::PerfectInterface, ::Type{T}) where {T} = zero(SMatrix{3, 3, T})

function _interface_Ks(itf::SurfaceConductiveInterface, ::Type{T}) where {T}
    z = zero(T)
    ks = convert(T, itf.conductance)
    return SMatrix{3, 3, T}(ks, z, z, z, ks, z, z, z, z)
end

_interface_Ks(::KapitzaInterface, ::Type{T}) where {T} = zero(SMatrix{3, 3, T})
_interface_Ks(::SpringInterface, ::Type{T}) where {T} = zero(SMatrix{3, 3, T})
_interface_Ks(::MembraneInterface, ::Type{T}) where {T} = zero(SMatrix{3, 3, T})

# ── Assembly over the cell ──────────────────────────────────────────────────

"""
    _interface_terms(lam, ::Type{T}, ::Val{4}) -> (P_int, C_surf)
    _interface_terms(lam, ::Type{T}, ::Val{2}) -> (P_int, K_surf)

Sum the interface contributions of a laminate, each weighted by `1/L`.
Returns a pair `(primal, dual)` ready for
`Core.laminate_stiffness` / `Core.laminate_conductivity`.

Short-circuits to a pair of zeros when every interface is perfect, so the
common case pays nothing and is **bit-identical** to the no-interface path.
"""
function _interface_terms(lam::Laminate, ::Type{T}, ::Val{4}) where {T}
    all(itf -> itf isa PerfectInterface, lam.interfaces) &&
        return (zero(SMatrix{6, 6, T}), zero(SMatrix{6, 6, T}))
    # The interfaces carry their own element type: differentiating with respect
    # to an interface compliance makes `kn` a `Dual` while the moduli stay
    # `Float64`, so `T` alone would truncate the perturbation.
    Tp = _interface_eltype(lam, T)
    Z = zero(SMatrix{6, 6, Tp})
    L = laminate_period(lam)
    P = Z
    Cs = Z
    for itf in lam.interfaces
        P = P + _interface_P(itf, Tp)
        Cs = Cs + _interface_Cs(itf, Tp)
    end
    return (P / L, Cs / L)
end

function _interface_terms(lam::Laminate, ::Type{T}, ::Val{2}) where {T}
    all(itf -> itf isa PerfectInterface, lam.interfaces) &&
        return (zero(SMatrix{3, 3, T}), zero(SMatrix{3, 3, T}))
    Tp = _interface_eltype(lam, T)
    Z = zero(SMatrix{3, 3, Tp})
    L = laminate_period(lam)
    P = Z
    Ks = Z
    for itf in lam.interfaces
        P = P + _interface_P2(itf, Tp)
        Ks = Ks + _interface_Ks(itf, Tp)
    end
    return (P / L, Ks / L)
end

"""
    _interface_eltype(lam, ::Type{T}) -> Type

Promotion of `T` (the element type of the layers and thicknesses) with the
element type of every interface, so that a `ForwardDiff.Dual` reaching an
interface compliance alone still propagates.
"""
function _interface_eltype(lam::Laminate, ::Type{T}) where {T}
    Tp = T
    for itf in lam.interfaces
        Tp = promote_type(Tp, eltype(itf))
    end
    return Tp
end
