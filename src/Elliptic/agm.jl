# ═════════════════════════════════════════════════════════════════════════════
#  Complete integrals via AGM (Abramowitz & Stegun 17.6, NIST DLMF 19.8)
# ═════════════════════════════════════════════════════════════════════════════

function _ell_K_agm(m::T) where {T <: Number}
    a = one(T)
    b = sqrt(one(T) - m)
    tol = _agm_tol(T)
    for _ in 1:60
        a_new = (a + b) / 2
        b_new = sqrt(a * b)
        a, b = a_new, b_new
        _agm_converged(T, a, b, tol) && break
    end
    return T(π) / (2 * a)
end

function _ell_E_agm(m::T) where {T <: Number}
    a = one(T)
    b = sqrt(one(T) - m)
    c = sqrt(m)
    s = c^2 / 2                              # (1/2) · 2⁰ · c₀²
    p = one(T)                               # running 2ⁿ
    tol = _agm_tol(T)
    for _ in 1:60
        a_new = (a + b) / 2
        b_new = sqrt(a * b)
        c_new = (a - b) / 2
        p *= 2
        s += p * c_new^2 / 2
        a, b, c = a_new, b_new, c_new
        _agm_converged(T, a, b, tol) && break
    end
    K_val = T(π) / (2 * a)
    return K_val * (one(T) - s)
end

# ─── Tolerances — tuned per scalar type ──────────────────────────────────────

function _agm_tol(::Type{T}) where {T <: Number}
    if T <: AbstractFloat
        return 10 * eps(T)
    elseif T <: Real
        return 10 * eps(Float64)
    else
        return 0.0                 # symbolic types — fall back to fixed max-iter
    end
end

function _agm_converged(::Type{T}, a, b, tol) where {T <: Number}
    T <: Real || return false
    denom = abs(a)
    return iszero(denom) ? abs(a - b) < tol : abs(a - b) < tol * denom
end
