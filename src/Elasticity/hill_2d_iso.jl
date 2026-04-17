# =============================================================================
#  hill_2d_iso.jl — Hill tensor for a 2-D ellipse in an isotropic matrix.
# =============================================================================

function _hill_2d_iso(ell::Ellipsoid{2, Circular}, C₀)
    T = promote_type(eltype(ell.semi_axes), eltype(C₀))
    α, β = C₀.data
    k = α / 3
    μ = β / 2
    local P1111, P1122
    if isa(k, AbstractFloat) && isinf(k)
        den   = 8μ
        P1111 = one(T) / den
        P1122 = -one(T) / den
    else
        den   = 8μ * (3k + 4μ)
        P1111 = (3k + 13μ) / den
        P1122 = -(3k + μ) / den
    end
    return TensND.TensISO{2}(T(P1111 + P1122), T(P1111 - P1122))
end

function _hill_2d_iso(ell::Ellipsoid{2, Elliptic}, C₀)
    T  = promote_type(eltype(ell.semi_axes), eltype(C₀))
    α, β = C₀.data
    k  = α / 3
    μ  = β / 2
    ρ  = T(ell.semi_axes[2] / ell.semi_axes[1])
    ρ2 = ρ * ρ
    up2 = (1 + ρ)^2
    local P1111, P2222, P1122, P1212
    if isa(k, AbstractFloat) && isinf(k)
        den   = 2μ * up2
        P1111 = ρ / den
        P1122 = -ρ / den
        P2222 = ρ / den
        P1212 = (ρ2 + 1) / (2 * den)
    else
        den   = 2μ * (3k + 4μ) * up2
        P1111 = ρ * (3k + 7μ + 6ρ*μ) / den
        P1122 = -ρ * (3k + μ) / den
        P2222 = (6μ + 7ρ*μ + 3ρ*k) / den
        P1212 = (3k + 4μ + 3ρ2*k + 4ρ2*μ + 6ρ*μ) / (2 * den)
    end
    P_arr = zeros(T, 2, 2, 2, 2)
    P_arr[1,1,1,1] = P1111
    P_arr[2,2,2,2] = P2222
    P_arr[1,1,2,2] = P_arr[2,2,1,1] = P1122
    for (a,b,c,d) in ((1,2,1,2),(1,2,2,1),(2,1,1,2),(2,1,2,1))
        P_arr[a,b,c,d] = P1212
    end
    return TensND.change_tens_canon(TensND.Tens(P_arr, ell.basis))
end
