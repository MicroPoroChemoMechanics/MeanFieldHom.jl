# =============================================================================
#  legendre.jl — associated Legendre functions of order m = 0, 1, first
#  (P) and second (Q) kind, by upward three-term recurrence.
#
#  Faithful transliteration of the seed values and recurrences used in
#  Barthélémy & Bignonnet (IJES 2020, eq:Legformc/d, eq:Leg2formc/d).
#
#  Two arguments are used throughout the spheroid solution:
#    - the "p" argument, real, |p| ≤ 1 (the angular / polar coordinate);
#    - the "q" argument, |q| > 1 for prolate (real) or purely imaginary
#      `q = iτ` for oblate (`τ` real) — carried generically as `Q<:Number`
#      (`T` for prolate, `Complex{T}` for oblate). No branch is taken on
#      the caller's side: `sqrt`, `atanh` dispatch to the right method
#      once `x` has the right static type (`Complex` for oblate).
#
#  Only ODD degrees 1, 3, …, 2𝒩−1 ever enter the axial/transverse
#  spheroid problem (by symmetry, see eq:Taxi/eq:Ttrans); the table
#  builders below return exactly those `𝒩` values (and derivatives),
#  never materializing the unused even-degree entries.
# =============================================================================

@inline _arccoth(x) = atanh(one(x) / x)

@inline function _legendre_next(n, m, x, Pnm, Pnm1)
    return ((2n + 1) * x * Pnm - (n + m) * Pnm1) / (n - m + 1)
end

@inline function _legendre_next_der(n, m, x, Pnm, dPnm, dPnm1)
    return ((2n + 1) * (Pnm + x * dPnm) - (n + m) * dPnm1) / (n - m + 1)
end

"""
    _legendre_grow!(tab, dtab, Nmax, m, x)

Grow the value/derivative tables `tab`, `dtab` (1-indexed, `tab[k+1]` =
degree-`k` value) up to degree `Nmax` (inclusive) by the upward
recurrence [`_legendre_next`](@ref) / [`_legendre_next_der`](@ref),
order `m`. `tab` and `dtab` must already hold their required seed
degrees (2 seeds for the standard case, 3 for [`_Q1_table`](@ref)'s
special low-degree closed forms).
"""
function _legendre_grow!(tab::Vector{Tx}, dtab::Vector{Tx}, Nmax::Int, m::Int, x) where {Tx}
    n = length(tab) - 1
    while n < Nmax
        Pn = tab[n + 1]
        Pnm1 = tab[n]
        Pnext = _legendre_next(n, m, x, Pn, Pnm1)
        push!(tab, Pnext)
        if length(dtab) == length(tab) - 1
            dPn = dtab[n + 1]
            dPnm1 = dtab[n]
            push!(dtab, _legendre_next_der(n, m, x, Pn, dPn, dPnm1))
        end
        n += 1
    end
    return tab, dtab
end

"""
    _P0_table(x, Nmax) -> (tab, dtab)

`Pₙ(x)`, `n = 0, …, Nmax`, plain Legendre polynomials (order `m = 0`).
Valid for any argument (the `p` branch, `|p| ≤ 1`, or the `q` branch).
"""
function _P0_table(x::Tx, Nmax::Int) where {Tx}
    tab = Tx[one(Tx), x]
    dtab = Tx[zero(Tx), one(Tx)]
    return _legendre_grow!(tab, dtab, Nmax, 0, x)
end

"""
    _Q0_table(x, Nmax) -> (tab, dtab)

`Qₙ(x)`, `n = 0, …, Nmax`, Legendre functions of the second kind
(order `m = 0`). Valid for the `q` branch (`|x| > 1`, real or the
oblate `iτ` substitute).
"""
function _Q0_table(x::Tx, Nmax::Int) where {Tx}
    ax = _arccoth(x)
    x2m1 = x^2 - one(Tx)
    tab = Tx[ax, x * ax - one(Tx)]
    dtab = Tx[-one(Tx) / x2m1, ax - x / x2m1]
    return _legendre_grow!(tab, dtab, Nmax, 0, x)
end

"""
    _P1p_table(x, Nmax) -> (tab, dtab)

`Pₙ¹(x)`, `n = 0, …, Nmax`, associated Legendre of the first kind,
order `m = 1`, on the `p` branch (`|p| ≤ 1`), seeded with
`P₁¹(p) = -√(1-p²)`.
"""
function _P1p_table(x::Tx, Nmax::Int) where {Tx}
    xb = -sqrt(one(Tx) - x^2)
    tab = Tx[zero(Tx), xb]
    dtab = Tx[zero(Tx), -x / xb]
    return _legendre_grow!(tab, dtab, Nmax, 1, x)
end

"""
    _P1_table(x, Nmax) -> (tab, dtab)

`Pₙ¹(x)`, `n = 0, …, Nmax`, associated Legendre of the first kind,
order `m = 1`, on the `q` branch (`|x| > 1`), seeded with
`P₁¹(q) = √(q²-1)`.
"""
function _P1_table(x::Tx, Nmax::Int) where {Tx}
    xb = sqrt(x^2 - one(Tx))
    tab = Tx[zero(Tx), xb]
    dtab = Tx[zero(Tx), x / xb]
    return _legendre_grow!(tab, dtab, Nmax, 1, x)
end

"""
    _Q1_table(x, Nmax) -> (tab, dtab)

`Qₙ¹(x)`, `n = 0, …, Nmax`, associated Legendre of the second kind,
order `m = 1`, on the `q` branch. The recurrence for `m = 1` is
singular at `n = 0`, so the degrees `0, 1, 2` are seeded from closed
forms and the upward recurrence resumes from `n = 2`.
"""
function _Q1_table(x::Tx, Nmax::Int) where {Tx}
    ax = _arccoth(x)
    xb = sqrt(x^2 - one(Tx))
    x2 = x^2
    x2m1 = x2 - one(Tx)
    tab = Tx[
        zero(Tx),
        xb * ax - x / xb,
        x * xb * (3 * ax - (3 * x2 - 2) / (x * x2m1)),
    ]
    dtab = Tx[
        zero(Tx),
        x / xb * (ax + (2 - x2) / (x * x2m1)),
        (2 * x2 - 1) / xb * (3 * ax - x * (6 * x2 - 7) / ((2 * x2 - 1) * x2m1)),
    ]
    return _legendre_grow!(tab, dtab, Nmax, 1, x)
end

"""
    legendre_odd(kind::Symbol, x, Nseries::Int) -> (vals, derivs)

Values and derivatives of the requested Legendre kind at the `Nseries`
ODD degrees `1, 3, …, 2·Nseries − 1`, as length-`Nseries` `Vector`s
(index `r` ↔ degree `2r − 1`).

`kind ∈ (:P0, :Q0, :P1, :P1p, :Q1)`:
- `:P0`  — `Pₙ(x)`   (m=0, any branch)
- `:Q0`  — `Qₙ(x)`   (m=0, q branch, |x|>1)
- `:P1`  — `Pₙ¹(x)`  (m=1, q branch, |x|>1)
- `:P1p` — `Pₙ¹(x)`  (m=1, p branch, |x|≤1)
- `:Q1`  — `Qₙ¹(x)`  (m=1, q branch, |x|>1)
"""
function legendre_odd(kind::Symbol, x, Nseries::Int)
    Nmax = 2 * Nseries - 1
    tab, dtab = if kind === :P0
        _P0_table(x, Nmax)
    elseif kind === :Q0
        _Q0_table(x, Nmax)
    elseif kind === :P1
        _P1_table(x, Nmax)
    elseif kind === :P1p
        _P1p_table(x, Nmax)
    elseif kind === :Q1
        _Q1_table(x, Nmax)
    else
        throw(ArgumentError("legendre_odd: unknown kind $kind"))
    end
    idx = 2:2:(2 * Nseries)
    return tab[idx], dtab[idx]
end
