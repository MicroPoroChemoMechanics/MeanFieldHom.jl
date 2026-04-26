# =============================================================================
#  trapezoidal.jl — discrete representation of a Stieltjes / Volterra
#  integral on a time grid via the trapezoidal rule of Sanahuja 2013.
#
#  Given a kernel `f(t, t')` (scalar, tensor, or 6×6 matrix), build the
#  lower-block-triangular matrix `M` of size `(B·n) × (B·n)` such that
#  the discrete relation `y = M · x` approximates the Stieltjes integral
#       y(t_i) = ∫_{t_0}^{t_i} f(t_i, τ) dx(τ)
#  on the time grid `t_0 < t_1 < … < t_{n-1}` (ECHOES `visco_law.cpp:80–103`).
#
#  Block-row index `i = 1..n` corresponds to time `t_{i-1}` (1-based).
#  Block-column index `j = 1..n` to the trapezoidal weight at time
#  `t_{j-1}`.  The block matrix layout (1-based) is:
#       M[1, 1]   = f(t_0, t_0)
#       M[i, i]   = 0.5 (f(t_{i-1}, t_{i-2}) + f(t_{i-1}, t_{i-1}))   (i ≥ 2)
#       M[i, 1]   = 0.5 (f(t_{i-1}, t_0)     - f(t_{i-1}, t_1))       (i ≥ 2)
#       M[i, j]   = 0.5 (f(t_{i-1}, t_{j-2}) - f(t_{i-1}, t_j))       (1 < j < i)
#       M[i, j]   = 0                                                  (j > i)
# =============================================================================

"""
    trapezoidal_matrix(law::ViscoLaw, times::AbstractVector{<:Real}) -> Matrix

Build the discrete `(B·n) × (B·n)` lower-block-triangular matrix
representing the Stieltjes integral
``y(t_i) = ∫_{t_0}^{t_i} f(t_i, τ) \\, dx(τ)``
on the time grid `times = (t_0, …, t_{n-1})` using the trapezoidal
rule of [@sanahuja2013].  `B = 1` for scalar-valued kernels and `B = 6`
for 4-tensor- or `6×6`-matrix-valued kernels (Mandel form).

The exact block layout is described at the top of `trapezoidal.jl`.
The output is type-stable in the element type returned by
`visco_eval(law, ...)`.
"""
function trapezoidal_matrix(law::ViscoLaw, times::AbstractVector{<:Real})
    n = length(times)
    n ≥ 1 || throw(ArgumentError("trapezoidal_matrix: times must be non-empty"))
    sample = visco_eval(law, times[1], times[1])
    return _trapezoidal_dispatch(law, times, sample)
end

# Scalar specialization: M is `n × n`.
function _trapezoidal_dispatch(law::ViscoLaw, times::AbstractVector{<:Real},
                               ::Number)
    n = length(times)
    T = typeof(visco_eval(law, times[1], times[1]))
    M = zeros(T, n, n)
    _fill_trapezoidal_scalar!(M, law, times)
    return M
end

# Tensor specialization (TensND.AbstractTens{4,3}): M is `6n × 6n`.
function _trapezoidal_dispatch(law::ViscoLaw, times::AbstractVector{<:Real},
                               ::TensND.AbstractTens{4, 3})
    n = length(times)
    sample = visco_eval(law, times[1], times[1])
    T = eltype(TensND.get_array(sample))
    M = zeros(T, 6 * n, 6 * n)
    _fill_trapezoidal_tensor!(M, law, times)
    return M
end

# Direct 6×6 Mandel matrix specialization.
function _trapezoidal_dispatch(law::ViscoLaw, times::AbstractVector{<:Real},
                               sample::AbstractMatrix)
    size(sample) == (6, 6) ||
        throw(ArgumentError("trapezoidal_matrix: only 6×6 matrices are supported"))
    n = length(times)
    T = eltype(sample)
    M = zeros(T, 6 * n, 6 * n)
    _fill_trapezoidal_mandel!(M, law, times)
    return M
end

# ── Scalar case ─────────────────────────────────────────────────────────────

@inline function _fill_trapezoidal_scalar!(M::AbstractMatrix, law::ViscoLaw,
                                           times::AbstractVector)
    n = length(times)
    @inbounds begin
        # Diagonal entry M[1,1] = f(t_0, t_0).
        M[1, 1] = visco_eval(law, times[1], times[1])
        for i in 2:n
            M[i, i] = (visco_eval(law, times[i], times[i - 1])
                       + visco_eval(law, times[i], times[i])) / 2
            M[i, 1] = (visco_eval(law, times[i], times[1])
                       - visco_eval(law, times[i], times[2])) / 2
            for j in 2:(i - 1)
                M[i, j] = (visco_eval(law, times[i], times[j - 1])
                           - visco_eval(law, times[i], times[j + 1])) / 2
            end
        end
    end
    return M
end

# ── Tensor / Mandel case (block-by-block) ───────────────────────────────────

# Convert a sample (4-tensor or 6×6 matrix) to a 6×6 Mandel array,
# preserving element type.
@inline _to_mandel(sample::AbstractMatrix) = sample
@inline function _to_mandel(sample::TensND.AbstractTens{4, 3})
    return _tens_to_mandel66(sample)
end

# Convert a TensND 4-tensor in {3,3,3,3} array form to a 6×6 Mandel matrix.
# Mandel convention: σ_4 = √2 σ_23, σ_5 = √2 σ_13, σ_6 = √2 σ_12.
function _tens_to_mandel66(C::TensND.AbstractTens{4, 3})
    arr = TensND.get_array(C)
    T = eltype(arr)
    M = zeros(T, 6, 6)
    sq2 = sqrt(T(2))
    voigt = ((1, 1), (2, 2), (3, 3), (2, 3), (1, 3), (1, 2))
    for I in 1:6, J in 1:6
        i, j = voigt[I]
        k, l = voigt[J]
        scale = (I ≥ 4 ? sq2 : one(T)) * (J ≥ 4 ? sq2 : one(T))
        M[I, J] = arr[i, j, k, l] * scale
    end
    return M
end

# Convert a 6×6 Mandel matrix back to a {3,3,3,3} array.  Used downstream
# (not strictly needed in `trapezoidal_matrix` itself).
function _mandel66_to_tens(M::AbstractMatrix)
    T = eltype(M)
    arr = zeros(T, 3, 3, 3, 3)
    sq2_inv = one(T) / sqrt(T(2))
    voigt = ((1, 1), (2, 2), (3, 3), (2, 3), (1, 3), (1, 2))
    for I in 1:6, J in 1:6
        i, j = voigt[I]
        k, l = voigt[J]
        scale = (I ≥ 4 ? sq2_inv : one(T)) * (J ≥ 4 ? sq2_inv : one(T))
        v = M[I, J] * scale
        arr[i, j, k, l] = v
        arr[j, i, k, l] = v
        arr[i, j, l, k] = v
        arr[j, i, l, k] = v
    end
    return arr
end

# Fill the (6×6)-block-structured matrix from a 4-tensor kernel.
@inline function _fill_trapezoidal_tensor!(M::AbstractMatrix, law::ViscoLaw,
                                           times::AbstractVector)
    n = length(times)
    @inbounds for i in 1:n
        for j in 1:i
            block = _block_value_tensor(law, times, i, j)
            _set_block!(M, i, j, block)
        end
    end
    return M
end

@inline function _fill_trapezoidal_mandel!(M::AbstractMatrix, law::ViscoLaw,
                                           times::AbstractVector)
    n = length(times)
    @inbounds for i in 1:n
        for j in 1:i
            block = _block_value_mandel(law, times, i, j)
            _set_block!(M, i, j, block)
        end
    end
    return M
end

# Compute the (i, j) trapezoidal weight as a 6×6 Mandel block from a
# 4-tensor kernel.  Indices are 1-based: `i` = 1..n, `j` = 1..i.
@inline function _block_value_tensor(law::ViscoLaw, times::AbstractVector, i::Int, j::Int)
    if i == 1 && j == 1
        return _to_mandel(visco_eval(law, times[1], times[1]))
    end
    if j == i
        a = _to_mandel(visco_eval(law, times[i], times[i - 1]))
        b = _to_mandel(visco_eval(law, times[i], times[i]))
        return (a .+ b) ./ 2
    end
    if j == 1
        a = _to_mandel(visco_eval(law, times[i], times[1]))
        b = _to_mandel(visco_eval(law, times[i], times[2]))
        return (a .- b) ./ 2
    end
    # Interior 1 < j < i.
    a = _to_mandel(visco_eval(law, times[i], times[j - 1]))
    b = _to_mandel(visco_eval(law, times[i], times[j + 1]))
    return (a .- b) ./ 2
end

@inline function _block_value_mandel(law::ViscoLaw, times::AbstractVector, i::Int, j::Int)
    if i == 1 && j == 1
        return copy(visco_eval(law, times[1], times[1]))
    end
    if j == i
        a = visco_eval(law, times[i], times[i - 1])
        b = visco_eval(law, times[i], times[i])
        return (a .+ b) ./ 2
    end
    if j == 1
        a = visco_eval(law, times[i], times[1])
        b = visco_eval(law, times[i], times[2])
        return (a .- b) ./ 2
    end
    a = visco_eval(law, times[i], times[j - 1])
    b = visco_eval(law, times[i], times[j + 1])
    return (a .- b) ./ 2
end

# Place a 6×6 block into the (i, j)-th block of the 6n×6n matrix.
@inline function _set_block!(M::AbstractMatrix, i::Int, j::Int, block::AbstractMatrix)
    rows = (6 * (i - 1) + 1):(6 * i)
    cols = (6 * (j - 1) + 1):(6 * j)
    @inbounds M[rows, cols] = block
    return M
end

# Read the (i, j)-th 6×6 block.
@inline function _get_block(M::AbstractMatrix, i::Int, j::Int)
    rows = (6 * (i - 1) + 1):(6 * i)
    cols = (6 * (j - 1) + 1):(6 * j)
    return @inbounds M[rows, cols]
end
