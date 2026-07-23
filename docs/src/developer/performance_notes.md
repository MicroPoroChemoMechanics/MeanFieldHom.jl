# Performance notes

- All hot-path helpers (`_inv3`, `_acoustic_tensor`, `_qnn_pair_components!`)
  use in-place `Matrix{T}` buffers and avoid LU factorization on 3×3
  matrices.
- `TensISO` / `TensTI{4}` / `TensOrtho` specializations return the
  *most specific* TensND type — this propagates to downstream
  homogenization schemes and avoids redundant symmetry checks.
- `ForwardDiff.Dual` propagation is honoured through nested `QuadGK`
  (no `PolynomialRoots` in the AD path).
