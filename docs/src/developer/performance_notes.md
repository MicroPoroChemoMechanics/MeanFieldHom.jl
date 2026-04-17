# Performance notes

- All hot-path helpers (`_inv3`, `_acoustic_tensor`, `_qnn_pair_components!`)
  use in-place `Matrix{T}` buffers and avoid LU factorisation on 3×3
  matrices.
- `TensISO` / `TensWalpole` / `TensOrtho` specialisations return the
  *most specific* TensND type — this propagates to downstream
  homogenisation schemes and avoids redundant symmetry checks.
- `ForwardDiff.Dual` propagation is honoured through nested `QuadGK`
  (no `PolynomialRoots` in the AD path).
