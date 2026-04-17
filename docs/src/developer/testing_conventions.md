# Testing conventions

- Unit tests live under `test/<SubModule>/`.
- Broader regression tests live under `test/regression/` and exercise
  representative Hill-tensor and crack cases across the full dispatch
  surface.
- `test/runtests.jl` aggregates all testsets.
- Every new public function must ship with at least one smoke test.
