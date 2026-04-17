# Adding a homogenisation scheme

`MeanFieldHom.Schemes` is currently a placeholder.  When a scheme is
ready:

1. Create `src/Schemes/<scheme_name>.jl`, include it from
   `src/Schemes/Schemes.jl`.
2. Consume `hill_tensor` / `cod_tensor` / `compliance_contribution` for
   the per-phase / per-crack ingredients.
3. Export the public API through `src/Schemes/Schemes.jl` and re-export
   from `src/MeanFieldHom.jl`.
4. Add a unit test under `test/Schemes/` and a manual chapter under
   `docs/src/manual/<scheme_name>.md`.
