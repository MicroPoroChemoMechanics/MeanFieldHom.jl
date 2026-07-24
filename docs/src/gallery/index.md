# Gallery

The [Tutorials](../tutorials/index.md) and [Applications](../applications/transport.md)
sections are hand-written narratives. This Gallery instead auto-generates a
doc page from a `scripts/` demo via [Literate.jl](https://github.com/fredrikekre/Literate.jl)
— **only** for scripts that cover a topic no tutorial or application already
does. Where a script duplicates an existing tutorial, it stays a plain
script (see `scripts/README.md`), avoiding two pages that say the same
thing.

Each script listed here also ships as a self-contained Jupyter notebook and
a cleaned standalone `.jl` file, generated alongside this page — see
`scripts/README.md` for how to fetch them from a local build
(`julia --project=docs docs/literate.jl`).

## Pages

- [n-layer sphere: volume-averaged localization tensors](generated/30_average_nlayers.md)
  — `LayeredSphere` / Hervé-Zaoui localization, no tutorial covers this.
- [Symmetrization showcase](generated/70_symmetrization_showcase.md) — exact
  rotation-group average vs. best-fit projection on a non-major-symmetric
  tensor.
- [n-layer confocal spheroid: imperfect-interface conductivity](generated/32_spheroid_nlayers_conductivity.md)
  — Mori-Tanaka effective conductivity vs. Kapitza interface parameter, the
  [Kushch, Sevostianov & Belyaev (2015)](@cite kushch2015) setting.
- [n-layer confocal spheroid: series truncation](generated/33_spheroid_series_convergence.md)
  — harmonic-series convergence vs. `𝒩`, and why the default `QuadGK`
  backend stays accurate in `Float64` where the original monomial series
  needs BigFloat.
- [n-layer confocal spheroid: equivalent conductivity](generated/34_spheroid_equivalent_conductivity.md)
  — the exact equivalent-particle conductivity `kᵉ𝑞 = ⟨B⟩·⟨A⟩⁻¹`
  ([Barthélémy & Bignonnet 2020](@cite barthelemyBignonnetIJES2020), §4)
  vs. aspect ratio.
- [n-layer confocal spheroid: local fields](generated/35_spheroid_local_fields.md)
  — pointwise temperature and flux maps across a 2-layer spheroid with an
  imperfect interface.
