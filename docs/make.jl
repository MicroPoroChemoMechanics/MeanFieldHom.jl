using Documenter
using DocumenterCitations
using MeanFieldHom

# GR needs a headless display driver on CI runners; without this the figures in
# the Applications pages fail to render.
ENV["GKSwstype"] = "100"

# Generates the Gallery pages (+ companion notebooks/scripts) from the
# curated `scripts/` demos before `makedocs` runs, so the generated markdown
# exists when `pages` below references it. Must run after the GKSwstype
# assignment above — the Literate notebook pass actually executes the
# scripts, including their Plots/GR calls.
include("literate.jl")

bib = CitationBibliography(
    joinpath(@__DIR__, "src", "references.bib");
    style = :numeric,
)

DocMeta.setdocmeta!(
    MeanFieldHom,
    :DocTestSetup,
    :(using MeanFieldHom);
    recursive = true,
)

makedocs(;
    clean    = false,
    modules  = [MeanFieldHom,
                MeanFieldHom.Elliptic,
                MeanFieldHom.Core,
                MeanFieldHom.Elasticity,
                MeanFieldHom.Cracks,
                MeanFieldHom.Conductivity,
                MeanFieldHom.LayeredSpheres,
                MeanFieldHom.LayeredSpheroids,
                MeanFieldHom.Schemes,
                MeanFieldHom.Viscoelasticity],
    remotes  = nothing,
    authors  = "Jean-François Barthélémy",
    sitename = "MeanFieldHom.jl",
    format   = Documenter.HTML(;
        canonical        = "https://MicroPoroChemoMechanics.github.io/MeanFieldHom.jl",
        repolink         = "https://github.com/MicroPoroChemoMechanics/MeanFieldHom.jl",
        edit_link        = "main",
        assets           = ["assets/favicon.ico", "assets/custom.css"],
        prettyurls       = (get(ENV, "CI", nothing) == "true"),
        collapselevel    = 1,
        mathengine       = Documenter.MathJax3(),
        # The interactive Plotly 3D percolation surfaces in the cement-paste
        # diffusion chapter embed their data inline, exceeding the 200 KiB
        # default; raise the ceiling for those pages.
        size_threshold        = 3_000_000,
        size_threshold_warn   = 1_500_000,
        # The interactive 3D surfaces embed their data as inline HTML; allow it.
        example_size_threshold = 2_000_000,
    ),
    plugins = [bib],
    pages = [
        "Home" => "index.md",
        "Theory"  => [
            "theory/overview.md",
            "theory/hill_tensors.md",
            "theory/cod_tensors.md",
            "theory/thermal_cracks.md",
            "theory/localization.md",
            "theory/homogenization.md",
            "theory/layered_sphere.md",
            "theory/layered_spheroid.md",
            "theory/viscoelasticity.md",
            "theory/elliptic_integrals.md",
            "theory/cylindrical_limits.md",
        ],
        "Manual"  => [
            "manual/installation.md",
            "manual/ellipsoidal_inclusions.md",
            "manual/cylindrical_inclusions.md",
            "manual/cracks.md",
            "manual/conductivity.md",
            "manual/schemes.md",
            "manual/viscoelasticity.md",
            "manual/sensitivities.md",
            "manual/elliptic_examples.md",
        ],
        "Tutorials" => [
            "tutorials/index.md",
            "tutorials/01_first_estimate.md",
            "tutorials/02_bounds_and_schemes.md",
            "tutorials/03_porous_materials.md",
            "tutorials/04_porous_benchmark.md",
            "tutorials/05_differential_paths.md",
            "tutorials/06_cracks.md",
            "tutorials/07_viscoelasticity.md",
            "tutorials/08_sensitivities.md",
            "tutorials/09_strength_criteria.md",
            "tutorials/10_from_echoes.md",
            "tutorials/11_symbolic_spheres.md",
            "tutorials/12_nonlinear_solvers.md",
            "tutorials/13_layered_spheroid.md",
        ],
        "Applications" => [
            "applications/transport.md",
            "applications/cement_paste.md",
            "applications/cement_paste_diffusion.md",
            "applications/strength.md",
            "applications/bituminous.md",
            "applications/ageing_creep.md",
        ],
        "Gallery" => [
            "gallery/index.md",
            "gallery/generated/30_average_nlayers.md",
            "gallery/generated/70_symmetrization_showcase.md",
            "gallery/generated/32_spheroid_nlayers_conductivity.md",
            "gallery/generated/33_spheroid_series_convergence.md",
            "gallery/generated/34_spheroid_equivalent_conductivity.md",
            "gallery/generated/35_spheroid_local_fields.md",
        ],
        "Developer" => [
            "developer/architecture.md",
            "developer/adding_inclusion.md",
            "developer/adding_algorithm.md",
            "developer/adding_scheme.md",
            "developer/testing_conventions.md",
            "developer/performance_notes.md",
            "developer/roadmap.md",
        ],
        "API" => [
            "api/elliptic.md",
            "api/core.md",
            "api/elasticity.md",
            "api/cracks.md",
            "api/conductivity.md",
            "api/localization.md",
            "api/layered_sphere.md",
            "api/layered_spheroid.md",
            "api/schemes.md",
            "api/viscoelasticity.md",
            "api/sensitivities.md",
        ],
        "References" => "references.md",
    ],
    warnonly = true,
)

deploydocs(;
    repo         = "github.com/MicroPoroChemoMechanics/MeanFieldHom.jl.git",
    devbranch    = "main",
    push_preview = false,
)
