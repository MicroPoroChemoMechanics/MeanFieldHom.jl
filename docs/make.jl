using Documenter
using DocumenterCitations
using MeanFieldHom

# GR needs a headless display driver on CI runners; without this the figures in
# the Applications pages fail to render.
ENV["GKSwstype"] = "100"

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
                MeanFieldHom.Schemes,
                MeanFieldHom.Viscoelasticity],
    remotes  = nothing,
    authors  = "Jean-François Barthélémy",
    sitename = "MeanFieldHom.jl",
    format   = Documenter.HTML(;
        canonical        = "https://MicroPoroChemoMechanics.github.io/MeanFieldHom.jl",
        repolink         = "https://github.com/MicroPoroChemoMechanics/MeanFieldHom.jl",
        edit_link        = "main",
        assets           = [
            "assets/favicon.ico",
            "assets/custom.css",
            # plotly.js for the interactive 3D percolation surfaces in the
            # cement-paste diffusion chapter (loaded globally so the inline
            # `Plotly.newPlot` divs render without requirejs/WebIO).
            Documenter.asset("https://cdn.plot.ly/plotly-2.35.2.min.js"; class = :js, islocal = false),
        ],
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
        "Applications" => [
            "applications/transport.md",
            "applications/cement_paste.md",
            "applications/cement_paste_diffusion.md",
            "applications/strength.md",
            "applications/bituminous.md",
            "applications/ageing_creep.md",
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
