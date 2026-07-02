using Documenter
using DocumenterCitations
using MeanFieldHom

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
        assets           = ["assets/favicon.ico", "assets/custom.css"],
        prettyurls       = (get(ENV, "CI", nothing) == "true"),
        collapselevel    = 1,
        mathengine       = Documenter.MathJax3(),
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
