using Documenter
using Documenter.Remotes
using MeanFieldHom

DocMeta.setdocmeta!(
    MeanFieldHom,
    :DocTestSetup,
    :(using MeanFieldHom);
    recursive = true,
)

makedocs(;
    clean    = false,
    modules  = [MeanFieldHom,
                MeanFieldHom.Core,
                MeanFieldHom.Elasticity,
                MeanFieldHom.Cracks,
                MeanFieldHom.Conductivity,
                MeanFieldHom.Schemes],
    authors  = "Jean-François Barthélémy",
    sitename = "MeanFieldHom.jl",
    remotes  = Dict(
        joinpath(@__DIR__, "..") => (Remotes.GitHub("MicMacTools", "MeanFieldHom.jl"), "main"),
    ),
    format   = Documenter.HTML(;
        canonical        = "https://MicMacTools.github.io/MeanFieldHom.jl",
        edit_link        = "main",
        assets           = ["assets/favicon.ico", "assets/custom.css"],
        prettyurls       = (get(ENV, "CI", nothing) == "true"),
        collapselevel    = 1,
        mathengine       = Documenter.MathJax3(),
    ),
    pages = [
        "Home" => "index.md",
        "Theory"  => [
            "theory/overview.md",
        ],
        "Manual"  => [
            "manual/installation.md",
            "manual/ellipsoidal_inclusions.md",
            "manual/cracks.md",
            "manual/conductivity.md",
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
            "api/core.md",
            "api/elasticity.md",
            "api/cracks.md",
            "api/conductivity.md",
        ],
    ],
    warnonly = true,
)

deploydocs(;
    repo      = "github.com/MicMacTools/MeanFieldHom.jl.git",
    devbranch = "main",
)
