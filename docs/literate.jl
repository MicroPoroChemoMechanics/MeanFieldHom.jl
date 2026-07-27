# docs/literate.jl
#
# Runs Literate.jl over the `scripts/` demos that are published as tutorial
# pages. Produces, per script, three artifacts:
#
#   - a Documenter-ready markdown page  -> docs/src/tutorials/generated/
#     (executed by Documenter itself via `@example` when `makedocs` runs —
#     Literate's `markdown(...; documenter=true)` leaves `execute=false`)
#   - a pre-run Jupyter notebook        -> docs/generated_notebooks/
#   - a cleaned standalone .jl script   -> docs/generated_scripts/
#     (markup-only lines stripped, `#jl` directives resolved)
#
# Called from `docs/make.jl` *before* `makedocs`, so the generated markdown
# exists when Documenter's `pages` list references it.
#
# There is no longer a separate "Gallery" section: a page generated from a
# script and a page written by hand are both just tutorials, and the reader has
# no reason to care which is which. What decides whether a script is published
# is the same as before — does it cover a topic no other tutorial does — with
# the classification kept in `Assets/plans/MFH_LITERATE_SCRIPTS.md`. Scripts
# that duplicate an existing tutorial stay plain scripts, and SymPy-heavy ones
# are never published (they are re-executed on every docs build).
#
# `PUBLISHED_SCRIPTS` maps each script to its **page name**. Scripts keep their
# numeric prefixes (they encode a running order), while pages carry thematic
# names, so inserting a tutorial never forces a renumbering. The mapping is what
# Literate's `name` option is for.

using Literate

const SCRIPTS_DIR = joinpath(@__DIR__, "..", "scripts")
const TUTORIAL_MD_DIR = joinpath(@__DIR__, "src", "tutorials", "generated")
const NOTEBOOK_DIR = joinpath(@__DIR__, "generated_notebooks")
const CLEAN_SCRIPT_DIR = joinpath(@__DIR__, "generated_scripts")

const PUBLISHED_SCRIPTS = [
    "30_average_nlayers.jl" => "layered_sphere",
    "32_spheroid_effective_conductivity.jl" => "layered_spheroid_effective",
    "35_spheroid_interfaces.jl" => "layered_spheroid_interfaces",
    "37_spheroid_hc_conductivity.jl" => "layered_spheroid_hc",
    "70_symmetrization_showcase.jl" => "symmetrization",
]

function build_tutorial_pages()
    mkpath(TUTORIAL_MD_DIR)
    mkpath(NOTEBOOK_DIR)
    mkpath(CLEAN_SCRIPT_DIR)
    for (script, page) in PUBLISHED_SCRIPTS
        src = joinpath(SCRIPTS_DIR, script)
        Literate.markdown(src, TUTORIAL_MD_DIR; documenter = true, name = page)
        Literate.notebook(src, NOTEBOOK_DIR; name = page)
        Literate.script(src, CLEAN_SCRIPT_DIR; name = page)
    end
    return nothing
end

build_tutorial_pages()
