# docs/literate.jl
#
# Runs Literate.jl over the curated list of `scripts/` demos that have been
# converted to the "dual-usage" contract (see `scripts/README.md`, §Literate
# convention) and are promoted into the doc site's Gallery section. Produces,
# per script, three artifacts:
#
#   - a Documenter-ready markdown page  -> docs/src/gallery/generated/
#     (executed by Documenter itself via `@example` when `makedocs` runs —
#     Literate's `markdown(...; documenter=true)` leaves `execute=false`)
#   - a pre-run Jupyter notebook        -> docs/generated_notebooks/
#   - a cleaned standalone .jl script   -> docs/generated_scripts/
#     (markup-only lines stripped, `#src`/`#nb`/`#md` directives resolved)
#
# Called from `docs/make.jl` *before* `makedocs`, so the generated markdown
# exists when Documenter's `pages` list references it.
#
# GALLERY_SCRIPTS is intentionally curated, not "every script in scripts/":
# only scripts classified "G" (gap-filler, no tutorial/application overlap)
# in `Assets/plans/MFH_LITERATE_SCRIPTS.md` belong here. Scripts that
# duplicate an existing tutorial ("D") stay plain scripts — converting them
# would put two pages saying the same thing in the site. See that file for
# the full classification and rationale.

using Literate

const SCRIPTS_DIR = joinpath(@__DIR__, "..", "scripts")
const GALLERY_MD_DIR = joinpath(@__DIR__, "src", "gallery", "generated")
const GALLERY_NB_DIR = joinpath(@__DIR__, "generated_notebooks")
const GALLERY_SCRIPT_DIR = joinpath(@__DIR__, "generated_scripts")

# Phase-1 pilot: only the two "G" (gap-filler) pilots are promoted to the
# gallery. `29_symbolic_schemes.jl` was also converted to the dual-usage
# contract as a pilot (to measure SymPy-under-`@example` build cost) but is
# a "D" (duplicate of tutorial 11) — deliberately NOT listed here, so it
# does not appear as a competing page in the built site.
const GALLERY_SCRIPTS = [
    "70_symmetrization_showcase.jl",
    "30_average_nlayers.jl",
]

function build_gallery()
    mkpath(GALLERY_MD_DIR)
    mkpath(GALLERY_NB_DIR)
    mkpath(GALLERY_SCRIPT_DIR)
    for name in GALLERY_SCRIPTS
        src = joinpath(SCRIPTS_DIR, name)
        Literate.markdown(src, GALLERY_MD_DIR; documenter = true)
        Literate.notebook(src, GALLERY_NB_DIR)
        Literate.script(src, GALLERY_SCRIPT_DIR)
    end
end

build_gallery()
