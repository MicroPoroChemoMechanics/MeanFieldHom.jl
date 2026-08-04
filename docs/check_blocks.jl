#!/usr/bin/env julia
# =============================================================================
#  check_blocks.jl — a partial documentation check.
#
#  A full `docs/make.jl` takes tens of minutes, most of it spent re-executing
#  pages nobody just edited. This script runs only the `@setup` / `@example` /
#  `@repl` blocks of the pages you name, the way Documenter does: one fresh
#  module per block group, blocks of the same group sharing state in page order.
#
#  Usage
#
#      julia --project=docs docs/check_blocks.jl                       # every page
#      julia --project=docs docs/check_blocks.jl theory/hill_tensors.md
#      julia --project=docs docs/check_blocks.jl manual/ theory/notation.md
#      julia --project=docs docs/check_blocks.jl $(git diff --name-only -- 'docs/src/*.md')
#
#  Arguments are paths to markdown pages or to directories, resolved relative to
#  the repository root, to `docs/src/`, or absolute — whichever matches first.
#
#  What it catches: undefined names, shadowed bindings (`strip = ...` over
#  `Base.strip`), state that a block silently relied on from another page,
#  method errors, and anything else that throws. Exit status is non-zero if any
#  block fails, so it can gate a commit.
#
#  What it does NOT catch — run the full build before publishing:
#    * cross-references (`@ref`), citations, the `pages` tree, the search index;
#    * page-size thresholds, since nothing is rendered to HTML here;
#    * whether a figure comes out *blank* — for the interactive Plotly scenes
#      that needs a real browser (see the note at the end of this file).
# =============================================================================

using Printf

const DOCS_SRC = joinpath(@__DIR__, "src")
const REPO = dirname(@__DIR__)

# GR needs a headless display driver, exactly as in `make.jl`.
ENV["GKSwstype"] = "100"

"""
    resolve(arg) -> Vector{String}

Turn one command-line argument into a list of markdown pages.
"""
function resolve(arg::AbstractString)
    for cand in (arg, joinpath(REPO, arg), joinpath(DOCS_SRC, arg))
        isfile(cand) && return [abspath(cand)]
        if isdir(cand)
            pages = String[]
            for (root, _, files) in walkdir(cand), f in files
                endswith(f, ".md") && push!(pages, abspath(joinpath(root, f)))
            end
            return sort(pages)
        end
    end
    error("no such page or directory: $arg")
end

"""
    block_groups(path) -> Vector{Pair{String, Vector{String}}}

Extract the executable blocks of a page, grouped by Documenter tag and kept in
page order. An untagged block gets a group of its own, which is how Documenter
sandboxes it.
"""
function block_groups(path)
    lines = readlines(path)
    groups = Pair{String, Vector{String}}[]
    index = Dict{String, Int}()
    i, anon = 1, 0
    while i ≤ length(lines)
        m = match(r"^```@(setup|example|repl)(\s+(\S+))?\s*$", lines[i])
        if m === nothing
            i += 1
            continue
        end
        tag = m.captures[3]
        if tag === nothing
            anon += 1
            tag = "__anonymous_$(anon)__"
        end
        j = i + 1
        body = String[]
        while j ≤ length(lines) && !startswith(lines[j], "```")
            push!(body, lines[j])
            j += 1
        end
        k = get(index, tag, 0)
        if k == 0
            push!(groups, tag => [join(body, "\n")])
            index[tag] = length(groups)
        else
            push!(groups[k].second, join(body, "\n"))
        end
        i = j + 1
    end
    return groups
end

function check(path)
    rel = relpath(path, DOCS_SRC)
    failures = 0
    for (tag, chunks) in block_groups(path)
        sandbox = Module(Symbol("Sandbox_", hash(rel * tag)))
        ## Documenter's sandboxes provide `include`; a bare `Module` does not,
        ## and several pages include `scripts/common/docviz.jl`.
        Core.eval(sandbox, :(include(p) = Base.include($sandbox, p)))
        for (n, code) in enumerate(chunks)
            try
                Base.include_string(sandbox, code, "$rel [$tag #$n]")
            catch err
                failures += 1
                msg = first(sprint(showerror, err), 300)
                @printf "  FAIL  %s  block [%s #%d]\n        %s\n" rel tag n msg
            end
        end
    end
    nblocks = sum(length(g.second) for g in block_groups(path); init = 0)
    if nblocks == 0
        @printf "  --    %s  (no executable block)\n" rel
    elseif failures == 0
        @printf "  ok    %s  (%d block%s)\n" rel nblocks (nblocks == 1 ? "" : "s")
    end
    return failures
end

function main(args)
    pages = isempty(args) ? resolve(DOCS_SRC) :
        unique(reduce(vcat, resolve.(args)))
    @printf "Checking %d page%s\n\n" length(pages) (length(pages) == 1 ? "" : "s")
    total = sum(check(p) for p in pages; init = 0)
    println()
    if total == 0
        println("All blocks ran.")
    else
        @printf "%d failing block%s.\n" total (total == 1 ? "" : "s")
    end
    return total == 0 ? 0 : 1
end

exit(main(ARGS))

# ── Checking that an interactive figure is not blank ─────────────────────────
#
# Documenter reports nothing when a Plotly scene renders empty, so the only real
# check is a browser. After a build:
#
#     python3 -m http.server 8000 --directory docs/build &
#     chrome --headless=new --no-sandbox --use-gl=angle --use-angle=swiftshader \
#            --enable-unsafe-swiftshader --window-size=1280,1600 \
#            --virtual-time-budget=15000 --screenshot=/tmp/page.png \
#            http://localhost:8000/manual/inclusion_gallery.html
#
# The software-WebGL flags matter: without them plotly.js reports "WebGL is not
# supported by your browser" and every 3-D figure comes out as a gray box, which
# looks like a bug in the figure and is not one. `file://` URLs do not work
# either — the require.js fetch of plotly.js from the CDN needs a real origin.
