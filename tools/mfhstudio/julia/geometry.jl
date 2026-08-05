# =============================================================================
#  geometry.jl — 3-D traces for the web view.
#
#  The parametrizations come from `scripts/common/docviz.jl`, the same code
#  that draws the figures in the documentation, so a shape looks the same in
#  the interface as it does in the manual.
#
#  What is *not* reused is docviz's trace builders: they emit JavaScript object
#  literals (unquoted keys), which are not JSON and would have to be injected
#  into the page as executable script. The numeric layer -- `ellipsoid_surface`,
#  `cylinder_surface`, `disc_surface`, `param_surface`, `_layer_colors` --
#  returns plain arrays, so real JSON is built here instead.
# =============================================================================

using MeanFieldHom
using LinearAlgebra

const DOCVIZ = joinpath(pkgdir(MeanFieldHom), "scripts", "common", "docviz.jl")
isfile(DOCVIZ) || error("mfhstudio: cannot find $DOCVIZ")
include(DOCVIZ)

const FACE = "#4a90d9"
const CRACK = "#8a6d3b"
const GUIDE = "#c0392b"

_grid(M) = [collect(Float64, @view M[i, :]) for i in axes(M, 1)]

"""
    surface(X, Y, Z; color, opacity, name) -> Dict

A Plotly `surface` trace as a JSON-ready dictionary.
"""
function surface(X, Y, Z; color = FACE, opacity = 0.55, name = "")
    d = Dict{String, Any}(
        "type" => "surface",
        "x" => _grid(X), "y" => _grid(Y), "z" => _grid(Z),
        "opacity" => opacity,
        "showscale" => false,
        "colorscale" => [[0, color], [1, color]],
        "contours" => Dict(
            "x" => Dict("highlight" => false),
            "y" => Dict("highlight" => false),
            "z" => Dict("highlight" => false),
        ),
        "hoverinfo" => "skip",
    )
    isempty(name) || (d["name"] = name; d["showlegend"] = true)
    return d
end

function polyline(pts; color = GUIDE, width = 4, dash = "solid", name = "")
    P = pts isa AbstractMatrix ? pts : reduce(hcat, pts)
    d = Dict{String, Any}(
        "type" => "scatter3d", "mode" => "lines",
        "x" => collect(Float64, @view P[1, :]),
        "y" => collect(Float64, @view P[2, :]),
        "z" => collect(Float64, @view P[3, :]),
        "line" => Dict("color" => color, "width" => width, "dash" => dash),
        "hoverinfo" => "skip",
        # Guides are scenery, not data: naming them fills the legend with
        # "trace 1", "trace 2", … and hides the layers that do matter.
        "showlegend" => false,
    )
    isempty(name) || (d["name"] = name; d["showlegend"] = true)
    d
end

function labels3d(pts, texts; color = GUIDE)
    P = pts isa AbstractMatrix ? pts : reduce(hcat, pts)
    return Dict{String, Any}(
        "type" => "scatter3d", "mode" => "text",
        "x" => collect(Float64, @view P[1, :]),
        "y" => collect(Float64, @view P[2, :]),
        "z" => collect(Float64, @view P[3, :]),
        "text" => collect(String, texts),
        "textfont" => Dict("color" => color, "size" => 13),
        "hoverinfo" => "skip", "showlegend" => false,
    )
end

# ── Axis guides ─────────────────────────────────────────────────────────────

function axis_guides(R, lengths, names)
    out = Dict{String, Any}[]
    tips = Vector{Float64}[]
    for (i, (L, nm)) in enumerate(zip(lengths, names))
        isfinite(L) || (L = 3.0)
        e = R[:, i] .* L
        push!(out, polyline(hcat(zeros(3), e); color = GUIDE, width = 3))
        push!(tips, e .* 1.12)
    end
    push!(out, labels3d(reduce(hcat, tips), collect(String, names)))
    return out
end

_rot(b) = Matrix{Float64}(I, 3, 3)  # canonical fallback
function _rot(ell::MeanFieldHom.Ellipsoid)
    try
        return Matrix{Float64}(MeanFieldHom.TensND.get_array(ell.basis))
    catch
        return Matrix{Float64}(I, 3, 3)
    end
end

# ── Per-geometry traces ─────────────────────────────────────────────────────

traces(x; kw...) = Dict{String, Any}[]

function traces(ell::MeanFieldHom.Ellipsoid{3}; guides::Bool = true, kw...)
    a, b, c = Float64.(ell.semi_axes)
    R = _rot(ell)
    X, Y, Z = ellipsoid_surface(a, b, c; R = R)
    out = Dict{String, Any}[surface(X, Y, Z)]
    guides && append!(out, axis_guides(R, (a, b, c), ("a", "b", "c")))
    return out
end

function traces(cyl::MeanFieldHom.Cylinder; guides::Bool = true, length_shown = 6.0, kw...)
    b, c = Float64.(cyl.semi_axes)
    R = try
        Matrix{Float64}(MeanFieldHom.TensND.get_array(cyl.basis))
    catch
        Matrix{Float64}(I, 3, 3)
    end
    X, Y, Z = cylinder_surface(b, c, length_shown; R = R)
    out = Dict{String, Any}[surface(X, Y, Z)]
    guides && append!(out, axis_guides(R, (length_shown / 2, b, c), ("L → ∞", "b", "c")))
    return out
end

function _crack_traces(a, b, R; guides::Bool = true, normal::Bool = true, kw...)
    X, Y, Z = disc_surface(a, b; R = R)
    out = Dict{String, Any}[surface(X, Y, Z; color = CRACK, opacity = 0.9)]
    if normal
        n = R[:, 3] .* max(a, b) .* 1.5
        push!(out, polyline(hcat(zeros(3), n); color = GUIDE, width = 5))
        push!(out, labels3d(reshape(n .* 1.12, 3, 1), ["n"]))
    end
    guides && append!(out, axis_guides(R, (a, b, max(a, b) / 2), ("a", "b", "")))
    return out
end

function traces(cr::MeanFieldHom.EllipticCrack; kw...)
    R = try
        Matrix{Float64}(MeanFieldHom.TensND.get_array(cr.basis))
    catch
        Matrix{Float64}(I, 3, 3)
    end
    return _crack_traces(Float64(cr.a), Float64(cr.b), R; kw...)
end

function traces(cr::MeanFieldHom.RibbonCrack; length_shown = 6.0, kw...)
    R = try
        Matrix{Float64}(MeanFieldHom.TensND.get_array(cr.basis))
    catch
        Matrix{Float64}(I, 3, 3)
    end
    return _crack_traces(length_shown / 2, Float64(cr.b), R; kw...)
end

"""
    traces(ls::LayeredSphere; cutaway = true)

Concentric shells. The cut-away view is what makes a layered inclusion
readable at all -- without it only the outermost shell is visible.
"""
function traces(ls::MeanFieldHom.LayeredSphere; cutaway::Bool = true, guides::Bool = false, kw...)
    radii = Float64.(collect(ls.radii))
    n = length(radii)
    colors = _layer_colors(n)
    out = Dict{String, Any}[]
    # A cut-away keeps half of each shell, so the inner ones stay visible.
    ulim = cutaway ? 1.0π : 2.0π
    for (i, r) in enumerate(radii)
        r > 0 || continue
        f = (u, v) -> (r * cos(v) * cos(u), r * cos(v) * sin(u), r * sin(v))
        X, Y, Z = param_surface(f, range(0, ulim; length = 41), range(-π / 2, π / 2; length = 21))
        push!(
            out,
            surface(
                X, Y, Z;
                color = colors[i], opacity = i == n ? 0.55 : 0.9,
                name = "layer $i (r = $(round(r; digits = 3)))",
            ),
        )
    end
    return out
end

function traces(ls::MeanFieldHom.LayeredSpheroid; cutaway::Bool = true, guides::Bool = false, kw...)
    out = Dict{String, Any}[]
    axis = Float64.(collect(MeanFieldHom.outer_semiaxes(ls)))
    n = MeanFieldHom.layer_count(ls)
    colors = _layer_colors(n)
    ulim = cutaway ? 1.0π : 2.0π
    for i in 1:n
        c, aeq = try
            Float64.(MeanFieldHom.layer_semiaxes(ls, i))
        catch
            (axis[1] * i / n, axis[end] * i / n)
        end
        f = (u, v) -> (aeq * cos(v) * cos(u), aeq * cos(v) * sin(u), c * sin(v))
        X, Y, Z = param_surface(f, range(0, ulim; length = 41), range(-π / 2, π / 2; length = 21))
        push!(out, surface(X, Y, Z; color = colors[i], opacity = i == n ? 0.55 : 0.9, name = "layer $i"))
    end
    return out
end

"""
    scene(geom; kw...) -> Dict

A complete Plotly payload: `{data, layout}`, ready for `Plotly.newPlot`.
"""
function scene(geom; kw...)
    data = traces(geom; kw...)
    return Dict(
        "data" => data,
        "layout" => Dict(
            "margin" => Dict("l" => 0, "r" => 0, "t" => 0, "b" => 0),
            "showlegend" => length(data) > 1,
            "legend" => Dict("x" => 0, "y" => 1),
            "scene" => Dict(
                "aspectmode" => "data",
                "xaxis" => Dict("title" => "x"),
                "yaxis" => Dict("title" => "y"),
                "zaxis" => Dict("title" => "z"),
            ),
            "paper_bgcolor" => "rgba(0,0,0,0)",
        ),
    )
end
