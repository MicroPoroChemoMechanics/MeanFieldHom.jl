# # What an imperfect interface actually does to the flux
#
# A single picture explains imperfect interfaces better than any formula: put an
# **insulating** particle in a conducting matrix, wrap it in a **highly
# conducting** surface layer, and watch the heat flux change its mind about
# where to go.
#
# With no surface layer the particle is an obstacle — the flux parts around it.
# Give the interface a surface conductance ``\beta`` and the particle's *skin*
# becomes a preferential path: flux is carried around the particle inside a
# layer of zero thickness. Past a certain ``\beta`` the **insulating** particle
# stops behaving like a hole and starts behaving like a conductor — its
# equivalent conductivity crosses that of the bare matrix.
#
# The two interface models available on a [`LayeredSpheroid`](@ref) are
# ```math
# \text{LC (Kapitza):}\quad [\![T]\!] = \rho\,q_n,
# \qquad
# \text{HC (surface-conductive):}\quad [\![q_n]\!] = -\beta\,\mathrm{div}_S(\nabla_S T),
# ```
# with ``\rho`` a genuine thermal resistance and ``\beta`` a genuine surface
# conductance. This script uses the **HC** one
# ([`SurfaceConductiveInterface`](@ref)); see
# [`32_spheroid_nlayers_conductivity.jl`](32_spheroid_nlayers_conductivity.md)
# for the LC counterpart at the level of effective properties.
#
# The configuration reproduces the one used in the ECHOES presentation of
# 06/07/2020: an **oblate** spheroid of aspect ratio ``1/2``, insulating core,
# unit matrix conductivity, HC interface swept over ``\beta \in [0, 3]``, with
# ``\mathcal{N} = 5`` series terms. Theory:
# [Layered spheroid](../../theory/layered_spheroid.md);
# reference [barthelemyBignonnetIJES2020](@cite), interface models
# [kushch2015](@cite), [miloh1999](@cite).

import Pkg                                                          #jl
Pkg.activate(joinpath(@__DIR__, ".."); io = devnull)                 #jl

using MeanFieldHom
using MeanFieldHom.LayeredSpheroids: local_flux, local_temperature
using TensND
using Printf
using Plots
gr()

# ## Geometry: an oblate spheroid and its confocal coordinates
#
# Oblate means the revolution semi-axis is the *shorter* one:
# `axis_radius < disk_radius`. The two share a focal **ring** of radius
# ``\bar c = \sqrt{\rho_t^2 - \rho_a^2}``, and the confocal parameter is
# ``q = \mathrm{i}\,\tau`` with ``\tau = \rho_a/\bar c`` real — the complex
# substitution that lets every prolate formula carry over unchanged.

const AXIS = 1.0                      # ρ_a, revolution semi-axis (short)
const DISK = 2.0                      # ρ_t, transverse semi-axis (long)
const CBAR = sqrt(DISK^2 - AXIS^2)    # focal ring radius
const KM = 1.0                        # matrix conductivity
const KCORE = 1.0e-6                  # insulating core (numerically safe zero)
const NSERIES = 5

# Meridian point `(x, z)` → oblate coordinates `(τ, p)`, inverting
# `x = c̄√(1+τ²)√(1-p²)`, `z = c̄ τ p`. Solving the quadratic in `u = τ²`:
# ``c̄²u² + u(c̄² - x² - z²) - z² = 0``.

function _xz_to_taup(x, z)
    s = x^2 + z^2 - CBAR^2
    u = (s + sqrt(s^2 + 4 * CBAR^2 * z^2)) / (2 * CBAR^2)
    τ = sqrt(max(u, 0.0))
    p = τ > 1.0e-12 ? clamp(z / (CBAR * τ), -1.0, 1.0) : 0.0
    return τ, p
end

# Build the particle for a given surface conductance `β`. A single layer (`N=1`)
# carries a single interface, at the particle/matrix boundary.

function _particle(β)
    return LayeredSpheroid(
        (AXIS,), (DISK,), (TensISO{3}(KCORE),);
        interfaces = (MeanFieldHom.SurfaceConductiveInterface(β),),
        Nseries = NSERIES,
    )
end

# ## Streamlines
#
# The flux field is integrated directly rather than drawn as arrows: a quiver
# plot of a field this strongly channelled is unreadable, because the arrows
# that matter (inside the interface skin) are exactly the ones a uniform arrow
# length hides. Streamlines are traced by RK4 on the **normalized** field
# ``\underline{q}/\|\underline{q}\|``, so the line spacing carries the geometry
# and the colour carries the magnitude.

const K0 = TensISO{3}(KM)

function _flux_xz(s, x, z)
    τ, p = _xz_to_taup(max(abs(x), 1.0e-9), z)
    u1, _, u3 = local_flux(s, K0, im * τ, p, 0.0; H_axial = 1.0, H_trans = 0.0)
    ux = real(u1) * sign(x == 0 ? 1.0 : x)
    return ux, real(u3)
end

function _unit(u, v)
    n = hypot(u, v)
    return n < 1.0e-12 ? (0.0, 0.0) : (u / n, v / n)
end

## `qmin` stops a line where the bulk flux dies out — inside the insulating
## core. That matters here: an HC interface carries a *surface* current, so the
## bulk flux inside the particle is numerically tiny but not exactly zero, and
## integrating the normalized field through it would draw straight lines across
## the particle that represent nothing physical.
function _streamline(s, x0, z0; ds = 0.05, nmax = 1200, xmax = 5.0, zmax = 4.0, qmin = 0.02)
    xs = Float64[x0]; zs = Float64[z0]; sp = Float64[hypot(_flux_xz(s, x0, z0)...)]
    x, z = x0, z0
    for _ in 1:nmax
        hypot(_flux_xz(s, x, z)...) < qmin && break
        k1 = _unit(_flux_xz(s, x, z)...)
        k1 == (0.0, 0.0) && break
        k2 = _unit(_flux_xz(s, x + 0.5ds * k1[1], z + 0.5ds * k1[2])...)
        k3 = _unit(_flux_xz(s, x + 0.5ds * k2[1], z + 0.5ds * k2[2])...)
        k4 = _unit(_flux_xz(s, x + ds * k3[1], z + ds * k3[2])...)
        x += ds * (k1[1] + 2k2[1] + 2k3[1] + k4[1]) / 6
        z += ds * (k1[2] + 2k2[2] + 2k3[2] + k4[2]) / 6
        (abs(x) ≤ xmax && abs(z) ≤ zmax) || break
        push!(xs, x); push!(zs, z); push!(sp, hypot(_flux_xz(s, x, z)...))
    end
    return xs, zs, sp
end

# Seeds are spread across the inflow edge, denser near the axis where the
# particle actually perturbs the field.

const SEEDS = vcat(-3.6:0.45:-0.2, 0.2:0.45:3.6)

function _streamplot(β; title = "")
    s = _particle(β)
    plt = plot(;
        aspect_ratio = :equal, xlims = (-4.2, 4.2), ylims = (-3.4, 3.4),
        xlabel = "x", ylabel = "z", title = title, legend = false,
        colorbar = true, colorbar_title = "‖q‖",
    )
    ## Seeds sit on the INFLOW edge: with `H_axial > 0` the gradient points
    ## along +z, so the flux `q = -k∇T` runs downward, entering from the top.
    for x0 in SEEDS
        sx, sz, sp = _streamline(s, x0, 3.3)
        length(sx) > 3 || continue
        plot!(plt, sx, sz; line_z = sp, color = :thermal, lw = 1.1, clims = (0.0, 2.2))
    end
    ## Particle outline, meridian trace: x²/DISK² + z²/AXIS² = 1. The interface
    ## is drawn with a width proportional to β — the ECHOES convention, and a
    ## reminder that the "layer" has zero thickness but finite conductance.
    tt = range(0, 2π; length = 300)
    plot!(plt, DISK .* cos.(tt), AXIS .* sin.(tt); color = :black, lw = 1 + 2.5β)
    return plt
end

# ## Before and after
#
# At `β = 0` the interface is perfect and the particle is a pure obstacle: the
# flux parts around it and **piles up at the equator**, the bright spots where
# ``\|\underline{q}\|`` is largest.
#
# At `β = 3` that concentration is gone and the streamlines run almost
# undisturbed past the particle. The skin is short-circuiting the obstacle: it
# accepts the flux tangentially, carries it around, and returns it downstream,
# so the far field barely registers that anything is there.
#
# Note what the lines do *not* show. The streamlines stop at the interface
# because they trace the **bulk** flux, and inside an insulating core there is
# none. The current carried by an HC interface is a genuine *surface* current,
# living on a layer of zero thickness — it cannot appear as a streamline. That
# is why the interface is drawn with a width proportional to ``\beta``, and why
# the quantitative statement below is needed to complete the picture.

p_off = _streamplot(0.0; title = "β = 0  (perfect interface, insulating particle)")
p_on = _streamplot(3.0; title = "β = 3  (highly conducting interface)")
p_pair = plot(p_off, p_on; layout = (1, 2), size = (1250, 520))
p_pair

# ## How much flux does the skin actually carry?
#
# A scalar summary of the same effect: the axial and transverse components of
# the equivalent particle conductivity ``\boldsymbol{k}^{eq} =
# \boldsymbol{B}_\Omega\cdot\boldsymbol{A}_\Omega^{-1}``, where
# ``\boldsymbol{A}_\Omega`` and ``\boldsymbol{B}_\Omega`` are the volume-averaged
# gradient and flux concentration tensors of the particle
# (`⟨∇T⟩_Ω = A_Ω·H`, `⟨K∇T⟩_Ω = B_Ω·H`).

βs = range(0.0, 3.0; length = 25)
keq_a = Float64[]; keq_t = Float64[]
for β in βs
    s = _particle(β)
    A = gradient_gradient_loc(s, TensISO{3}(KCORE), K0)
    B = flux_gradient_loc(s, TensISO{3}(KCORE), K0)
    push!(keq_t, real(B[1, 1] / A[1, 1]))
    push!(keq_a, real(B[3, 3] / A[3, 3]))
end

println("HC interface on an insulating oblate spheroid (ρ_a=$AXIS, ρ_t=$DISK)")
println("─"^70)
@printf "β = 0.0 :  k_t^eq/k_m = %6.3f   k_a^eq/k_m = %6.3f\n" keq_t[1]/KM keq_a[1]/KM
@printf "β = 3.0 :  k_t^eq/k_m = %6.3f   k_a^eq/k_m = %6.3f\n" keq_t[end]/KM keq_a[end]/KM
println()

p_keq = plot(
    βs, keq_t ./ KM; label = "transverse  k_t^eq / k_m", lw = 2.5, color = :steelblue,
    xlabel = "surface conductance β", ylabel = "k^eq / k_m",
    title = "Equivalent particle conductivity vs. interface conductance",
    legend = :topleft, size = (760, 460),
)
plot!(p_keq, βs, keq_a ./ KM; label = "axial  k_a^eq / k_m", lw = 2.5, color = :darkorange)
hline!(p_keq, [1.0]; color = :grey, ls = :dash, label = "matrix (k_m)")
p_keq

# ## The sweep, animated
#
# The same streamline picture as ``\beta`` grows from 0 to 3 — the frames of the
# ECHOES animation, regenerated from the Julia solution.

anim = @animate for β in range(0.0, 3.0; length = 10)
    _streamplot(β; title = @sprintf("HC interface,  β = %.2f", β))
end
## No output path here on purpose: `gif(anim; fps)` writes to a temporary file
## that Documenter and the notebook embed directly. Passing a path built from
## `@__DIR__` would break both, because in the generated notebook `@__DIR__` is
## `docs/generated_notebooks/`, which has no `figures/` subdirectory — ffmpeg
## then fails with a bare exit 254. The standalone copy is written in the `#jl`
## block at the end of the script instead.
gif(anim; fps = 2)

# ## Interactive 3-D view
#
# The meridian streamlines revolved around the axis, as an interactive Plotly
# scene: drag to rotate, scroll to zoom. The grey surface is the particle; the
# coloured curves are flux lines, coloured by ``\|\underline{q}\|``.

function _plotly_streamlines(β; nrev = 8, uid = "spheroid-3d")
    s = _particle(β)
    traces = String[]
    for x0 in (-3.2, -2.2, -1.3, -0.6, 0.6, 1.3, 2.2, 3.2), krev in 0:(nrev - 1)
        sx, sz, sp = _streamline(s, x0, 3.0; ds = 0.12, nmax = 300)
        length(sx) > 4 || continue
        θ = 2π * krev / nrev
        X = [x * cos(θ) for x in sx]
        Y = [x * sin(θ) for x in sx]
        push!(
            traces, """{type:"scatter3d",mode:"lines",
            x:[$(join(round.(X; digits=3), ","))],
            y:[$(join(round.(Y; digits=3), ","))],
            z:[$(join(round.(sz; digits=3), ","))],
            line:{width:3,color:[$(join(round.(sp; digits=3), ","))],
                  colorscale:"Hot",cmin:0,cmax:2.2},showlegend:false}"""
        )
    end
    ## Particle surface (oblate spheroid), coarse parametric mesh.
    nu, nv = 24, 16
    us = range(0, 2π; length = nu); vs = range(-π / 2, π / 2; length = nv)
    Xs = [DISK * cos(v) * cos(u) for v in vs, u in us]
    Ys = [DISK * cos(v) * sin(u) for v in vs, u in us]
    Zs = [AXIS * sin(v) for v in vs, _ in us]
    jsmat(M) = "[" * join(("[" * join(round.(M[i, :]; digits = 3), ",") * "]" for i in axes(M, 1)), ",") * "]"
    push!(
        traces, """{type:"surface",x:$(jsmat(Xs)),y:$(jsmat(Ys)),z:$(jsmat(Zs)),
        opacity:0.35,showscale:false,colorscale:[[0,"#888"],[1,"#888"]]}"""
    )
    return Base.HTML(
        """
        <div id="$uid" style="width:100%;height:560px;"></div>
        <script>
        (function () {
          var data = [$(join(traces, ","))];
          var layout = {height:560, margin:{l:0,r:0,t:30,b:0},
            title:{text:"Flux lines around an insulating spheroid, β = $β"},
            scene:{aspectmode:"data", xaxis:{title:"x"}, yaxis:{title:"y"}, zaxis:{title:"z"}}};
          function draw(Plotly) { Plotly.newPlot("$uid", data, layout); }
          if (window.Plotly) { draw(window.Plotly); }
          else if (window.require) {
            require.config({paths: {plotly_mfh: "https://cdn.plot.ly/plotly-2.35.2.min"}});
            require(["plotly_mfh"], draw);
          }
        })();
        </script>
        """
    )
end

_plotly_streamlines(3.0; uid = "spheroid-hc-3d")

const figdir = joinpath(@__DIR__, "figures")                        #jl
isdir(figdir) || mkdir(figdir)                                       #jl
savefig(p_pair, joinpath(figdir, "36_interface_effect_pair.png"))     #jl
savefig(p_keq, joinpath(figdir, "36_interface_effect_keq.png"))       #jl
gif(anim, joinpath(figdir, "36_interface_effect.gif"); fps = 2, show_msg = false)  #jl
display(p_pair)                                                       #jl
@printf "\nSaved : %s\n" joinpath(figdir, "36_interface_effect_*.png") #jl
