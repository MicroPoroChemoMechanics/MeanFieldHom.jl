# # n-layer confocal spheroid: local temperature and flux fields
#
# Pointwise reconstruction of the temperature and heat-flux fields
# in/around a 2-layer confocal prolate spheroid under a remote axial
# gradient. The reference local-field maps were cross-checked against a
# finite-element solution. The analogous script for the composite SPHERE
# is `31_local_nlayers.jl`.
#
# The API exercised here:
# - [`local_temperature`](@ref MeanFieldHom.LayeredSpheroids.local_temperature)`(s,
#   k₀, q, p, φ; H_axial, H_trans)`;
# - [`local_gradient`](@ref MeanFieldHom.LayeredSpheroids.local_gradient),
#   [`local_flux`](@ref MeanFieldHom.LayeredSpheroids.local_flux) — same
#   signature;
# - [`get_layer`](@ref MeanFieldHom.LayeredSpheroids.get_layer)`(s, q)` to
#   locate a query point.
#
# `(q, p, φ)` are the spheroid's OWN confocal coordinates (revolution
# axis ≡ `ê₃`); a meridian-plane point `(x, z)` (`φ = 0`, `x ≥ 0`) is
# recovered from the two-foci distances
# `q = (r₁+r₂)/(2c)`, `p = (r₁-r₂)/(2c)` (`r₁, r₂` distances to
# `(0,0,±c)`), the standard prolate-spheroidal-coordinate definition.

import Pkg                                                          #jl
Pkg.activate(joinpath(@__DIR__, ".."); io = devnull)                 #jl

using MeanFieldHom
using MeanFieldHom.LayeredSpheroids: local_temperature, local_gradient, local_flux
using TensND
using Printf
using Plots
gr()

# ## Setup — 2-layer confocal prolate spheroid, Kapitza interface at the core

const K1 = TensISO{3}(2.0)    # core
const K2 = TensISO{3}(15.0)   # shell
const K0 = TensISO{3}(1.0)    # matrix
const ρ = 0.15                # Kapitza resistance at the core/shell interface

const a_out, b_out = 3.0, 1.0
const focal2 = a_out^2 - b_out^2
const a_in = 2.9   # must exceed the shared focal distance √(a_out²-b_out²) for a real confocal shell
const b_in = sqrt(a_in^2 - focal2)
const c = sqrt(focal2)

s = LayeredSpheroid(
    (a_in, a_out), (b_in, b_out), (K1, K2);
    interfaces = (MeanFieldHom.KapitzaInterface(ρ), MeanFieldHom.PerfectInterface()),
    Nseries = 8,
)

# `(x, z)` (meridian half-plane, `φ=0`) → `(q, p)` via the two-foci
# distances (`r₁` to the `+c` focus, `r₂` to the `-c` focus): `q =
# (r₁+r₂)/(2c)`, `p = (r₂-r₁)/(2c)` — the sign on `p` is fixed by
# `z = c·q·p` (eq:xLeg), checked on-axis: `z > c ⟹ r₁ = z-c, r₂ = z+c
# ⟹ p = 1`, matching `z = c·q·1 = z`.
function _xz_to_qp(x, z)
    r1 = sqrt(x^2 + (z - c)^2)
    r2 = sqrt(x^2 + (z + c)^2)
    q = (r1 + r2) / (2c)
    p = (r2 - r1) / (2c)
    return q, clamp(p, -1.0, 1.0)
end

# ## Field grid

const nx, nz = 60, 90
const xs = range(1.0e-3, 4.5; length = nx)
const zs = range(-4.5, 4.5; length = nz)

T = Matrix{Float64}(undef, nz, nx)
for (j, z) in enumerate(zs), (i, x) in enumerate(xs)
    q, p = _xz_to_qp(x, z)
    T[j, i] = local_temperature(s, K0, q, p, 0.0; H_axial = 1.0, H_trans = 0.0)
end

# ## Flux streamlines
#
# The flux is shown as **integrated streamlines**, not as arrows on a grid. A
# quiver plot of this field is unreadable: normalizing every arrow to a common
# length throws away the magnitude — which is the whole story near an imperfect
# interface — while an arrow long enough to be visible is longer than the grid
# spacing, so the arrows overlap and hide the geometry they are meant to reveal.
#
# Streamlines avoid both problems: they are traced by RK4 on the **normalized**
# field ``\underline{q}/\|\underline{q}\|``, so line *spacing* carries the
# geometry, and the line *colour* carries ``\|\underline{q}\|``.

function _flux_xz(x, z)
    q, p = _xz_to_qp(max(abs(x), 1.0e-9), z)
    u1, _, u3 = local_flux(s, K0, q, p, 0.0; H_axial = 1.0, H_trans = 0.0)
    return u1 * (x < 0 ? -1.0 : 1.0), u3
end

_unit(u, v) = (n = hypot(u, v); n < 1.0e-12 ? (0.0, 0.0) : (u / n, v / n))

function _streamline(x0, z0; ds = 0.06, nmax = 1400, xmax = 4.6, zmax = 4.5)
    xs = Float64[x0]; zs = Float64[z0]; sp = Float64[hypot(_flux_xz(x0, z0)...)]
    x, z = x0, z0
    for _ in 1:nmax
        k1 = _unit(_flux_xz(x, z)...)
        k1 == (0.0, 0.0) && break
        k2 = _unit(_flux_xz(x + 0.5ds * k1[1], z + 0.5ds * k1[2])...)
        k3 = _unit(_flux_xz(x + 0.5ds * k2[1], z + 0.5ds * k2[2])...)
        k4 = _unit(_flux_xz(x + ds * k3[1], z + ds * k3[2])...)
        x += ds * (k1[1] + 2k2[1] + 2k3[1] + k4[1]) / 6
        z += ds * (k1[2] + 2k2[2] + 2k3[2] + k4[2]) / 6
        (abs(x) ≤ xmax && abs(z) ≤ zmax) || break
        push!(xs, x); push!(zs, z); push!(sp, hypot(_flux_xz(x, z)...))
    end
    return xs, zs, sp
end

const SEEDS = vcat(-4.0:0.4:-0.2, 0.2:0.4:4.0)

println("Local fields — 2-layer confocal prolate spheroid (a_out=$a_out, b_out=$b_out)")
println("─"^70)
@printf "core (a=%.2f,b=%.2f), Kapitza ρ=%.2f  |  shell (a=%.2f,b=%.2f), perfect\n" a_in b_in ρ a_out b_out
@printf "T at center (q≈1,p=1): %.4f   T far away (x=0,z=10): %.4f\n" local_temperature(
    s, K0, 1.0 + 1.0e-6, 1.0, 0.0; H_axial = 1.0,
) local_temperature(s, K0, _xz_to_qp(1.0e-6, 10.0)..., 0.0; H_axial = 1.0)
println()

# ## Graphical output — temperature map and flux streamlines
#
# Both panels show the full meridian section: the field is axisymmetric, so the
# `x < 0` half is the mirror image of the computed `x > 0` half.

const xs_full = vcat(-reverse(xs), xs)
const T_full = hcat(reverse(T; dims = 2), T)

p1 = heatmap(
    xs_full, zs, T_full; xlabel = "x", ylabel = "z",
    title = "Temperature field (axial remote gradient)",
    color = :thermal, aspect_ratio = :equal, colorbar_title = "T",
)
# Confocal interface boundaries, meridian trace: x²/b² + z²/a² = 1.
for (a, b) in ((a_in, b_in), (a_out, b_out))
    tt = range(0, 2π; length = 300)
    plot!(p1, b .* sin.(tt), a .* cos.(tt); color = :white, lw = 1.5, label = false)
end

p2 = plot(;
    xlabel = "x", ylabel = "z", title = "Heat flux streamlines (colour = ‖q‖)",
    aspect_ratio = :equal, legend = false, xlims = (-4.6, 4.6), ylims = (-4.5, 4.5),
    colorbar = true, colorbar_title = "‖q‖",
)
## Seeds sit on the INFLOW edge: with `H_axial > 0` the gradient points along
## +z, so the flux `q = -k∇T` runs downward and the lines enter from the top.
for x0 in SEEDS
    sx, sz, sp = _streamline(x0, 4.4)
    length(sx) > 3 || continue
    plot!(p2, sx, sz; line_z = sp, color = :thermal, lw = 1.1)
end
for (a, b) in ((a_in, b_in), (a_out, b_out))
    tt = range(0, 2π; length = 300)
    plot!(p2, b .* sin.(tt), a .* cos.(tt); color = :black, lw = 1.5, label = false)
end

p_full = plot(p1, p2; layout = (1, 2), size = (1200, 600))
p_full

const figdir = joinpath(@__DIR__, "figures")                        #jl
isdir(figdir) || mkdir(figdir)                                       #jl
figpath = joinpath(figdir, "35_spheroid_local_fields.png")            #jl
savefig(p_full, figpath)                                              #jl
display(p_full)                                                       #jl
@printf "\nSaved : %s\n" figpath                                      #jl
