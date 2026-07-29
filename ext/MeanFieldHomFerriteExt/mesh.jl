# =============================================================================
#  mesh.jl — the elliptical-crack-in-a-ball mesh.
#
#  Ported from the gmsh model of the `SifAniso` study (FEniCSx implementation
#  of the same corrected finite Eshelby cell), including the two traps that
#  cost real debugging time there and are re-checked here:
#
#   * the crack disc is `embed`ed in the ball, never `fragment`ed — fragment
#     would cut the ball into two half-balls, i.e. two disjoint volumes;
#   * the `Crack` plugin duplicates the nodes of the crack *front* as well as
#     those of the lips, in spite of `OpenBoundaryPhysicalGroup` (still true in
#     gmsh 4.15). Left alone, the crack is effectively half an element larger
#     than asked for and the opening is overestimated by 10-20 %. The front is
#     therefore welded back together after import — see `weld_crack_front`.
# =============================================================================

const TAG_MATRIX = 100
const TAG_SPHERE_EXT = 200
const TAG_CRACK = 300
const TAG_CRACK_FRONT = 400

const SET_OUTER = "outer"
const SET_CRACK = "crack"

"""
    _build_gmsh_crack_model(a, b, R, htipdiv)

Populate the current gmsh session with the crack-in-a-ball model: a ball of
radius `R` centred on an elliptical crack of semi-axes `(a, b)` lying in the
`z = 0` plane, refined to `min(a,b)/htipdiv` in a torus hugging the crack front
and coarsening to `R/3` at the outer boundary.

The caller owns `gmsh.initialize()` / `gmsh.finalize()`.
"""
function _build_gmsh_crack_model(a::Float64, b::Float64, R::Float64, htipdiv::Float64)
    gmsh.model.add("mfh_crack_in_ball")

    bmin = min(a, b)
    h_tip = bmin / htipdiv
    h_crack = bmin / 2
    h_far = R / 3
    r_flat = 0.4bmin          # radius of the fully-refined disc around the front
    r_tip = 0.75bmin           # end of the tip → crack-zone transition
    r_dense = min(0.4R, 3.0bmin)

    sphere = gmsh.model.occ.addSphere(0, 0, 0, R)
    disk = gmsh.model.occ.addDisk(0, 0, 0, a, b)
    gmsh.model.occ.synchronize()

    front = [abs(t) for (_, t) in gmsh.model.getBoundary([(2, disk)], false, false, false)]
    outer = [abs(t) for (_, t) in gmsh.model.getBoundary([(3, sphere)], false, false, false)]

    # `embed` constrains the mesher to place nodes on the disc while leaving the
    # ball a single volume; `occ.fragment` would split it in two.
    gmsh.model.mesh.embed(2, [disk], 3, sphere)

    gmsh.model.addPhysicalGroup(3, [sphere], TAG_MATRIX, "matrix")
    gmsh.model.addPhysicalGroup(2, outer, TAG_SPHERE_EXT, SET_OUTER)
    gmsh.model.addPhysicalGroup(2, [disk], TAG_CRACK, SET_CRACK)
    gmsh.model.addPhysicalGroup(1, front, TAG_CRACK_FRONT, "front")

    # Size field 1: distance to the crack front.  `SizeMax` must be `h_far`,
    # not `h_crack`, otherwise the `Min` combination clamps the whole domain.
    f1 = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(f1, "CurvesList", Float64.(front))
    gmsh.model.mesh.field.setNumber(f1, "Sampling", 200)
    f2 = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(f2, "InField", f1)
    gmsh.model.mesh.field.setNumber(f2, "SizeMin", h_tip)
    gmsh.model.mesh.field.setNumber(f2, "SizeMax", h_far)
    gmsh.model.mesh.field.setNumber(f2, "DistMin", r_flat)
    gmsh.model.mesh.field.setNumber(f2, "DistMax", r_tip)

    # Size field 2: radial coarsening.  A `MathEval` distance rather than a
    # `Distance`-to-a-point field, so that no stray geometric point is added —
    # it would survive as a node with no cell attached, i.e. a zero row in the
    # stiffness matrix.
    f3 = gmsh.model.mesh.field.add("MathEval")
    gmsh.model.mesh.field.setString(f3, "F", "Sqrt(x*x + y*y + z*z)")
    f4 = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(f4, "InField", f3)
    gmsh.model.mesh.field.setNumber(f4, "SizeMin", h_crack)
    gmsh.model.mesh.field.setNumber(f4, "SizeMax", h_far)
    gmsh.model.mesh.field.setNumber(f4, "DistMin", r_dense)
    gmsh.model.mesh.field.setNumber(f4, "DistMax", R)

    fmin = gmsh.model.mesh.field.add("Min")
    gmsh.model.mesh.field.setNumbers(fmin, "FieldsList", Float64[f2, f4])
    gmsh.model.mesh.field.setAsBackgroundMesh(fmin)

    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.Algorithm3D", 10)      # HXT

    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.optimize("Gmsh")

    # Split the lips: duplicate the nodes of the crack surface.
    gmsh.plugin.setNumber("Crack", "Dimension", 2)
    gmsh.plugin.setNumber("Crack", "PhysicalGroup", TAG_CRACK)
    gmsh.plugin.setNumber("Crack", "OpenBoundaryPhysicalGroup", TAG_CRACK_FRONT)
    gmsh.plugin.setNumber("Crack", "NormalX", 0.0)
    gmsh.plugin.setNumber("Crack", "NormalY", 0.0)
    gmsh.plugin.setNumber("Crack", "NormalZ", 1.0)
    gmsh.plugin.run("Crack")
    return nothing
end

"""
    weld_crack_front(grid, a, b) -> (grid, nweld)

Merge the coincident node pairs the `Crack` plugin leaves **on the crack
front**, so that the two lips close there and the discontinuity is bounded by
the ellipse the user asked for.

Nodes of the interior of the lips must stay doubled — that *is* the crack — so
the pairs to weld are identified geometrically: coincident, in the plane
`z = 0`, and on the ellipse `(x/a)² + (y/b)² = 1`. Duplicates produced by the
plugin are bit-identical copies, so the coincidence test is exact.

Cell ordering is untouched, so every facet set survives unchanged.
"""
function weld_crack_front(grid::Ferrite.Grid, a::Real, b::Real; ell_tol = 1.0e-6, z_tol = 1.0e-9)
    nodes = Ferrite.getnodes(grid)
    nn = length(nodes)

    rep = Dict{NTuple{3, Float64}, Int}()
    remap = collect(1:nn)
    nweld = 0
    for i in 1:nn
        x = nodes[i].x
        (abs(x[3]) < z_tol && abs((x[1] / a)^2 + (x[2] / b)^2 - 1) < ell_tol) || continue
        key = (x[1], x[2], x[3])
        j = get(rep, key, 0)
        if j == 0
            rep[key] = i
        else
            remap[i] = j
            nweld += 1
        end
    end
    nweld == 0 && return grid, 0

    keep = [i for i in 1:nn if remap[i] == i]
    newid = zeros(Int, nn)
    for (k, i) in enumerate(keep)
        newid[i] = k
    end
    for i in 1:nn
        newid[i] = newid[remap[i]]
    end

    CT = eltype(grid.cells)
    welded = Ferrite.Grid(
        [CT(map(n -> newid[n], c.nodes)) for c in grid.cells],
        [nodes[i] for i in keep];
        cellsets = Ferrite.getcellsets(grid),
        facetsets = Ferrite.getfacetsets(grid),
    )
    return welded, nweld
end

"""
    split_crack_lips(grid) -> (upper, lower)

Split the single imported `"crack"` facet set — which carries **both** lips —
into the two of them, by the sign of `z` of the centroid of the adjacent cell.
"""
function split_crack_lips(grid::Ferrite.Grid)
    up, dn = Set{Ferrite.FacetIndex}(), Set{Ferrite.FacetIndex}()
    for fi in Ferrite.getfacetset(grid, SET_CRACK)
        coords = Ferrite.getcoordinates(grid, fi[1])
        zc = sum(c[3] for c in coords) / length(coords)
        push!(zc >= 0 ? up : dn, fi)
    end
    return up, dn
end

"""
    facetset_area(grid, set, ip_geo, qr) -> Float64

Area of a facet set, by quadrature. Used to check the mesh against the exact
crack area `πab`.
"""
function facetset_area(grid::Ferrite.Grid, set, ip_geo, qr)
    fv = Ferrite.FacetValues(qr, ip_geo, ip_geo)
    A = 0.0
    for fi in set
        Ferrite.reinit!(fv, Ferrite.getcoordinates(grid, fi[1]), fi[2])
        for q in 1:Ferrite.getnquadpoints(fv)
            A += Ferrite.getdetJdV(fv, q)
        end
    end
    return A
end
