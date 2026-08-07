"""Tests for MFH Studio.

    python3 tests/test_studio.py            model, codegen, graph — no Julia
    python3 tests/test_studio.py --julia    adds the sidecar-backed tests

The load-bearing test is `test_preserves_every_demo_script`: the interface may
only be pointed at somebody's existing work if opening and saving cannot damage
it.
"""

from __future__ import annotations

import glob
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))

from mfhstudio.codegen import extract_embedded, generate, render_cell  # noqa: E402
from mfhstudio.model import (  # noqa: E402
    Cell,
    Geometry,
    Lens,
    Model,
    Param,
    Phase,
    Property,
    Sweep,
    default_model,
)

REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
WITH_JULIA = "--julia" in sys.argv


# ---------------------------------------------------------------------------
# Conventions the interface exists to remove
# ---------------------------------------------------------------------------


def test_matrix_has_no_amount():
    """MFH derives it as 1 - Σ f and raises if it is set."""
    src = generate(default_model(), embed_model=False)
    add_matrix = next(l for l in src.splitlines() if "add_matrix!" in l)
    assert "fraction" not in add_matrix


def test_physical_moduli_not_raw_tensiso():
    """`iso_stiffness(k, μ)` takes physical moduli; TensISO{3} takes (3k, 2μ)."""
    src = generate(default_model(), embed_model=False)
    assert "iso_stiffness(72.0, 32.0)" in src
    assert "TensISO{3}(216" not in src


def test_solver_options_attach_to_the_scheme():
    m = default_model()
    m.sweep.schemes = [{
        "name": "SelfConsistent",
        "options": {"abstol": 1e-10, "maxiters": 300, "select_best": True},
    }]
    src = generate(m, embed_model=False)
    assert "SelfConsistent(; abstol = 1.0e-10, maxiters = 300, select_best = true)" in src
    assert "homogenize(cell, scheme, :C)" in src


def test_several_schemes_share_one_figure():
    m = default_model()
    m.sweep.schemes = [
        {"name": "MoriTanaka", "options": {}},
        {"name": "Voigt", "options": {}},
    ]
    src = generate(m, embed_model=False)
    assert '("MoriTanaka", MoriTanaka())' in src
    assert '("Voigt", Voigt())' in src
    assert "for (name, scheme) in SCHEMES" in src


def test_single_point_uses_the_amounts_as_entered():
    """The answer to "just compute with the fractions I typed"."""
    m = default_model()
    m.sweep.mode = "single"
    src = generate(m, embed_model=False)
    assert "One homogenization with the amounts entered" in src
    assert "set_param" not in src
    assert "homogenize(cell, MoriTanaka(), :C)" in src


def test_kelvin_mandel_output_needs_no_isotropy():
    """`k_mu` has a method for TensISO alone; an oriented inclusion with no
    orientation average does not give one, and the run dies with a
    MethodError deep inside. Components are defined whatever the symmetry."""
    m = default_model()
    m.sweep.projection = "none"
    m.sweep.outputs = [{"kind": "km", "i": 1, "j": 1}, {"kind": "km", "i": 3, "j": 3}]
    src = generate(m, embed_model=False)
    assert "KM(C)[1, 1]" in src and "KM(C)[3, 3]" in src
    assert "k_mu" not in src


def test_isotropic_only_output_without_a_projection_is_flagged():
    m = default_model()
    m.sweep.projection = "none"
    m.sweep.outputs = [{"kind": "k"}]
    for c in m.cells:
        for ph in c.phases:
            ph.symmetrize = "none"
    assert any("isotropic result" in p for p in m.validate())


def test_viscoelastic_laws_use_the_real_signatures():
    """`maxwell_iso` takes two relaxation times, not one."""
    from mfhstudio.codegen import CodeGen

    g = CodeGen(Model())
    assert g._prop_expr(Property(
        builder="maxwell_iso",
        args={"k": 10.0, "mu": 5.0, "eta_k": 2.0, "eta_mu": 3.0},
    )) == "maxwell_iso(10.0, 5.0, 2.0, 3.0)"
    assert g._prop_expr(Property(
        builder="kelvin_iso",
        args={"k0": 10.0, "mu0": 5.0, "k1": 20.0, "mu1": 10.0,
              "tau_k": 1.0, "tau_mu": 2.0},
    )) == "kelvin_iso(10.0, 5.0, [20.0], [10.0], [1.0], [2.0])"


def test_anisotropic_conductivity_forms():
    from mfhstudio.codegen import CodeGen

    g = CodeGen(Model())
    # The OUTER constructor: `TensTI{2, Float64, 2}(data, n)` demands a 3-tuple
    # axis and rejects the vector an oriented frame yields.
    assert g._prop_expr(Property(builder="TensTI2", args={"kt": 1.0, "ka": 5.0})) == (
        "TensTI{2}(1.0, 5.0, (0.0, 0.0, 1.0))"
    )
    assert g._prop_expr(Property(
        builder="TensDiag2", args={"k1": 1.0, "k2": 2.0, "k3": 5.0}
    )) == "Tens([1.0 0.0 0.0; 0.0 2.0 0.0; 0.0 0.0 5.0])"


def test_anisotropic_properties_carry_their_own_frame():
    """The frame a tensor's constants are written in is not the shape's.

    A transversely isotropic tensor takes an axis (the third vector of the
    frame); an orthotropic one takes the basis itself, its components then
    being read *in* that basis.
    """
    from mfhstudio.codegen import CodeGen

    g = CodeGen(Model())
    ang = ["pi/4", 0.7, 0.0]

    ti = g._prop_expr(Property(builder="TensTI2", args={"kt": 1.0, "ka": 5.0},
                               euler_angles=ang))
    assert ti == "TensTI{2}(1.0, 5.0, vecbasis(RotatedBasis(pi/4, 0.7, 0.0))[:, 3])"

    # `hoenig_stiffness` declares no five-argument method: the axis is required.
    hoenig = g._prop_expr(Property(builder="hoenig_stiffness", euler_angles=ang))
    assert hoenig.startswith("hoenig_stiffness(")
    assert hoenig.endswith("vecbasis(RotatedBasis(pi/4, 0.7, 0.0))[:, 3])")
    assert g._prop_expr(Property(builder="hoenig_stiffness")).endswith("(0.0, 0.0, 1.0))")

    ortho = g._prop_expr(Property(builder="TensDiag2",
                                  args={"k1": 1.0, "k2": 2.0, "k3": 5.0},
                                  euler_angles=ang))
    assert ortho.endswith(", RotatedBasis(pi/4, 0.7, 0.0))")

    # `RotatedBasis` with fewer than three angles builds a 2-D basis, so the
    # frame is always spelled out in full.
    assert g._prop_expr(Property(builder="TensOrtho", euler_angles=[0.3])) \
        .endswith("RotatedBasis(0.3, 0.0, 0.0))")


def test_hoenig_defaults_are_not_the_isotropic_point():
    """h = 1 with ν₁ = ν₂ and γ = 1 is isotropy wearing a TI type."""
    from mfhstudio.catalog import PROPERTIES

    form = next(f for f in PROPERTIES if f["name"] == "ti_hoenig")
    d = {f["name"]: f["default"] for f in form["fields"]}
    assert not (d["h"] == 1.0 and d["nu1"] == d["nu2"] and d["gamma"] == 1.0)


def test_alv_curve_follows_the_documented_extraction():
    m = default_model()
    m.alv.enabled = True
    src = generate(m, embed_model=False)
    assert "volterra_inverse(R; block_size = 6)" in src
    assert "homogenize_alv(" in src
    assert "using Plots" in src, "the ALV run plots too"


def _layered_spheroid_model() -> Model:
    def layer(fr, k):
        return {
            "fraction": fr,
            "property": {
                "key": ":K", "source": "builder", "builder": "TensISO{3}",
                "form": "iso_conduction", "args": {"k": k},
            },
        }

    g = Geometry(
        kind="layered_spheroid",
        args={"omega": 0.5, "radius": 1.0, "Nseries": 5},
        layers=[layer(0.3, 1.0), layer(0.7, 5.0)],
    )
    c = Cell(name="r", matrix_name="M", phases=[
        Phase(name="M", is_matrix=True, properties=[
            Property(key=":K", source="builder", builder="TensISO{3}",
                     form="iso_conduction", args={"k": 2.0})]),
        Phase(name="I", amount=0.2, geometry=g, properties=[]),
    ])
    return Model(cells=[c])


def test_layered_spheroid_orientation_reaches_the_generated_call():
    """A spheroid of revolution is orientable and the solver honors it —
    `scheme_integration.jl` returns `TensTI{2}(αt, αa, s.axis)`. The angles
    used to be dropped for this shape alone, leaving the axis at its (0,0,1)
    default: the interface's orientation fields did nothing, in the
    computation as much as in the 3-D view.

    `axis::Tuple` is the declared type, so the vector `vecbasis(...)[:, 3]`
    returns has to be wrapped — the trap the `TensTI{2}` builder documents.
    """
    m = _layered_spheroid_model()
    m.cells[0].phases[1].geometry.euler_angles = [0.7, 1.1]
    src = generate(m, embed_model=False)
    assert "axis = Tuple(vecbasis(RotatedBasis(0.7, 1.1" in src
    # …and no axis keyword at all when the shape is left unrotated.
    m2 = _layered_spheroid_model()
    assert "axis" not in generate(m2, embed_model=False)


def test_layered_spheroid_uses_the_fraction_constructor():
    """The raw constructor demands confocal layers, which typed-in radii are
    not: it threw, and the shape drew nothing."""
    src = generate(_layered_spheroid_model(), embed_model=False)
    assert "layered_spheroid_from_fractions(0.5, 1.0, (0.3, 0.7)" in src
    assert "LayeredSpheroid(" not in src


def test_conductivity_builder_has_the_right_arity():
    """`TensISO{dim}` — one argument is the 2nd-order form. `TensISO{2, 3}`
    named neither the right dimension nor the right order and threw."""
    src = generate(_layered_spheroid_model(), embed_model=False)
    assert "TensISO{3}(2.0)" in src
    assert "TensISO{2, 3}" not in src


def test_orientation_reaches_the_generated_call():
    c = Cell(name="r", matrix_name="M", phases=[
        Phase(name="M", is_matrix=True, properties=[Property()]),
        Phase(name="I", amount=0.2, geometry=Geometry(
            kind="spheroid", args={"omega": 0.3}, euler_angles=[0.7, 1.1])),
    ])
    src = generate(Model(cells=[c]), embed_model=False)
    assert "Spheroid(0.3; euler_angles = (0.7, 1.1))" in src


def test_a_single_angle_gets_the_tuple_comma():
    c = Cell(name="r", matrix_name="M", phases=[
        Phase(name="M", is_matrix=True, properties=[Property()]),
        Phase(name="I", amount=0.2, geometry=Geometry(
            kind="spheroid", args={"omega": 0.3}, euler_angles=[1.2])),
    ])
    assert "euler_angles = (1.2,)" in generate(Model(cells=[c]), embed_model=False)


def test_angles_are_floats_like_every_other_size():
    """A bare `0` next to `1.1` would make the tuple `Tuple{Int, Float64}`."""
    c = Cell(name="r", matrix_name="M", phases=[
        Phase(name="M", is_matrix=True, properties=[Property()]),
        Phase(name="I", amount=0.2, geometry=Geometry(
            kind="ellipsoid", args={"a": 2, "b": 1, "c": 0.5},
            euler_angles=[0, 1.1, 0.3])),
    ])
    src = generate(Model(cells=[c]), embed_model=False)
    assert "euler_angles = (0.0, 1.1, 0.3)" in src


def test_no_angles_means_no_keyword():
    src = generate(default_model(), embed_model=False)
    assert "euler_angles" not in src


def test_geometry_sizes_are_floats():
    """An NTuple mixing Int and Float64 fails to dispatch."""
    g = Geometry(kind="layered_sphere", layers=[
        {"radius": 0.6, "property": {"key": ":C", "source": "builder",
                                     "builder": "iso_stiffness", "args": {"k": 1, "mu": 1}}},
        {"radius": 1, "property": {"key": ":C", "source": "builder",
                                   "builder": "iso_stiffness", "args": {"k": 30, "mu": 12}}},
    ])
    c = Cell(name="r", matrix_name="M", phases=[
        Phase(name="M", is_matrix=True, properties=[Property()]),
        Phase(name="I", amount=1, geometry=g),
    ])
    src = generate(Model(cells=[c]), embed_model=False)
    assert "LayeredSphere((0.6, 1.0)" in src
    assert "(0.6, 1)" not in src


# ---------------------------------------------------------------------------
# Multiscale
# ---------------------------------------------------------------------------


def _two_scale() -> Model:
    inner = Cell(name="foam", matrix_name="SOLID", phases=[
        Phase(name="SOLID", is_matrix=True,
              properties=[Property(args={"k": 72.0, "mu": 32.0})]),
        Phase(name="PORE", amount=0.3,
              properties=[Property(args={"k": 1e-6, "mu": 1e-6})]),
    ])
    outer = Cell(name="paste", matrix_name="FOAM", phases=[
        Phase(name="FOAM", is_matrix=True, properties=[
            Property(key=":C", source="cell", cell=inner.id, scheme="SelfConsistent",
                     scheme_options={"abstol": 1e-10}),
        ]),
        Phase(name="CLINKER", amount=0.2,
              properties=[Property(args={"k": 100.0, "mu": 50.0})]),
    ])
    return Model(title="ms", cells=[inner, outer], root_cell=outer.id)


def test_seam_emits_homogenized():
    src = generate(_two_scale(), embed_model=False)
    assert "Homogenized(build_foam(), SelfConsistent(; abstol = 1.0e-10))" in src


def test_inner_scale_is_emitted_first():
    src = generate(_two_scale(), embed_model=False)
    assert src.index("function build_foam") < src.index("function build_paste")


def test_cycle_is_refused_at_construction():
    m = _two_scale()
    foam, paste = m.cells
    foam.phases[0].properties[0] = Property(key=":C", source="cell", cell=paste.id)
    problems = m.validate()
    assert any("cycle" in p for p in problems)
    assert "foam" in problems[0] and "paste" in problems[0]


def test_nested_lens_crosses_scales():
    m = _two_scale()
    m.sweep = Sweep(
        enabled=True, cell=m.cells[1].id,
        lens=Lens(kind="nested", member="FOAM", property=":C",
                  inner=Lens(kind="amount", phase="PORE").to_dict()),
    )
    src = generate(m, embed_model=False)
    assert "nested(:FOAM, :C, amount(:PORE))" in src


def test_alv_and_multiscale_are_refused_together():
    """MFH cannot re-express a homogenized inner result as a ViscoLaw."""
    m = _two_scale()
    m.alv.enabled = True
    assert any("viscoelast" in p.lower() for p in m.validate())


# ---------------------------------------------------------------------------
# Round-trip
# ---------------------------------------------------------------------------


def test_embedded_model_reopens_exactly():
    m = default_model()
    src = generate(m)
    back = extract_embedded(src)
    assert back is not None
    assert generate(Model.from_dict(back)) == src


def test_graph_positions_survive_the_round_trip():
    m = _two_scale()
    m.cells[0].ui = {"x": 123, "y": 456}
    back = Model.from_dict(extract_embedded(generate(m)))
    assert back.cells[0].ui == {"x": 123, "y": 456}


def test_untouched_parameter_keeps_its_original_text():
    m = default_model()
    m.params.append(Param(name="T", value="[1, 2]", origin="const T = [\n    1,\n    2,\n]"))
    src = generate(m, embed_model=False)
    assert "const T = [\n    1,\n    2,\n]" in src


def test_edited_parameter_is_regenerated():
    m = default_model()
    m.params.append(Param(name="T", value="99.0", origin="const T = 1.0", edited=True))
    assert "const T = 99.0" in generate(m, embed_model=False)


# ---------------------------------------------------------------------------
# Working without Julia
#
# The interface must come up whether or not the sidecar does. When it did not,
# `S.model` stayed null in the browser and every control threw a TypeError
# nobody sees — the whole thing looked broken with no clue why.
# ---------------------------------------------------------------------------


def test_catalog_is_complete_without_julia():
    from mfhstudio import catalog as catalog_module

    cat = catalog_module.base_catalog()
    assert cat["introspected"] is False
    for key in ("schemes", "geometries", "properties", "symmetrize",
                "projections", "interfaces", "lenses", "visco"):
        assert cat[key], f"{key} is empty without Julia"


def test_introspected_schemes_replace_the_fallback_wholesale():
    """A scheme MeanFieldHomogenization drops must disappear, not linger from a merge."""
    from mfhstudio import catalog as catalog_module

    merged = catalog_module.merge({
        "schemes": [{"name": "OnlyOne", "options": [], "singleton": True}],
        "mfh_version": "9.9.9", "julia_version": "1.x",
    })
    assert [s["name"] for s in merged["schemes"]] == ["OnlyOne"]
    assert merged["introspected"] is True
    assert merged["geometries"], "form definitions must survive the merge"


def test_session_serves_a_catalog_when_the_sidecar_is_dead():
    from mfhstudio.server import Session

    s = Session()
    s.bridge.julia = "/nonexistent-julia"
    cat = s.catalog()
    assert cat["introspected"] is False
    assert cat["schemes"] and cat["geometries"]
    assert s.catalog_error, "the failure must be reported, not swallowed"
    # and the model still generates a script
    assert "add_matrix!" in s.script()


def test_startup_failure_is_diagnosed_not_dumped():
    """A stack trace says what happened; the user needs to know what to do."""
    from mfhstudio.juliabridge import _diagnose

    log = ("ERROR: LoadError: ArgumentError: Package JSON3 [0f8b85d8] is "
           "required but does not seem to be installed:\n"
           " - Run `Pkg.instantiate()` to install all recorded dependencies.")
    msg = _diagnose(log)
    assert "instantiate" in msg
    assert "julia --project=" in msg
    assert log in msg, "the original error must still be there"


# ---------------------------------------------------------------------------
# Julia-backed
# ---------------------------------------------------------------------------


def _bridge():
    from mfhstudio.juliabridge import Bridge

    b = Bridge()
    b.start()
    return b


def test_catalog_covers_every_exported_scheme():
    """The interface must not fall behind MeanFieldHomogenization."""
    b = _bridge()
    try:
        cat = b.catalog()
        names = {s["name"] for s in cat["schemes"]}
        src = open(os.path.join(REPO, "src", "Schemes", "scheme_types.jl")).read()
        for expected in ("Voigt", "Reuss", "MoriTanaka", "SelfConsistent",
                         "AsymmetricSelfConsistent", "DifferentialScheme",
                         "Dilute", "DiluteDual", "Maxwell",
                         "PonteCastanedaWillis", "Laminated"):
            assert f"struct {expected}" in src or expected in src
            assert expected in names, f"{expected} missing from the catalog"
    finally:
        b.stop()


def test_self_consistent_offers_only_what_it_reads():
    """The kwargs bag accepts anything, so the option list must not come from
    probing the constructor."""
    b = _bridge()
    try:
        cat = b.catalog()
        sc = next(s for s in cat["schemes"] if s["name"] == "SelfConsistent")
        editable = {o["name"] for o in sc["options"] if o["editable"]}
        assert "abstol" in editable and "select_best" in editable
        assert "nsteps" not in editable, "nsteps is meaningless for SelfConsistent"
        diff = next(s for s in cat["schemes"] if s["name"] == "DifferentialScheme")
        deditable = {o["name"] for o in diff["options"] if o["editable"]}
        assert "nsteps" in deditable
        assert "select_best" not in deditable
    finally:
        b.stop()


def test_preserves_every_demo_script():
    """Open and save must not lose a line of anybody's script."""
    from mfhstudio.readback import model_from_script

    b = _bridge()
    try:
        files = sorted(glob.glob(os.path.join(REPO, "scripts", "*.jl")))
        assert files, "no demo scripts found"
        lossy = []
        for f in files:
            src = open(f, encoding="utf-8", errors="replace").read()
            model, _ = model_from_script(src, b)
            out = generate(model, embed_model=False)
            missing = [
                l.strip() for l in src.splitlines()
                if l.strip() and not l.strip().startswith("#") and l.strip() not in out
            ]
            if missing:
                lossy.append((os.path.basename(f), missing[:2]))
        assert not lossy, f"lost lines in {len(lossy)} script(s): {lossy[:3]}"
    finally:
        b.stop()


def test_generated_script_matches_the_echoes_reference():
    """The porous benchmark must reproduce the captured Echoes 1.0 values."""
    b = _bridge()
    try:
        m = default_model()
        m.sweep.start, m.sweep.stop, m.sweep.length = 0.3, 0.3, 2
        m.sweep.plot = False
        r = b.run(generate(m, embed_model=False), timeout=300)
        assert r["ok"], r.get("error")
        got = {}
        for line in r["stdout"].splitlines():
            mm = re.match(r"\s*MoriTanaka (\w+)\s+first = ([-\d.eE+]+)", line)
            if mm:
                got[mm.group(1)] = float(mm.group(2))
        assert got, r["stdout"]
        assert abs(got["k"] - 33.460582) < 1e-5, got
        assert abs(got["mu"] - 17.626742) < 1e-5, got
    finally:
        b.stop()


def test_traces_come_back_as_real_json():
    b = _bridge()
    try:
        sc = b.traces("Spheroid(0.4)")
        assert set(sc) == {"data", "layout"}
        assert sc["data"] and sc["data"][0]["type"] == "surface"
        # guides must not clutter the legend
        assert all(t.get("showlegend") is not True for t in sc["data"])
    finally:
        b.stop()


def test_a_tilted_layered_spheroid_is_drawn_tilted():
    """`LayeredSpheroid` stores a unit revolution axis, not a `.basis`, and
    `inclusion_basis` returns the canonical one whatever that axis is. The
    trace builder parametrized every layer with the axis hard-coded along z,
    so a tilted spheroid drew upright — the picture disagreeing with the
    script, exactly the failure the `_rot(::Ellipsoid)` comment records.
    """
    import math

    b = _bridge()
    try:
        mods = "(TensISO{3}(1.0), TensISO{3}(5.0))"
        base = f"layered_spheroid_from_fractions(0.5, 1.0, (0.3, 0.7), {mods}; Nseries = 5"
        th, ph = 0.9, 0.4
        axis = (
            math.sin(th) * math.cos(ph),
            math.sin(th) * math.sin(ph),
            math.cos(th),
        )
        upright = b.traces(base + ")")
        tilted = b.traces(base + f", axis = {axis})")

        def revolution_dir(scene):
            """Shortest principal direction of the outer layer's point cloud.

            ω = 0.5 is oblate, so the *short* semi-axis is the revolution one.
            """
            tr = scene["data"][-1]
            pts = [
                (x, y, z)
                for xs, ys, zs in zip(tr["x"], tr["y"], tr["z"])
                for x, y, z in zip(xs, ys, zs)
            ]
            n = len(pts)
            cov = [[sum(p[i] * p[j] for p in pts) / n for j in range(3)] for i in range(3)]
            # Power iteration on (tr(C)·I − C) converges to C's *smallest*
            # eigenvector, avoiding a numpy dependency in the test suite.
            tr_c = sum(cov[i][i] for i in range(3))
            v = [0.3, 0.5, 0.81]
            for _ in range(400):
                w = [
                    sum((tr_c * (i == j) - cov[i][j]) * v[j] for j in range(3))
                    for i in range(3)
                ]
                nrm = math.sqrt(sum(c * c for c in w)) or 1.0
                v = [c / nrm for c in w]
            return v

        d_up = revolution_dir(upright)
        d_tl = revolution_dir(tilted)
        # A principal direction has no sign, so compare |cos|.
        assert abs(d_up[2]) > 0.99, d_up
        assert abs(sum(a * b_ for a, b_ in zip(d_tl, axis))) > 0.99, d_tl
        assert abs(d_tl[2]) < 0.9, d_tl
    finally:
        b.stop()


JULIA_TESTS = {
    "test_catalog_covers_every_exported_scheme",
    "test_self_consistent_offers_only_what_it_reads",
    "test_preserves_every_demo_script",
    "test_generated_script_matches_the_echoes_reference",
    "test_traces_come_back_as_real_json",
    "test_a_tilted_layered_spheroid_is_drawn_tilted",
}


# ---------------------------------------------------------------------------

if __name__ == "__main__":
    fns = [(k, v) for k, v in sorted(globals().items()) if k.startswith("test_")]
    failed = skipped = 0
    for name, fn in fns:
        if name in JULIA_TESTS and not WITH_JULIA:
            skipped += 1
            print(f"  skip  {name}  (pass --julia to run)")
            continue
        try:
            fn()
            print(f"  ok    {name}")
        except AssertionError as e:
            failed += 1
            print(f"  FAIL  {name}  {e}")
        except Exception as e:  # noqa: BLE001
            failed += 1
            print(f"  ERROR {name}  {type(e).__name__}: {e}")
    total = len(fns) - skipped
    print(f"\n{total - failed}/{total} passed" + (f", {skipped} skipped" if skipped else ""))
    sys.exit(1 if failed else 0)
