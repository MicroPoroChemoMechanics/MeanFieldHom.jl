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
    m.sweep.scheme = "SelfConsistent"
    m.sweep.scheme_options = {"abstol": 1e-10, "maxiters": 300, "select_best": True}
    src = generate(m, embed_model=False)
    assert "SelfConsistent(; abstol = 1.0e-10, maxiters = 300, select_best = true)" in src
    assert "homogenize(cell, scheme, :C)" in src


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
# Julia-backed
# ---------------------------------------------------------------------------


def _bridge():
    from mfhstudio.juliabridge import Bridge

    b = Bridge()
    b.start()
    return b


def test_catalogue_covers_every_exported_scheme():
    """The interface must not fall behind MeanFieldHom."""
    b = _bridge()
    try:
        cat = b.catalogue()
        names = {s["name"] for s in cat["schemes"]}
        src = open(os.path.join(REPO, "src", "Schemes", "scheme_types.jl")).read()
        for expected in ("Voigt", "Reuss", "MoriTanaka", "SelfConsistent",
                         "AsymmetricSelfConsistent", "DifferentialScheme",
                         "Dilute", "DiluteDual", "Maxwell",
                         "PonteCastanedaWillis", "Laminated"):
            assert f"struct {expected}" in src or expected in src
            assert expected in names, f"{expected} missing from the catalogue"
    finally:
        b.stop()


def test_self_consistent_offers_only_what_it_reads():
    """The kwargs bag accepts anything, so the option list must not come from
    probing the constructor."""
    b = _bridge()
    try:
        cat = b.catalogue()
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
        r = b.run(generate(m, embed_model=False), timeout=300)
        assert r["ok"], r.get("error")
        rows = [l.split() for l in r["stdout"].splitlines() if l.strip().startswith("0.3")]
        assert rows, r["stdout"]
        k, mu = float(rows[0][1]), float(rows[0][2])
        assert abs(k - 33.460582) < 1e-5, k
        assert abs(mu - 17.626742) < 1e-5, mu
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


JULIA_TESTS = {
    "test_catalogue_covers_every_exported_scheme",
    "test_self_consistent_offers_only_what_it_reads",
    "test_preserves_every_demo_script",
    "test_generated_script_matches_the_echoes_reference",
    "test_traces_come_back_as_real_json",
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
