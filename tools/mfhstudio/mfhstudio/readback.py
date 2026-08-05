"""Reading an existing script back into the model.

Two levels, and the second one carries the guarantee that matters:

1. A script the studio wrote embeds its model, so it reopens exactly.
2. Any other script is parsed by the Julia sidecar (`Meta.parse`, the real
   parser) and matched against the vocabulary the generator emits. Everything
   else is kept as an **opaque block, verbatim**.

The second rule is the contract: what the studio does not understand, it does
not touch. A script somebody hand-tuned can be opened, edited in the parts the
interface does understand, and written back with the rest byte-identical.
"""

from __future__ import annotations

import re
from typing import Optional, Tuple

from .codegen import extract_embedded
from .model import (
    Cell,
    Geometry,
    Model,
    OpaqueBlock,
    Param,
    Phase,
    Property,
)

# ---------------------------------------------------------------------------
# Julia expression → the form fields the interface edits
# ---------------------------------------------------------------------------

_GEOM_PATTERNS = [
    ("spheroid", re.compile(r"^Spheroid\(\s*([^;,)]+)")),
    ("ellipsoid", re.compile(r"^Ellipsoid\(\s*([^;,)]+),\s*([^;,)]+),\s*([^;,)]+)")),
    ("cylinder", re.compile(r"^Cylinder\(\s*([^;,)]+),\s*([^;,)]+)")),
    ("penny_crack", re.compile(r"^PennyCrack\(\s*([^;,)]+)")),
    ("elliptic_crack", re.compile(r"^EllipticCrack\(\s*([^;,)]+),\s*([^;,)]+)")),
    ("ribbon_crack", re.compile(r"^RibbonCrack\(\s*([^;,)]+)")),
]

_GEOM_FIELDS = {
    "spheroid": ["omega"],
    "ellipsoid": ["a", "b", "c"],
    "cylinder": ["b", "c"],
    "penny_crack": ["a"],
    "elliptic_crack": ["a", "b"],
    "ribbon_crack": ["b"],
}


def _number(text: str):
    """A literal becomes a float; anything else stays a Julia expression."""
    t = text.strip()
    try:
        return float(t)
    except ValueError:
        return t


def parse_geometry(expr: str) -> Optional[Geometry]:
    e = expr.strip()
    for kind, pat in _GEOM_PATTERNS:
        m = pat.match(e)
        if m is None:
            continue
        names = _GEOM_FIELDS[kind]
        args = {n: _number(v) for n, v in zip(names, m.groups())}
        g = Geometry(kind=kind, args=args)
        ang = re.search(r"euler_angles\s*=\s*\(([^)]*)\)", e)
        if ang:
            g.euler_angles = [
                _number(x) for x in ang.group(1).split(",") if x.strip()
            ]
        return g
    return None


_PROP_PATTERNS = [
    ("iso_stiffness", ["k", "mu"], re.compile(r"^iso_stiffness\(\s*([^,)]+),\s*([^,)]+)\)")),
    (
        "iso_stiffness_E_nu", ["E", "nu"],
        re.compile(r"^iso_stiffness_E_nu\(\s*([^,)]+),\s*([^,)]+)\)"),
    ),
    (
        "hoenig_stiffness", ["E1", "h", "nu1", "nu2", "gamma"],
        re.compile(
            r"^hoenig_stiffness\(\s*([^,)]+),\s*([^,)]+),\s*([^,)]+),\s*([^,)]+),\s*([^,)]+)\)"
        ),
    ),
]


def parse_property(key: str, expr: str) -> Property:
    e = expr.strip()

    m = re.match(r"^Homogenized\(\s*build_([A-Za-z_][\w]*)\s*\(", e)
    if m:
        # The multiscale seam. The inner cell is resolved by name afterwards,
        # once every cell has been read.
        scheme = "MoriTanaka"
        sm = re.search(r",\s*([A-Za-z][\w]*)\s*\(", e[m.end():])
        if sm:
            scheme = sm.group(1)
        p = Property(key=key, source="cell", scheme=scheme)
        p.cell = "@name:" + m.group(1)
        return p

    for builder, names, pat in _PROP_PATTERNS:
        mm = pat.match(e)
        if mm:
            return Property(
                key=key, source="builder", builder=builder,
                args={n: _number(v) for n, v in zip(names, mm.groups())},
            )

    # Anything else is kept as the expression the author wrote.
    return Property(key=key, source="expr", expr=e)


def parse_properties(expr: str) -> list:
    """`Dict(:C => iso_stiffness(72.0, 32.0), :K => …)` → a list of Property."""
    e = expr.strip()
    if not e.startswith("Dict"):
        return []
    inner = e[e.find("(") + 1: e.rfind(")")]
    out = []
    depth = 0
    start = 0
    parts = []
    for i, ch in enumerate(inner):
        if ch in "([{":
            depth += 1
        elif ch in ")]}":
            depth -= 1
        elif ch == "," and depth == 0:
            parts.append(inner[start:i])
            start = i + 1
    parts.append(inner[start:])
    for part in parts:
        if "=>" not in part:
            continue
        key, _, val = part.partition("=>")
        key = key.strip()
        if not key.startswith(":"):
            continue
        out.append(parse_property(key, val))
    return out


# ---------------------------------------------------------------------------
# Nodes → model
# ---------------------------------------------------------------------------


def model_from_script(source: str, bridge) -> Tuple[Model, dict]:
    """Read a script. Returns the model and a report on what was understood."""

    embedded = extract_embedded(source)
    if embedded is not None:
        return Model.from_dict(embedded), {
            "exact": True,
            "recognized": None,
            "opaque": 0,
            "note": "reopened from the model the studio embedded in the file",
        }

    try:
        parsed = bridge.parse(source)
    except Exception as exc:  # noqa: BLE001
        # Without the sidecar nothing can be read safely, and guessing with a
        # regex would be exactly the silent-corruption risk this design avoids.
        return Model(), {
            "exact": False, "error": str(exc),
            "note": "the Julia sidecar is required to read a script the studio "
                    "did not write",
        }

    if "error" in parsed:
        return Model(), {"exact": False, "error": parsed["error"]}

    model = Model(title="opened")
    cells_by_name: dict = {}
    order = 0

    for node in parsed.get("nodes", []):
        kind = node.get("kind")
        src = node.get("source", "")

        if kind == "trivia":
            # Blank lines and comments between statements: keep them, but a
            # run of pure whitespace carries nothing worth restoring.
            if src.strip():
                model.opaque.append(OpaqueBlock(source=src, order=order))
                order += 1
            continue

        if kind == "param":
            model.params.append(
                Param(name=node["name"], value=node["value"], origin=src)
            )
            order += 1
            continue

        if kind == "rve_builder":
            why = _unrepresentable(node)
            if why is None:
                why = _fails_to_reproduce(node, src)
            if why is not None:
                # Recognizing a builder the model cannot reproduce faithfully
                # is worse than not recognizing it: regenerating would change
                # the function's signature and break its callers. Keep it.
                model.opaque.append(
                    OpaqueBlock(
                        source=src, order=order,
                        note=f"kept as written: {why}",
                    )
                )
                order += 1
                continue
            cell = _cell_from_node(node)
            model.cells.append(cell)
            cells_by_name[cell.name] = cell
            order += 1
            continue

        # `homogenize`, and everything else, is preserved as written. The
        # interface will not rewrite a result section it cannot fully model.
        model.opaque.append(OpaqueBlock(source=src, order=order))
        order += 1

    # Resolve the multiscale seams now that every cell has a name.
    unresolved = []
    for c in model.cells:
        for ph in c.phases:
            for pr in ph.properties:
                if pr.source == "cell" and isinstance(pr.cell, str) and pr.cell.startswith("@name:"):
                    name = pr.cell[len("@name:"):]
                    target = cells_by_name.get(name)
                    if target is None:
                        unresolved.append(name)
                        pr.source = "expr"
                        pr.expr = f"#= unresolved inner cell `{name}` =#"
                    else:
                        pr.cell = target.id

    if model.cells:
        model.root_cell = model.cells[-1].id

    report = {
        "exact": False,
        "recognized": parsed.get("recognized", 0),
        "opaque": len(model.opaque),
        "cells": [c.name for c in model.cells],
        "unresolved_cells": unresolved,
        "note": (
            "the parts shown as preserved were not recognized and will be "
            "written back unchanged"
        ),
    }
    return model, report


def _fails_to_reproduce(node: dict, original: str) -> Optional[str]:
    """Would writing this cell back change the file?

    The structural checks below catch the cases we know about; this one
    catches the ones we do not. The cell is rendered exactly as it would be
    saved and compared against the source it came from — if the two differ,
    the studio does not really understand it, whatever the parser said, and
    the original text is kept.

    Comparison ignores indentation and blank lines only: those the generator
    normalizes on purpose, and they carry no meaning in Julia.
    """
    from .codegen import render_cell

    def norm(text: str) -> list:
        return [ln.strip() for ln in text.splitlines() if ln.strip()]

    try:
        cell = _cell_from_node(node)
        rendered = render_cell(Model(cells=[cell]), cell)
    except Exception as exc:  # noqa: BLE001
        return f"rebuilding it raised {type(exc).__name__}"

    if norm(rendered) != norm(original):
        return "writing it back would not reproduce the original text"
    return None


_IDENT = re.compile(r"^[A-Za-z_À-￿][\wÀ-￿!]*$")


def _unrepresentable(node: dict) -> Optional[str]:
    """Why this builder cannot be regenerated faithfully, or None if it can.

    The model describes a builder as a plain list of positional parameter
    names. A signature with type annotations, defaults or keyword arguments
    carries meaning the model would drop, and its callers depend on that
    meaning.
    """
    for p in node.get("params") or []:
        s = str(p).strip()
        if not _IDENT.match(s):
            return f"the builder signature uses `{s}`, which the studio does not model"

    for ph in node.get("phases") or []:
        geom = ph.get("geometry", "")
        if parse_geometry(geom) is None:
            return f"geometry `{geom}` is not one the studio can rebuild"
        props = ph.get("properties", "")
        if props.strip() and not props.strip().startswith("Dict"):
            return f"properties `{props}` are not a Dict the studio can rebuild"
        for opt, val in (ph.get("options") or {}).items():
            if opt == "symmetrize" and not (
                "IsoSymmetrize" in val or "TISymmetrize" in val or val.strip() == "nothing"
            ):
                return f"symmetrize = `{val}` is not one the studio can rebuild"
    return None


def _cell_from_node(node: dict) -> Cell:
    fname = node.get("function", "build_rve")
    name = fname[len("build_"):] if fname.startswith("build_") else fname
    cell = Cell(
        name=name,
        matrix_name=node.get("matrix") or "MATRIX",
        params=list(node.get("params") or []),
        # Keep the function's own name: regenerating it as `build_<name>`
        # would rename a function its callers refer to.
        builder_name=fname,
        rve_options=dict(node.get("rve_options") or {}),
    )
    for p in node.get("phases", []):
        geom = parse_geometry(p.get("geometry", "")) or Geometry()
        ph = Phase(
            name=p.get("name", "PHASE"),
            is_matrix=bool(p.get("is_matrix")),
            geometry=geom,
            properties=parse_properties(p.get("properties", "")),
        )
        opts = p.get("options") or {}
        if "fraction" in opts:
            ph.amount_kind, ph.amount = "fraction", _number(opts["fraction"])
        elif "density" in opts:
            ph.amount_kind, ph.amount = "density", _number(opts["density"])
        sym = opts.get("symmetrize", "")
        if "IsoSymmetrize" in sym:
            ph.symmetrize = "iso"
        elif "TISymmetrize" in sym:
            ph.symmetrize = "ti"
        cell.phases.append(ph)
    return cell
