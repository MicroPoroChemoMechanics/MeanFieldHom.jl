"""The authoring model.

This is what the interface edits and what the code generator reads. It is a
*graph of cells*, not a single RVE, because MeanFieldHom chains scales
declaratively: a phase property may hold a `Homogenized(inner_cell, scheme)`
instead of a tensor, and the outer scheme resolves the inner scale when it
reads the property.

Everything is plain dataclasses with `to_dict`/`from_dict`, so the whole model
serializes to JSON — for the browser, for the embedded round-trip block, and
for the test fixtures.
"""

from __future__ import annotations

import uuid
from dataclasses import dataclass, field, asdict
from typing import Any, Optional


def _new_id() -> str:
    return uuid.uuid4().hex[:8]


# ---------------------------------------------------------------------------
# Values: either a literal, a named parameter, or a nested scale
# ---------------------------------------------------------------------------


@dataclass
class Property:
    """One entry of a phase's property dictionary.

    `source` selects where the value comes from:

    - ``"builder"``  — a tensor built from moduli, e.g. ``iso_stiffness(k, μ)``
    - ``"expr"``     — a raw Julia expression the user typed
    - ``"cell"``     — **the multiscale seam**: the value is the effective
      property of another cell, emitted as ``Homogenized(cell, scheme)``
    """

    key: str = ":C"
    source: str = "builder"
    builder: str = "iso_stiffness"
    #: which catalog entry produced this. Several entries share one builder
    #: (`iso_stiffness` backs both the plain isotropic form and the near-zero
    #: pore preset), so the builder alone cannot identify the form.
    form: str = "iso_kmu"
    args: dict = field(default_factory=lambda: {"k": 10.0, "mu": 5.0})
    expr: str = ""
    # multiscale seam
    cell: Optional[str] = None
    scheme: Optional[str] = None
    scheme_options: dict = field(default_factory=dict)
    # viscoelastic law instead of a tensor
    visco: Optional[dict] = None

    def to_dict(self) -> dict:
        return asdict(self)

    @staticmethod
    def from_dict(d: dict) -> "Property":
        return Property(**{k: v for k, v in d.items() if k in Property.__annotations__})


@dataclass
class Geometry:
    kind: str = "spheroid"
    args: dict = field(default_factory=lambda: {"omega": 1.0})
    euler_angles: list = field(default_factory=list)
    #: layered inclusions only: list of {radius, property, interface}
    layers: list = field(default_factory=list)

    def to_dict(self) -> dict:
        return asdict(self)

    @staticmethod
    def from_dict(d: dict) -> "Geometry":
        return Geometry(**{k: v for k, v in d.items() if k in Geometry.__annotations__})


@dataclass
class Phase:
    name: str = "PHASE"
    is_matrix: bool = False
    geometry: Geometry = field(default_factory=Geometry)
    properties: list = field(default_factory=list)  # list[Property]
    #: ("fraction" | "density", value-or-parameter-name); ignored for the matrix,
    #: whose amount MFH derives as 1 - Σ f_inclusions and refuses to be set.
    amount_kind: str = "fraction"
    amount: Any = 0.1
    symmetrize: str = "none"

    def to_dict(self) -> dict:
        d = asdict(self)
        d["geometry"] = self.geometry.to_dict()
        d["properties"] = [
            p.to_dict() if isinstance(p, Property) else p for p in self.properties
        ]
        return d

    @staticmethod
    def from_dict(d: dict) -> "Phase":
        return Phase(
            name=d.get("name", "PHASE"),
            is_matrix=bool(d.get("is_matrix", False)),
            geometry=Geometry.from_dict(d.get("geometry", {})),
            properties=[Property.from_dict(p) for p in d.get("properties", [])],
            amount_kind=d.get("amount_kind", "fraction"),
            amount=d.get("amount", 0.1),
            symmetrize=d.get("symmetrize", "none"),
        )


@dataclass
class Cell:
    """One scale: an RVE with a matrix and inclusion phases."""

    id: str = field(default_factory=_new_id)
    name: str = "rve"
    matrix_name: str = "MATRIX"
    phases: list = field(default_factory=list)  # list[Phase]
    #: parameters this cell's builder takes (discovered from the sweep)
    params: list = field(default_factory=list)
    #: the builder's function name. Read-back keeps whatever the file used, so
    #: a script whose builder is called `_ec_equiv` keeps that name and its
    #: callers keep working; only new cells get the `build_<name>` convention.
    builder_name: Optional[str] = None
    #: extra keywords on the RVE constructor, e.g. `T = ComplexF64`
    rve_options: dict = field(default_factory=dict)
    #: where the user dragged this scale in the graph view. Purely cosmetic,
    #: but it rides along in the embedded model so a reopened file looks the
    #: way it was left.
    ui: dict = field(default_factory=lambda: {"x": 40, "y": 40})

    @property
    def builder(self) -> str:
        return self.builder_name or f"build_{self.name}"

    def matrix(self) -> Optional[Phase]:
        return next((p for p in self.phases if p.is_matrix), None)

    def inclusions(self) -> list:
        return [p for p in self.phases if not p.is_matrix]

    def to_dict(self) -> dict:
        return {
            "id": self.id,
            "name": self.name,
            "matrix_name": self.matrix_name,
            "phases": [p.to_dict() for p in self.phases],
            "params": list(self.params),
            "builder_name": self.builder_name,
            "rve_options": dict(self.rve_options),
            "ui": dict(self.ui),
        }

    @staticmethod
    def from_dict(d: dict) -> "Cell":
        return Cell(
            id=d.get("id") or _new_id(),
            name=d.get("name", "rve"),
            matrix_name=d.get("matrix_name", "MATRIX"),
            phases=[Phase.from_dict(p) for p in d.get("phases", [])],
            params=list(d.get("params", [])),
            builder_name=d.get("builder_name"),
            rve_options=dict(d.get("rve_options", {})),
            ui=dict(d.get("ui") or {"x": 40, "y": 40}),
        )


# ---------------------------------------------------------------------------
# Parameters, sweeps, outputs
# ---------------------------------------------------------------------------


@dataclass
class Param:
    """A named constant emitted as `const name = value`.

    `origin` holds the exact text the parameter had in the file it was read
    from. As long as the user has not edited it, that text is what gets written
    back — so a carefully laid-out multi-line constant survives instead of
    being collapsed into one line by the AST round-trip. Re-formatting code the
    interface did not author is a form of damage, even when the meaning is
    preserved.
    """

    name: str = "k"
    value: str = "1.0"
    comment: str = ""
    origin: Optional[str] = None
    edited: bool = False

    def to_dict(self) -> dict:
        return asdict(self)

    @staticmethod
    def from_dict(d: dict) -> "Param":
        return Param(**{k: v for k, v in d.items() if k in Param.__annotations__})


@dataclass
class Lens:
    """A sensitivity/sweep lens, mirroring `src/Schemes/parameters.jl`.

    `nested` wraps another lens to reach into an inner scale, which is how a
    sweep crosses scales without any hand-written closure.
    """

    kind: str = "amount"  # amount | property | geometry | shape_param | nested
    phase: str = ""
    property: str = ":C"
    field_name: str = "semi_axes"
    index: int = 1
    member: str = ""
    inner: Optional[dict] = None

    def to_dict(self) -> dict:
        return asdict(self)

    @staticmethod
    def from_dict(d: dict) -> "Lens":
        return Lens(**{k: v for k, v in d.items() if k in Lens.__annotations__})


@dataclass
class Sweep:
    """A parametric sweep over one lens, plotted."""

    enabled: bool = False
    variable: str = "φ"
    start: float = 0.0
    stop: float = 1.0
    length: int = 21
    lens: Lens = field(default_factory=Lens)
    #: which cell the sweep drives
    cell: Optional[str] = None
    scheme: str = "MoriTanaka"
    scheme_options: dict = field(default_factory=dict)
    property: str = ":C"
    #: reporting projection applied to the result: none | iso | ti | ortho
    projection: str = "none"
    #: which scalars to extract and plot
    outputs: list = field(default_factory=lambda: ["k", "mu"])
    plot: bool = True

    def to_dict(self) -> dict:
        d = asdict(self)
        d["lens"] = self.lens.to_dict()
        return d

    @staticmethod
    def from_dict(d: dict) -> "Sweep":
        s = Sweep(
            **{
                k: v
                for k, v in d.items()
                if k in Sweep.__annotations__ and k != "lens"
            }
        )
        s.lens = Lens.from_dict(d.get("lens", {}))
        return s


@dataclass
class Alv:
    """Ageing linear viscoelasticity settings."""

    enabled: bool = False
    t_start: float = 0.0
    t_stop: float = 10.0
    length: int = 41
    log_time: bool = False
    cell: Optional[str] = None
    scheme: str = "MoriTanaka"
    property: str = ":C"

    def to_dict(self) -> dict:
        return asdict(self)

    @staticmethod
    def from_dict(d: dict) -> "Alv":
        return Alv(**{k: v for k, v in d.items() if k in Alv.__annotations__})


@dataclass
class OpaqueBlock:
    """Source the studio did not recognize, preserved verbatim.

    This is the whole reason the interface is safe to point at a hand-written
    script: what it does not understand, it does not touch. The block is
    re-emitted byte for byte and shown read-only in the UI.
    """

    source: str = ""
    #: where it sat in the original file, so ordering survives
    order: int = 0
    note: str = ""

    def to_dict(self) -> dict:
        return asdict(self)

    @staticmethod
    def from_dict(d: dict) -> "OpaqueBlock":
        return OpaqueBlock(
            **{k: v for k, v in d.items() if k in OpaqueBlock.__annotations__}
        )


# ---------------------------------------------------------------------------
# The whole model
# ---------------------------------------------------------------------------


@dataclass
class Model:
    title: str = "model"
    description: str = ""
    params: list = field(default_factory=list)   # list[Param]
    cells: list = field(default_factory=list)    # list[Cell]
    sweep: Sweep = field(default_factory=Sweep)
    alv: Alv = field(default_factory=Alv)
    opaque: list = field(default_factory=list)   # list[OpaqueBlock]
    #: id of the cell that carries the final result
    root_cell: Optional[str] = None

    # -- lookups ----------------------------------------------------------

    def cell(self, cid: Optional[str]) -> Optional[Cell]:
        if cid is None:
            return None
        return next((c for c in self.cells if c.id == cid), None)

    def cell_by_name(self, name: str) -> Optional[Cell]:
        return next((c for c in self.cells if c.name == name), None)

    def root(self) -> Optional[Cell]:
        return self.cell(self.root_cell) or (self.cells[-1] if self.cells else None)

    # -- the multiscale graph --------------------------------------------

    def dependencies(self, cell: Cell) -> list:
        """The cells this one reads through a `Homogenized` seam."""
        out = []
        for ph in cell.phases:
            for pr in ph.properties:
                if pr.source == "cell" and pr.cell:
                    out.append(pr.cell)
        return out

    def topological_order(self) -> list:
        """Cells ordered so every cell comes after the ones it depends on.

        Raises `ValueError` naming the cycle when the graph has one — a scale
        cannot be built from itself, and catching it here means the interface
        refuses at construction rather than emitting a script that recurses
        forever.
        """
        state: dict = {}
        order: list = []

        def visit(cid: str, trail: list) -> None:
            st = state.get(cid)
            if st == "done":
                return
            if st == "visiting":
                names = [self.cell(x).name if self.cell(x) else x for x in trail + [cid]]
                start = names.index(names[-1])
                raise ValueError(
                    "multiscale cycle: " + " → ".join(names[start:])
                )
            c = self.cell(cid)
            if c is None:
                return
            state[cid] = "visiting"
            for dep in self.dependencies(c):
                visit(dep, trail + [cid])
            state[cid] = "done"
            order.append(c)

        for c in self.cells:
            visit(c.id, [])
        return order

    def uses_multiscale(self) -> bool:
        return any(self.dependencies(c) for c in self.cells)

    def validate(self) -> list:
        """Problems worth blocking on, each as a human-readable string."""
        problems = []
        try:
            self.topological_order()
        except ValueError as exc:
            problems.append(str(exc))

        names = [c.name for c in self.cells]
        for n in set(names):
            if names.count(n) > 1:
                problems.append(f"two cells are both named `{n}`")

        for c in self.cells:
            if c.matrix() is None:
                problems.append(f"cell `{c.name}` has no matrix phase")
            pn = [p.name for p in c.phases]
            for n in set(pn):
                if pn.count(n) > 1:
                    problems.append(f"cell `{c.name}` has two phases named `{n}`")

        # Documented MFH constraint: an inner Homogenized cannot sit inside an
        # ageing-viscoelastic chain, because the inner result would have to be
        # re-expressible as a ViscoLaw (src/Core/cells.jl).
        if self.alv.enabled and self.uses_multiscale():
            problems.append(
                "ageing viscoelasticity cannot be combined with a nested scale: "
                "MeanFieldHom cannot re-express a homogenized inner result as a "
                "ViscoLaw"
            )
        return problems

    # -- serialization ----------------------------------------------------

    def to_dict(self) -> dict:
        return {
            "title": self.title,
            "description": self.description,
            "params": [p.to_dict() for p in self.params],
            "cells": [c.to_dict() for c in self.cells],
            "sweep": self.sweep.to_dict(),
            "alv": self.alv.to_dict(),
            "opaque": [o.to_dict() for o in self.opaque],
            "root_cell": self.root_cell,
        }

    @staticmethod
    def from_dict(d: dict) -> "Model":
        m = Model(
            title=d.get("title", "model"),
            description=d.get("description", ""),
            params=[Param.from_dict(p) for p in d.get("params", [])],
            cells=[Cell.from_dict(c) for c in d.get("cells", [])],
            opaque=[OpaqueBlock.from_dict(o) for o in d.get("opaque", [])],
            root_cell=d.get("root_cell"),
        )
        m.sweep = Sweep.from_dict(d.get("sweep", {}))
        m.alv = Alv.from_dict(d.get("alv", {}))
        return m


def default_model() -> Model:
    """The porous benchmark, which is the shortest useful thing to open on."""
    solid = Phase(
        name="SOLID", is_matrix=True,
        geometry=Geometry(kind="spheroid", args={"omega": 1.0}),
        properties=[
            Property(key=":C", builder="iso_stiffness", form="iso_kmu",
                     args={"k": 72.0, "mu": 32.0})
        ],
    )
    pore = Phase(
        name="PORE", is_matrix=False,
        geometry=Geometry(kind="spheroid", args={"omega": 1.0}),
        properties=[
            Property(key=":C", builder="iso_stiffness", form="void",
                     args={"k": 1.0e-6, "mu": 1.0e-6})
        ],
        amount_kind="fraction", amount=0.1,
    )
    cell = Cell(name="rve", matrix_name="SOLID", phases=[solid, pore])
    m = Model(title="porous_benchmark", cells=[cell], root_cell=cell.id)
    m.sweep = Sweep(
        enabled=True, variable="φ", start=0.0, stop=0.9, length=19,
        lens=Lens(kind="amount", phase="PORE"), cell=cell.id,
        scheme="MoriTanaka", projection="iso", outputs=["k", "mu"],
    )
    return m
