"""Model → Julia.

The conventions here are the ones `tools/echoes2mfh/echoes2mfh/emit.py` already
establishes, and the two generators are kept deliberately consistent: a script
written by the studio and one translated from Echoes should read the same way.

Three of those conventions are load-bearing and worth restating, because the
whole point of the interface is that the user never has to remember them:

* `iso_stiffness(k, μ)` takes *physical* moduli, while the raw `TensISO{3}(a, b)`
  constructor takes `(3k, 2μ)`. Only the former is ever emitted.
* The matrix amount is derived (`1 − Σ f`) and assigning it raises. The matrix
  phase therefore has no amount in the generated code.
* Solver options belong to the **scheme instance**, not to `homogenize`.
"""

from __future__ import annotations

import json
from typing import Optional

from .model import Cell, Model, Phase, Property

IND = "    "
RULE = "# " + "=" * 77
MODEL_TAG = "mfhstudio-model v1"


def _rule(title: str) -> str:
    bar = "─" * max(3, 72 - len(title))
    return f"# ── {title} {bar}"


def _tuple(items) -> str:
    """A Julia tuple; only a 1-tuple needs the trailing comma."""
    parts = [p for p in items if p]
    if not parts:
        return "()"
    if len(parts) == 1:
        return f"({parts[0]},)"
    return "(" + ", ".join(parts) + ")"


def _fnum(x) -> str:
    """Like `_num`, but never yields a bare integer.

    Geometry constructors take an `NTuple{N, T}`: mixing `1` and `0.6` gives a
    `Tuple{Float64, Int64}` and no method matches. Sizes are therefore always
    written as floats, while genuine counts (`nsteps`, `length`) keep `_num`.
    """
    if isinstance(x, str):
        return x
    if isinstance(x, bool):
        return _num(x)
    if isinstance(x, int):
        return _num(float(x))
    return _num(x)


def _num(x) -> str:
    """Render a value as Julia, keeping float-ness explicit."""
    if isinstance(x, bool):
        return "true" if x else "false"
    if isinstance(x, int):
        return str(x)
    if isinstance(x, float):
        s = repr(x)
        if "e" in s or "E" in s:
            mant, _, exp = s.partition("e")
            if "." not in mant:
                mant += ".0"
            return f"{mant}e{int(exp)}"
        return s if "." in s else s + ".0"
    return str(x)


# The activation preamble is shared with echoes2mfh: a generated script must
# run from wherever it is saved, not only from `scripts/`.
_ACTIVATE = """import Pkg
let d = @__DIR__
    while true
        pt = joinpath(d, "Project.toml")
        if isfile(pt) && occursin("MeanFieldHom", read(pt, String))
            Pkg.activate(d; io = devnull)
            break
        end
        parent = dirname(d)
        parent == d && break
        d = parent
    end
end"""


class CodeGen:
    def __init__(self, model: Model, embed_model: bool = True):
        self.m = model
        self.embed = embed_model
        self.out: list = []

    # -- plumbing ---------------------------------------------------------

    def w(self, line: str = "") -> None:
        self.out.append(line)

    def blank(self) -> None:
        if self.out and self.out[-1] != "":
            self.out.append("")

    # -- entry point ------------------------------------------------------

    def generate(self) -> str:
        problems = self.m.validate()
        self._header(problems)
        self._preamble()
        self._params()
        self._builders()
        self._opaque()
        self._main()
        self._embedded_model()
        return "\n".join(self.out).rstrip() + "\n"

    # -- sections ---------------------------------------------------------

    def _header(self, problems: list) -> None:
        self.w(RULE)
        self.w(f"#  {self.m.title}.jl")
        self.w("#")
        self.w("#  Built with MFH Studio. Editing this file by hand is fine: the")
        self.w("#  studio reads it back and preserves anything it does not")
        self.w("#  recognize.")
        if self.m.description:
            self.w("#")
            for line in self.m.description.splitlines():
                self.w(f"#  {line}")
        if problems:
            self.w("#")
            self.w("#  UNRESOLVED PROBLEMS — this script will not run as-is:")
            for p in problems:
                self.w(f"#    * {p}")
        self.w(RULE)
        self.w()

    def _preamble(self) -> None:
        for line in _ACTIVATE.splitlines():
            self.w(line)
        self.blank()
        self.w("using MeanFieldHom")
        self.w("using TensND")
        self.w("using LinearAlgebra")
        self.w("using Printf")
        if self.m.sweep.enabled and self.m.sweep.plot:
            self.w("using Plots")
            self.w("gr()")
        self.blank()

    def _params(self) -> None:
        if not self.m.params:
            return
        self.w(_rule("Parameters"))
        for p in self.m.params:
            # A parameter read from a file keeps its original text until the
            # user actually changes it: re-printing from the AST would collapse
            # deliberate multi-line layout for no benefit.
            if p.origin and not p.edited:
                for line in p.origin.splitlines():
                    self.w(line)
                continue
            line = f"const {p.name} = {p.value}"
            if p.comment:
                line += f"  # {p.comment}"
            self.w(line)
        self.blank()

    # -- cells ------------------------------------------------------------

    def _builders(self) -> None:
        if not self.m.cells:
            return
        try:
            order = self.m.topological_order()
        except ValueError:
            # The cycle is already reported in the header; emit in declaration
            # order so the rest of the file is still readable.
            order = list(self.m.cells)

        if self.m.uses_multiscale():
            self.w(_rule("Scales"))
            self.w("#")
            self.w("# Inner scales are defined first. A phase property holding a")
            self.w("# `Homogenized(cell, scheme)` is the seam: the outer scheme")
            self.w("# resolves the inner scale when it reads that property, and")
            self.w("# memoizes it for the duration of the call.")
            self.blank()

        for c in order:
            self._builder(c)

    def _builder(self, c: Cell) -> None:
        args = ", ".join(c.params)
        self.w(f"function {c.builder}({args})")
        opts = c.rve_options or {}
        tail = ("; " + ", ".join(f"{k} = {v}" for k, v in sorted(opts.items()))) if opts else ""
        self.w(f"{IND}rve = RVE(:{c.matrix_name}{tail})")
        matrix = c.matrix()
        if matrix is not None:
            self._emit_phase(matrix, c)
        for ph in c.inclusions():
            self._emit_phase(ph, c)
        self.w(f"{IND}return rve")
        self.w("end")
        self.blank()

    def _emit_phase(self, ph: Phase, c: Cell) -> None:
        geom = self._geometry(ph)
        props = self._properties(ph)
        opts = []
        if ph.symmetrize and ph.symmetrize != "none":
            opts.append(
                "symmetrize = "
                + {"iso": "IsoSymmetrize()", "ti": "TISymmetrize()"}[ph.symmetrize]
            )

        if ph.is_matrix:
            # No amount: MFH derives the matrix fraction as 1 - Σ f_inclusions
            # and raises if it is set explicitly.
            tail = ("; " + ", ".join(opts)) if opts else ""
            self._call(f"add_matrix!(rve, {geom}, {props}{tail})")
            return

        opts.insert(0, f"{ph.amount_kind} = {self._amount(ph)}")
        self._call(f"add_phase!(rve, :{ph.name}, {geom}, {props}; " + ", ".join(opts) + ")")

    def _call(self, line: str) -> None:
        full = IND + line
        if len(full) <= 92:
            self.w(full)
            return
        head, _, tail = line.partition("(")
        body, _, _close = tail.rpartition(")")
        # The split must happen at the call's *own* keyword separator. A naive
        # search finds the one inside a nested constructor
        # (`SelfConsistent(; abstol = …)`) and produces unbalanced Julia, so
        # only depth-zero separators count.
        cut = _toplevel_kw_split(body)
        if cut is None:
            # nothing safe to break on — a long line beats a broken one
            self.w(full)
            return
        pos, kw = body[:cut], body[cut + 2:]
        self.w(f"{IND}{head}(")
        self.w(f"{IND}{IND}{pos};")
        self.w(f"{IND}{IND}{kw}")
        self.w(f"{IND})")

    def _amount(self, ph: Phase) -> str:
        return ph.amount if isinstance(ph.amount, str) else _fnum(ph.amount)

    def _geometry(self, ph: Phase) -> str:
        g = ph.geometry
        a = g.args
        v = lambda k, d=0.0: (a[k] if isinstance(a.get(k), str) else _fnum(a.get(k, d)))
        ang = ""
        if g.euler_angles:
            # Angles are physical quantities: a bare `0` next to `0.5` would
            # make the tuple `Tuple{Int, Float64}`.
            vals = ", ".join(
                x if isinstance(x, str) else _fnum(x) for x in g.euler_angles
            )
            ang = f"; euler_angles = ({vals},)" if len(g.euler_angles) == 1 else f"; euler_angles = ({vals})"

        k = g.kind
        if k == "spheroid":
            return f"Spheroid({v('omega', 1.0)}{ang})"
        if k == "ellipsoid":
            return f"Ellipsoid({v('a', 1.0)}, {v('b', 1.0)}, {v('c', 1.0)}{ang})"
        if k == "cylinder":
            return f"Cylinder({v('b', 1.0)}, {v('c', 1.0)}{ang})"
        if k == "penny_crack":
            return f"PennyCrack({v('a', 1.0)}{ang})"
        if k == "elliptic_crack":
            return f"EllipticCrack({v('a', 1.0)}, {v('b', 0.5)}{ang})"
        if k == "ribbon_crack":
            return f"RibbonCrack({v('b', 1.0)}{ang})"
        if k == "layered_sphere":
            return self._layered_sphere(g)
        if k == "layered_spheroid":
            return self._layered_spheroid(g)
        return "Spheroid(1.0)"

    def _layered_sphere(self, g) -> str:
        radii = ", ".join(
            (l["radius"] if isinstance(l.get("radius"), str) else _fnum(l.get("radius", 1.0)))
            for l in g.layers
        )
        moduli = ", ".join(self._prop_expr(Property.from_dict(l["property"])) for l in g.layers)
        radii_t = _tuple(radii.split(", ")) if radii else "()"
        mod_t = _tuple(moduli.split(", ")) if moduli else "()"
        ifaces = [l.get("interface") for l in g.layers]
        if any(i and i.get("kind", "PerfectInterface") != "PerfectInterface" for i in ifaces):
            parts = []
            for i in ifaces:
                i = i or {"kind": "PerfectInterface", "args": {}}
                kind = i.get("kind", "PerfectInterface")
                args = ", ".join(_num(x) for x in (i.get("args") or {}).values())
                parts.append(f"{kind}({args})" if args else f"{kind}()")
            it = f"({parts[0]},)" if len(parts) == 1 else "(" + ", ".join(parts) + ")"
            return f"LayeredSphere({radii_t}, {mod_t}; interfaces = {it})"
        return f"LayeredSphere({radii_t}, {mod_t})"

    def _layered_spheroid(self, g) -> str:
        """`layered_spheroid_from_fractions`, not the raw constructor.

        `LayeredSpheroid(axis_radii, disk_radii, …)` demands that every layer
        share one focal distance, and radii typed in layer by layer essentially
        never do — the old form here scaled one radius list by ω, which is not
        confocal and threw. The fraction constructor takes the outer aspect
        ratio and size and solves for the confocal inner radii itself, which is
        the only form a person can drive.
        """
        fractions = _tuple(
            (
                l["fraction"] if isinstance(l.get("fraction"), str)
                else _fnum(l.get("fraction", 1.0))
            )
            for l in g.layers
        )
        moduli = _tuple(
            self._prop_expr(Property.from_dict(l["property"])) for l in g.layers
        )
        omega = g.args.get("omega", 0.5)
        radius = g.args.get("radius", 1.0)
        ns = g.args.get("Nseries", 5)
        return (
            f"layered_spheroid_from_fractions({_fnum(omega)}, {_fnum(radius)}, "
            f"{fractions}, {moduli}; Nseries = {_num(ns)})"
        )

    def _properties(self, ph: Phase) -> str:
        entries = [f"{p.key} => {self._prop_expr(p)}" for p in ph.properties]
        if not entries:
            return "Dict{Symbol, Any}()"
        return "Dict(" + ", ".join(entries) + ")"

    def _prop_expr(self, p: Property) -> str:
        if p.source == "expr":
            return p.expr or "one(TensISO{3})"

        if p.source == "cell":
            # The multiscale seam.
            inner = self.m.cell(p.cell)
            if inner is None:
                return "#= missing inner cell =# one(TensISO{3})"
            call = f"{inner.builder}(" + ", ".join(inner.params) + ")"
            scheme = self._scheme(p.scheme or "MoriTanaka", p.scheme_options or {})
            return f"Homogenized({call}, {scheme})"

        if p.source == "visco" or p.visco:
            return self._visco_expr(p.visco or {})

        a = p.args
        num = lambda k, d=0.0: (a[k] if isinstance(a.get(k), str) else _fnum(a.get(k, d)))
        b = p.builder
        if b == "iso_stiffness":
            return f"iso_stiffness({num('k', 1.0)}, {num('mu', 1.0)})"
        if b == "iso_stiffness_E_nu":
            return f"iso_stiffness_E_nu({num('E', 1.0)}, {num('nu', 0.2)})"
        if b == "hoenig_stiffness":
            return (
                f"hoenig_stiffness({num('E1', 1.0)}, {num('h', 1.0)}, "
                f"{num('nu1', 0.2)}, {num('nu2', 0.2)}, {num('gamma', 1.0)})"
            )
        if b == "TensISO{3}":
            # One argument to TensISO{dim} is the 2nd-order (conductivity) form.
            return f"TensISO{{3}}({num('k', 1.0)})"
        return f"iso_stiffness({num('k', 1.0)}, {num('mu', 1.0)})"

    def _visco_expr(self, v: dict) -> str:
        kind = v.get("kind", "maxwell_iso")
        a = v.get("args", {})
        num = lambda k, d=1.0: _num(a.get(k, d))
        if kind == "maxwell_iso":
            return f"maxwell_iso({num('k', 10.0)}, {num('mu', 5.0)}, {num('tau')})"
        if kind == "kelvin_iso":
            return f"kelvin_iso({num('k', 10.0)}, {num('mu', 5.0)}, {num('tau')})"
        if kind == "heaviside":
            return (
                f"heaviside_law(iso_stiffness({num('k', 10.0)}, {num('mu', 5.0)}))"
            )
        expr = a.get("expr", "1.0")
        mode = v.get("mode", "creep")
        return f"ViscoLaw((t, t′) -> {expr}, :{mode})"

    # -- schemes ----------------------------------------------------------

    def _scheme(self, name: str, options: dict) -> str:
        opts = {k: v for k, v in (options or {}).items() if v is not None}
        if not opts:
            return f"{name}()"
        # Solver options attach to the scheme instance, not to `homogenize`.
        parts = ", ".join(f"{k} = {_num(v)}" for k, v in sorted(opts.items()))
        return f"{name}(; {parts})"

    # -- opaque -----------------------------------------------------------

    def _opaque(self) -> None:
        if not self.m.opaque:
            return
        self.w(_rule("Preserved from the original script"))
        self.w("#")
        self.w("# MFH Studio did not recognize the code below, so it is kept")
        self.w("# exactly as it was rather than rewritten.")
        self.blank()
        for blk in sorted(self.m.opaque, key=lambda b: b.order):
            if blk.note:
                self.w(f"# {blk.note}")
            for line in blk.source.splitlines():
                self.w(line)
            self.blank()

    # -- main -------------------------------------------------------------

    def _main(self) -> None:
        root = self.m.root()
        if root is None:
            return
        self.w(_rule("Result"))

        sw = self.m.sweep
        alv = self.m.alv

        if alv.enabled:
            self._alv_main(root)
            return

        if not sw.enabled:
            call = f"{root.builder}(" + ", ".join(root.params) + ")"
            scheme = self._scheme(sw.scheme, sw.scheme_options)
            self.w(f"C = homogenize({call}, {scheme}, {sw.property})")
            self._report("C", sw)
            return

        self._sweep_main(root, sw)

    def _report(self, var: str, sw) -> None:
        proj = {"iso": "best_fit_iso", "ti": "best_fit_ti", "ortho": "best_fit_ortho"}
        v = var
        if sw.projection in proj:
            self.w(f"{var}_fit = {proj[sw.projection]}({var})")
            v = f"{var}_fit"
        if sw.property == ":C":
            self.w(f'@printf "k = %.6f   μ = %.6f\\n" k_mu({v})[1] k_mu({v})[2]')
        else:
            self.w(f"@show {v}")

    def _sweep_main(self, root: Cell, sw) -> None:
        lens = self._lens_expr(sw.lens)
        scheme = self._scheme(sw.scheme, sw.scheme_options)
        proj = {"iso": "best_fit_iso", "ti": "best_fit_ti", "ortho": "best_fit_ortho"}

        self.w(f"const {sw.variable}s = range({_num(sw.start)}, {_num(sw.stop)}; length = {sw.length})")
        self.blank()
        self.w("#")
        self.w("# `set_param` returns a *new* cell, leaving the original intact,")
        self.w("# so the sweep is a pure map rather than a mutation.")
        self.w(f"const base_cell = {root.builder}(" + ", ".join(root.params) + ")")
        self.w(f"const scheme = {scheme}")
        self.blank()
        self.w(f"function evaluate({sw.variable})")
        self.w(f"{IND}cell = set_param(base_cell, {lens}, {sw.variable})")
        self.w(f"{IND}C = homogenize(cell, scheme, {sw.property})")
        if sw.projection in proj:
            self.w(f"{IND}C = {proj[sw.projection]}(C)")
        outs = []
        if sw.property == ":C":
            for o in sw.outputs:
                if o == "k":
                    outs.append("k_mu(C)[1]")
                elif o == "mu":
                    outs.append("k_mu(C)[2]")
                elif o == "E":
                    outs.append("E_nu(C)[1]")
                elif o == "nu":
                    outs.append("E_nu(C)[2]")
        if not outs:
            outs = ["tr(Array(C)) / 3"]
        self.w(f"{IND}return ({', '.join(outs)})")
        self.w("end")
        self.blank()
        self.w(f"const results = [evaluate({sw.variable}) for {sw.variable} in {sw.variable}s]")
        for i, o in enumerate(sw.outputs or ["value"]):
            self.w(f"const {o}s = [r[{i + 1}] for r in results]")
        self.blank()

        # A machine-readable handle for the interface to plot without parsing
        # the script's stdout.
        names = sw.outputs or ["value"]
        pairs = ", ".join(f'"{o}" => {o}s' for o in names)
        self.w("# Published for MFH Studio; harmless when the script runs alone.")
        self.w(
            f'const MFHSTUDIO_RESULTS = Dict("x" => collect({sw.variable}s), '
            f'"xlabel" => "{sw.variable}", {pairs})'
        )
        self.blank()

        self.w('@printf "%10s' + '%14s' * len(names) + '\\n" "' + sw.variable + '" ' +
               " ".join(f'"{o}"' for o in names))
        self.w(f"for (i, {sw.variable}) in enumerate({sw.variable}s)")
        self.w(
            f'{IND}@printf "%10.4f' + "%14.6f" * len(names) + '\\n" '
            + sw.variable + " " + " ".join(f"{o}s[i]" for o in names)
        )
        self.w("end")

        if sw.plot:
            self.blank()
            self.w("p = plot(; xlabel = \"" + sw.variable + "\", ylabel = \"effective property\", legend = :best)")
            for o in names:
                self.w(f'plot!(p, {sw.variable}s, {o}s; label = "{o}", lw = 2)')
            self.w("display(p)")

    def _alv_main(self, root: Cell) -> None:
        alv = self.m.alv
        cell = self.m.cell(alv.cell) or root
        call = f"{cell.builder}(" + ", ".join(cell.params) + ")"
        scheme = self._scheme(alv.scheme, {})
        if alv.log_time:
            self.w(
                f"const times = 10 .^ range({_num(alv.t_start)}, {_num(alv.t_stop)}; "
                f"length = {alv.length})"
            )
        else:
            self.w(
                f"const times = range({_num(alv.t_start)}, {_num(alv.t_stop)}; "
                f"length = {alv.length})"
            )
        self.blank()
        self.w(f"const R = homogenize_alv({call}, {scheme}, {alv.property}; times = times)")
        self.w("@show size(R)")

    # -- the embedded model ----------------------------------------------

    def _embedded_model(self) -> None:
        if not self.embed:
            return
        self.blank()
        self.w("#=" + f" {MODEL_TAG}")
        self.w("The studio reads this block to reopen the model exactly as it was.")
        self.w("Deleting it costs nothing but a best-effort re-reading of the code.")
        self.w(json.dumps(self.m.to_dict(), indent=1, sort_keys=True))
        self.w("=#")

    # -- lenses -----------------------------------------------------------

    def _lens_expr(self, lens) -> str:
        k = lens.kind
        if k == "amount":
            return f"amount(:{lens.phase})"
        if k == "property":
            return f"property(:{lens.phase}, {lens.property}, {lens.index})"
        if k == "geometry":
            return f"geometry(:{lens.phase}, :{lens.field_name}, {lens.index})"
        if k == "shape_param":
            return f"shape_param(:{lens.field_name}, {lens.index})"
        if k == "nested":
            from .model import Lens as _L

            inner = _L.from_dict(lens.inner or {})
            return f"nested(:{lens.member}, {lens.property}, {self._lens_expr(inner)})"
        return f"amount(:{lens.phase})"


def _toplevel_kw_split(body: str) -> Optional[int]:
    """Index of the `; ` that separates positional from keyword arguments.

    Only a separator at bracket depth zero belongs to this call; any deeper one
    belongs to a nested constructor.
    """
    depth = 0
    for i, ch in enumerate(body):
        if ch in "([{":
            depth += 1
        elif ch in ")]}":
            depth -= 1
        elif ch == ";" and depth == 0:
            return i
    return None


def generate(model: Model, embed_model: bool = True) -> str:
    return CodeGen(model, embed_model).generate()


def render_cell(model: Model, cell: Cell) -> str:
    """Just one builder, as it would appear in the script.

    Used by the reader to check that a cell it *thinks* it understood would be
    written back unchanged. A construct the studio can parse but not reproduce
    is one it should leave alone.
    """
    g = CodeGen(model, embed_model=False)
    g._builder(cell)
    return "\n".join(g.out).rstrip()


def extract_embedded(source: str) -> Optional[dict]:
    """The model a studio-written script carries, if any."""
    i = source.find("#=" + f" {MODEL_TAG}")
    if i < 0:
        return None
    j = source.find("\n=#", i)
    if j < 0:
        return None
    body = source[i:j]
    brace = body.find("{")
    if brace < 0:
        return None
    try:
        return json.loads(body[brace:])
    except json.JSONDecodeError:
        return None
