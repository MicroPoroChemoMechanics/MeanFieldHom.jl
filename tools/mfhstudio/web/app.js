/* MFH Studio — the browser side.
 *
 * The model lives on the server; this file edits a local copy, POSTs it back
 * and re-renders from the answer. That keeps one source of truth and means the
 * Julia script shown on the right is always the one that would be saved.
 */

const $ = (s) => document.querySelector(s);
const el = (tag, attrs = {}, ...kids) => {
  const n = document.createElement(tag);
  for (const [k, v] of Object.entries(attrs)) {
    if (k === "class") n.className = v;
    else if (k === "html") n.innerHTML = v;
    else if (k.startsWith("on")) n.addEventListener(k.slice(2), v);
    else if (v === true) n.setAttribute(k, "");
    else if (v !== false && v != null) n.setAttribute(k, v);
  }
  for (const kid of kids.flat()) {
    if (kid == null) continue;
    n.append(kid.nodeType ? kid : document.createTextNode(String(kid)));
  }
  return n;
};

const S = {
  model: null,
  catalog: null,
  cellId: null,
  phaseIdx: 0,
  problems: [],
  keptReport: null,
};

/* ── server ─────────────────────────────────────────────────────── */

async function api(path, body) {
  const opt = body === undefined
    ? {}
    : { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(body) };
  const r = await fetch(path, opt);
  const j = await r.json().catch(() => ({ error: `${r.status} ${r.statusText}` }));
  if (j && j.error) throw new Error(j.error);
  return j;
}

function toast(msg, bad = false) {
  const t = el("div", { class: bad ? "bad" : "" }, msg);
  $("#toast").append(t);
  setTimeout(() => t.remove(), bad ? 9000 : 3500);
}

/** Push the edited model and re-render from the server's answer. */
async function push() {
  try {
    const st = await api("/api/model", { model: S.model });
    apply(st);
  } catch (e) {
    toast(e.message, true);
  }
}

function apply(st) {
  S.model = st.model;
  S.problems = st.problems || [];
  if (!S.cellId || !S.model.cells.some((c) => c.id === S.cellId)) {
    S.cellId = S.model.cells.length ? S.model.cells[0].id : null;
  }
  $("#code").textContent = st.script;
  $("#problems").textContent = S.problems.length
    ? `${S.problems.length} problem(s)`
    : "";
  $("#problems").className = S.problems.length ? "problem" : "muted";
  if (st.path) $("#path").value = st.path;
  render();
  draw3d();
}

/* ── helpers over the model ─────────────────────────────────────── */

const cell = () => S.model.cells.find((c) => c.id === S.cellId) || null;
const phases = () => (cell() ? cell().phases : []);
const phase = () => phases()[S.phaseIdx] || null;

const geomForm = (kind) =>
  (S.catalog.geometries || []).find((g) => g.kind === kind) || null;
const propForm = (name) =>
  (S.catalog.properties || []).find((p) => p.name === name) || null;

/** Render a value as a Julia *float* literal.
 *
 * Geometry constructors take an NTuple of one element type, so a stray
 * integer makes `LayeredSphere((0.6, 1), …)` a `Tuple{Float64, Int64}` and
 * the call fails to dispatch. JavaScript has no int/float distinction, so
 * the decimal point has to be put back explicitly.
 */
function jnum(v) {
  if (typeof v === "string") return v.trim() === "" ? "0.0" : v;
  if (!isFinite(v)) return String(v);
  return Number.isInteger(v) ? v.toFixed(1) : String(v);
}

/** The Julia expression for a phase's geometry, for the 3-D preview. */
function geomExpr(ph) {
  if (!ph) return "";
  const g = ph.geometry, a = g.args || {};
  const n = (k, d) => {
    const v = a[k];
    return jnum(v === undefined || v === "" ? d : v);
  };
  switch (g.kind) {
    case "spheroid": return `Spheroid(${n("omega", 1)})`;
    case "ellipsoid": return `Ellipsoid(${n("a", 1)}, ${n("b", 1)}, ${n("c", 1)})`;
    case "cylinder": return `Cylinder(${n("b", 1)}, ${n("c", 1)})`;
    case "penny_crack": return `PennyCrack(${n("a", 1)})`;
    case "elliptic_crack": return `EllipticCrack(${n("a", 1)}, ${n("b", 0.5)})`;
    case "ribbon_crack": return `RibbonCrack(${n("b", 1)})`;
    case "layered_sphere": {
      const r = (g.layers || []).map((l) => jnum(l.radius));
      if (!r.length) return "";
      const mods = (g.layers || []).map(() => "iso_stiffness(1.0, 1.0)");
      const t = (xs) => (xs.length === 1 ? `(${xs[0]},)` : `(${xs.join(", ")})`);
      return `LayeredSphere(${t(r)}, ${t(mods)})`;
    }
    default: return "";
  }
}

/* ── rendering ──────────────────────────────────────────────────── */

function render() {
  renderScales();
  renderSweep();
  renderAlv();
  renderParams();
  renderKept();
}

function renderScales() {
  const sel = $("#cell-select");
  sel.replaceChildren(
    ...S.model.cells.map((c) =>
      el("option", { value: c.id, selected: c.id === S.cellId }, c.name)
    )
  );
  const c = cell();
  $("#cell-name").value = c ? c.name : "";
  $("#cell-matrix").value = c ? c.matrix_name : "";

  // The graph only earns its place once there is more than one scale: with a
  // single RVE there is nothing to connect and it would be noise.
  const multi = S.model.cells.length > 1;
  $("#graph").hidden = !multi;
  $("#graph-head").hidden = !multi;
  if (multi) drawGraph();

  $("#phases").replaceChildren(...phases().map(phaseCard));
}

function phaseCard(ph, i) {
  const seam = (ph.properties || []).some((p) => p.source === "cell");
  const head = el(
    "header", {},
    el("b", {}, ph.name || "(unnamed)"),
    ph.is_matrix ? el("span", { class: "tag matrix" }, "matrix") : null,
    seam ? el("span", { class: "tag seam" }, "nested scale") : null,
    el("button", {
      class: "small",
      title: "Remove this phase",
      onclick: (e) => { e.stopPropagation(); cell().phases.splice(i, 1); push(); },
    }, "−")
  );

  const body = el("div", {});
  if (i === S.phaseIdx) {
    body.append(
      field("Name", input(ph.name, (v) => { ph.name = v; push(); })),
      field("Role", select(
        [["inclusion", "inclusion"], ["matrix", "matrix"]],
        ph.is_matrix ? "matrix" : "inclusion",
        (v) => {
          ph.is_matrix = v === "matrix";
          if (ph.is_matrix) {
            // MFH derives the matrix amount as 1 − Σ f and refuses to be told
            // otherwise, so the field simply goes away.
            for (const o of cell().phases) if (o !== ph) o.is_matrix = false;
            cell().matrix_name = ph.name;
          }
          push();
        }
      )),
      ph.is_matrix
        ? el("div", { class: "note" },
            "The matrix amount is derived as 1 − Σ f of the inclusions; MeanFieldHom raises if it is set.")
        : el("div", { class: "grid2" },
            field("Amount", select(
              [["fraction", "volume fraction"], ["density", "crack density"]],
              ph.amount_kind, (v) => { ph.amount_kind = v; push(); }
            )),
            field("Value", input(ph.amount, (v) => {
              ph.amount = isFinite(+v) && v.trim() !== "" ? +v : v;
              push();
            }))
          ),
      geometryEditor(ph),
      field("Orientation average", select(
        (S.catalog.symmetrize || []).map((s) => [s.name, s.label]),
        ph.symmetrize, (v) => { ph.symmetrize = v; push(); }
      )),
      el("h3", {}, "Properties",
        el("button", {
          class: "small",
          onclick: () => {
            ph.properties.push({
              key: ":C", source: "builder", builder: "iso_stiffness",
              form: "iso_kmu", args: { k: 10, mu: 5 }, scheme_options: {},
            });
            push();
          },
        }, "+")),
      ...(ph.properties || []).map((pr, j) => propertyEditor(ph, pr, j))
    );
  }

  return el("div", {
    class: "card" + (i === S.phaseIdx ? " selected" : ""),
    onclick: () => { S.phaseIdx = i; render(); draw3d(); },
  }, head, body);
}

function geometryEditor(ph) {
  const g = ph.geometry;
  const form = geomForm(g.kind);
  const box = el("div", {},
    field("Shape", select(
      (S.catalog.geometries || []).map((x) => [x.kind, x.name]),
      g.kind,
      (v) => {
        g.kind = v;
        const f = geomForm(v);
        g.args = {};
        for (const fl of (f && f.fields) || []) g.args[fl.name] = fl.default;
        if (f && f.layered && !(g.layers || []).length) {
          g.layers = [
            { radius: 0.6, property: { key: ":C", builder: "iso_stiffness", form: "void", args: { k: 1e-6, mu: 1e-6 }, source: "builder", scheme_options: {} }, interface: { kind: "PerfectInterface", args: {} } },
            { radius: 1.0, property: { key: ":C", builder: "iso_stiffness", form: "iso_kmu", args: { k: 30, mu: 12 }, source: "builder", scheme_options: {} }, interface: { kind: "PerfectInterface", args: {} } },
          ];
        }
        push();
      }
    ))
  );
  if (form && form.doc) box.append(el("div", { class: "note" }, form.doc));
  const fields = (form && form.fields) || [];
  if (fields.length) {
    box.append(el("div", { class: fields.length > 2 ? "grid3" : "grid2" },
      ...fields.map((f) =>
        field(f.label, input(g.args[f.name] ?? f.default, (v) => {
          g.args[f.name] = isFinite(+v) && v.trim() !== "" ? +v : v;
          push();
        }))
      )));
  }
  if (form && form.layered) box.append(layersEditor(g));
  return box;
}

function layersEditor(g) {
  const box = el("div", {},
    el("h3", {}, "Layers",
      el("button", {
        class: "small",
        onclick: () => {
          g.layers = g.layers || [];
          const last = g.layers.length ? g.layers[g.layers.length - 1].radius : 0.5;
          g.layers.push({
            radius: (+last || 0.5) + 0.5,
            property: { key: ":C", source: "builder", builder: "iso_stiffness", form: "iso_kmu", args: { k: 10, mu: 5 }, scheme_options: {} },
            interface: { kind: "PerfectInterface", args: {} },
          });
          push();
        },
      }, "+")),
    el("div", { class: "note" }, "Ascending radii, r = 0 implicit at the center; layer 1 is the core.")
  );
  (g.layers || []).forEach((l, i) => {
    box.append(el("div", { class: "card" },
      el("header", {}, el("b", {}, `layer ${i + 1}`),
        el("button", { class: "small", onclick: () => { g.layers.splice(i, 1); push(); } }, "−")),
      el("div", { class: "grid3" },
        field("outer r", input(l.radius, (v) => { l.radius = +v || v; push(); })),
        field("k", input(l.property.args.k, (v) => { l.property.args.k = +v || v; push(); })),
        field("μ", input(l.property.args.mu, (v) => { l.property.args.mu = +v || v; push(); }))
      ),
      field("Interface with the next layer", select(
        (S.catalog.interfaces || []).map((x) => [x.name, x.label]),
        (l.interface && l.interface.kind) || "PerfectInterface",
        (v) => {
          const f = (S.catalog.interfaces || []).find((x) => x.name === v);
          const args = {};
          for (const fl of (f && f.fields) || []) args[fl.name] = fl.default;
          l.interface = { kind: v, args };
          push();
        }
      ))
    ));
  });
  return box;
}

function propertyEditor(ph, pr, j) {
  const others = S.model.cells.filter((c) => c.id !== S.cellId);
  const sourceOptions = [
    ["builder", "from moduli"],
    ["expr", "Julia expression"],
  ];
  if (others.length) sourceOptions.push(["cell", "another scale (nested)"]);

  const box = el("div", { class: "card" },
    el("header", {},
      el("b", {}, pr.key),
      el("button", { class: "small", onclick: () => { ph.properties.splice(j, 1); push(); } }, "−")),
    el("div", { class: "grid2" },
      field("Key", input(pr.key, (v) => { pr.key = v.startsWith(":") ? v : ":" + v; push(); })),
      field("Source", select(sourceOptions, pr.source, (v) => { pr.source = v; push(); }))
    )
  );

  if (pr.source === "builder") {
    box.append(field("Parametrization", select(
      (S.catalog.properties || []).map((p) => [p.name, p.label]),
      pr.form || "iso_kmu",
      (v) => {
        const f = propForm(v);
        pr.form = v;
        pr.builder = (f && f.builder) || "iso_stiffness";
        pr.args = {};
        for (const fl of (f && f.fields) || []) pr.args[fl.name] = fl.default;
        push();
      }
    )));
    const f = propForm(pr.form || "iso_kmu");
    const fields = (f && f.fields) || [];
    box.append(el("div", { class: fields.length > 2 ? "grid3" : "grid2" },
      ...fields.map((fl) =>
        field(fl.label, input(pr.args[fl.name] ?? fl.default, (v) => {
          pr.args[fl.name] = isFinite(+v) && v.trim() !== "" ? +v : v;
          push();
        }))
      )));
  } else if (pr.source === "expr") {
    box.append(field("Julia", input(pr.expr || "", (v) => { pr.expr = v; push(); })));
  } else if (pr.source === "cell") {
    box.append(
      el("div", { class: "note" },
        "The seam: this property is the effective value of another scale, "
        + "resolved by the outer scheme when it reads the key."),
      el("div", { class: "grid2" },
        field("Inner scale", select(
          others.map((c) => [c.id, c.name]),
          pr.cell || (others[0] && others[0].id),
          (v) => { pr.cell = v; push(); }
        )),
        field("Homogenized with", select(
          (S.catalog.schemes || []).map((s) => [s.name, s.name]),
          pr.scheme || "MoriTanaka",
          (v) => { pr.scheme = v; push(); }
        ))
      ),
      schemeOptions(pr.scheme || "MoriTanaka", pr.scheme_options || (pr.scheme_options = {}))
    );
  }
  return box;
}

/** Editable solver options for a scheme, straight from the catalog. */
function schemeOptions(name, target) {
  const s = (S.catalog.schemes || []).find((x) => x.name === name);
  const opts = ((s && s.options) || []).filter((o) => o.editable);
  if (!opts.length) return el("div", { class: "note" }, "No solver options.");
  return el("div", { class: "grid2" },
    ...opts.map((o) =>
      field(o.name, typeof o.default === "boolean"
        ? checkbox(target[o.name] ?? o.default, (v) => { target[o.name] = v; push(); })
        : input(target[o.name] ?? o.default, (v) => {
            target[o.name] = v.trim() === "" ? null : +v;
            push();
          }))
    ));
}

/* ── sweep ──────────────────────────────────────────────────────── */

function renderSweep() {
  const sw = S.model.sweep;
  const t = $("#tab-sweep");
  const inner = sw.lens.inner || (sw.lens.inner = { kind: "amount", phase: "", property: ":C", field_name: "semi_axes", index: 1, member: "", inner: null });

  t.replaceChildren(
    field("", checkboxLabel("Sweep a parameter", sw.enabled, (v) => { sw.enabled = v; push(); })),
    el("div", { class: "grid3" },
      field("Variable", input(sw.variable, (v) => { sw.variable = v || "x"; push(); })),
      field("From", input(sw.start, (v) => { sw.start = +v; push(); })),
      field("To", input(sw.stop, (v) => { sw.stop = +v; push(); }))
    ),
    el("div", { class: "grid2" },
      field("Points", input(sw.length, (v) => { sw.length = Math.max(2, +v | 0); push(); })),
      field("Scheme", select(
        (S.catalog.schemes || []).map((s) => [s.name, s.name]),
        sw.scheme, (v) => { sw.scheme = v; push(); }
      ))
    ),
    schemeOptions(sw.scheme, sw.scheme_options || (sw.scheme_options = {})),
    el("h3", {}, "What varies"),
    field("Lens", select(
      (S.catalog.lenses || []).map((l) => [l.name, l.label]),
      sw.lens.kind, (v) => { sw.lens.kind = v; push(); }
    )),
    lensDoc(sw.lens.kind),
    ...lensFields(sw.lens, inner),
    el("h3", {}, "Output"),
    el("div", { class: "grid2" },
      field("Property", input(sw.property, (v) => { sw.property = v.startsWith(":") ? v : ":" + v; push(); })),
      field("Report as", select(
        (S.catalog.projections || []).map((p) => [p.name, p.label]),
        sw.projection, (v) => { sw.projection = v; push(); }
      ))
    ),
    field("Plot", checkboxLabel("draw a figure", sw.plot, (v) => { sw.plot = v; push(); }))
  );
}

function lensDoc(kind) {
  const l = (S.catalog.lenses || []).find((x) => x.name === kind);
  return l && l.doc ? el("div", { class: "note" }, l.doc) : el("span");
}

function lensFields(lens, inner) {
  const names = phases().map((p) => [p.name, p.name]);
  const out = [];
  if (lens.kind === "nested") {
    out.push(
      el("div", { class: "grid2" },
        field("Through phase", select(names, lens.member || (names[0] && names[0][0]), (v) => { lens.member = v; push(); })),
        field("Key", input(lens.property, (v) => { lens.property = v.startsWith(":") ? v : ":" + v; push(); }))
      ),
      el("h3", {}, "Inside that scale"),
      field("Lens", select(
        (S.catalog.lenses || []).filter((l) => l.name !== "nested").map((l) => [l.name, l.label]),
        inner.kind, (v) => { inner.kind = v; push(); }
      )),
      ...lensFields(inner, {})
    );
    return out;
  }
  if (lens.kind === "amount") {
    out.push(field("Phase", select(names, lens.phase || (names[0] && names[0][0]), (v) => { lens.phase = v; push(); })));
  } else if (lens.kind === "property") {
    out.push(el("div", { class: "grid3" },
      field("Phase", select(names, lens.phase, (v) => { lens.phase = v; push(); })),
      field("Key", input(lens.property, (v) => { lens.property = v; push(); })),
      field("Index", input(lens.index, (v) => { lens.index = +v | 0 || 1; push(); }))
    ));
  } else if (lens.kind === "geometry") {
    out.push(el("div", { class: "grid3" },
      field("Phase", select(names, lens.phase, (v) => { lens.phase = v; push(); })),
      field("Field", input(lens.field_name, (v) => { lens.field_name = v; push(); })),
      field("Index", input(lens.index, (v) => { lens.index = +v | 0 || 1; push(); }))
    ));
  } else if (lens.kind === "shape_param") {
    out.push(el("div", { class: "grid2" },
      field("Field", input(lens.field_name, (v) => { lens.field_name = v; push(); })),
      field("Index", input(lens.index, (v) => { lens.index = +v | 0 || 1; push(); }))
    ));
  }
  return out;
}

/* ── viscoelasticity ────────────────────────────────────────────── */

function renderAlv() {
  const a = S.model.alv;
  const blocked = S.model.cells.some((c) =>
    c.phases.some((p) => (p.properties || []).some((x) => x.source === "cell")));

  $("#tab-alv").replaceChildren(
    blocked
      ? el("div", { class: "note problem" },
          "Ageing viscoelasticity cannot be combined with a nested scale: "
          + "MeanFieldHom cannot re-express a homogenized inner result as a ViscoLaw. "
          + "Remove the seam first.")
      : el("span"),
    field("", checkboxLabel("Ageing linear viscoelastic run", a.enabled, (v) => { a.enabled = v; push(); })),
    el("div", { class: "grid3" },
      field("t from", input(a.t_start, (v) => { a.t_start = +v; push(); })),
      field("t to", input(a.t_stop, (v) => { a.t_stop = +v; push(); })),
      field("Steps", input(a.length, (v) => { a.length = Math.max(2, +v | 0); push(); }))
    ),
    field("", checkboxLabel("logarithmic time", a.log_time, (v) => { a.log_time = v; push(); })),
    field("Scheme", select(
      (S.catalog.schemes || []).map((s) => [s.name, s.name]),
      a.scheme, (v) => { a.scheme = v; push(); }
    )),
    el("div", { class: "note" },
      "Give a phase a viscoelastic property by choosing “Julia expression” and "
      + "writing e.g. maxwell_iso(10.0, 5.0, 1.0) or ViscoLaw((t, t′) -> …, :creep).")
  );
}

/* ── constants and kept blocks ──────────────────────────────────── */

function renderParams() {
  const t = $("#tab-params");
  t.replaceChildren(
    el("h3", {}, "Constants",
      el("button", {
        class: "small",
        onclick: () => { S.model.params.push({ name: "x", value: "1.0", comment: "", origin: null, edited: true }); push(); },
      }, "+")),
    ...S.model.params.map((p, i) =>
      el("div", { class: "card" },
        el("header", {}, el("b", {}, p.name),
          p.origin && !p.edited ? el("span", { class: "tag" }, "as written") : null,
          el("button", { class: "small", onclick: () => { S.model.params.splice(i, 1); push(); } }, "−")),
        el("div", { class: "grid2" },
          field("Name", input(p.name, (v) => { p.name = v; p.edited = true; push(); })),
          field("Value", input(p.value, (v) => { p.value = v; p.edited = true; push(); }))
        ))
    ),
    S.model.params.length ? el("div", { class: "note" },
      "A constant read from a file keeps its original text until you edit it, "
      + "so multi-line layout is not collapsed.") : el("span")
  );
}

function renderKept() {
  const t = $("#tab-kept");
  const kept = S.model.opaque || [];
  const rep = S.keptReport;
  t.replaceChildren(
    el("div", { class: "note" },
      "Code MFH Studio did not recognize. It is written back unchanged — "
      + "editing the model never rewrites it."),
    rep && rep.exact
      ? el("div", { class: "note" }, "This file was written by the studio and reopened exactly.")
      : el("span"),
    ...kept.map((o) =>
      el("div", {},
        o.note ? el("div", { class: "muted" }, o.note) : null,
        el("div", { class: "kept" }, o.source))
    ),
    kept.length ? null : el("div", { class: "muted" }, "Nothing kept.")
  );
}

/* ── small widgets ──────────────────────────────────────────────── */

function field(label, node) {
  return el("div", { class: "field" }, label ? el("label", {}, label) : null, node);
}
function input(value, on) {
  const n = el("input", { type: "text", spellcheck: "false" });
  n.value = value ?? "";
  n.addEventListener("change", () => on(n.value));
  return n;
}
function select(pairs, value, on) {
  const n = el("select", {}, ...pairs.map(([v, t]) =>
    el("option", { value: v, selected: String(v) === String(value) }, t)));
  n.addEventListener("change", () => on(n.value));
  return n;
}
function checkbox(value, on) {
  const n = el("input", { type: "checkbox" });
  n.checked = !!value;
  n.addEventListener("change", () => on(n.checked));
  return n;
}
function checkboxLabel(text, value, on) {
  return el("label", { class: "inline" }, checkbox(value, on), " " + text);
}

/* ── 3-D ────────────────────────────────────────────────────────── */

let lastExpr = null;
async function draw3d() {
  const ph = phase();
  const expr = geomExpr(ph);
  $("#shape-label").textContent = expr || "";
  if (!expr) { Plotly.purge("view3d"); lastExpr = null; return; }
  const key = expr + "|" + $("#cutaway").checked;
  if (key === lastExpr) return;
  lastExpr = key;
  try {
    const sc = await api("/api/traces", { expr, cutaway: $("#cutaway").checked });
    Plotly.react("view3d", sc.data, sc.layout, { displayModeBar: false, responsive: true });
  } catch (e) {
    // A shape the sidecar cannot build is a modeling error worth showing.
    $("#shape-label").textContent = expr + " — " + e.message.split("\n")[0];
  }
}

/* ── run ────────────────────────────────────────────────────────── */

async function run() {
  const btn = $("#run");
  btn.disabled = true;
  btn.textContent = "Running…";
  $("#stdout").textContent = "";
  try {
    const r = await api("/api/run", {});
    $("#stdout").textContent = (r.stdout || "") + (r.error ? "\n" + r.error : "");
    if (r.results) plotResults(r.results);
    if (r.ok) toast("Ran cleanly.");
    else toast(r.timeout ? "Timed out." : "The script raised — see the output.", true);
  } catch (e) {
    toast(e.message, true);
  } finally {
    btn.disabled = false;
    btn.textContent = "Run";
  }
}

function plotResults(res) {
  const x = res.x || [];
  const traces = Object.entries(res)
    .filter(([k]) => k !== "x" && k !== "xlabel")
    .map(([k, y]) => ({ x, y, name: k, type: "scatter", mode: "lines+markers" }));
  if (!traces.length) return;
  Plotly.react("plot", traces, {
    margin: { l: 48, r: 10, t: 10, b: 38 },
    xaxis: { title: res.xlabel || "x" },
    paper_bgcolor: "rgba(0,0,0,0)",
    plot_bgcolor: "rgba(0,0,0,0)",
    legend: { orientation: "h" },
  }, { displayModeBar: false, responsive: true });
}

/* ── files ──────────────────────────────────────────────────────── */

async function openFile() {
  const path = $("#path").value.trim();
  if (!path) return toast("Give a path first.", true);
  try {
    const st = await api("/api/open", { path });
    S.keptReport = st.read_report || null;
    apply(st);
    const r = st.read_report || {};
    toast(r.exact
      ? "Reopened from the model embedded in the file."
      : `Read: ${r.recognized || 0} construct(s) understood, ${r.opaque || 0} kept as written.`);
  } catch (e) {
    toast(e.message, true);
  }
}

async function saveFile() {
  const path = $("#path").value.trim();
  if (!path) return toast("Give a path first.", true);
  try {
    const r = await api("/api/save", { path });
    toast(`Saved ${r.bytes} bytes to ${r.path}`);
  } catch (e) {
    toast(e.message, true);
  }
}

/* ── sidecar status ─────────────────────────────────────────────── */

async function pollSidecar() {
  try {
    const s = await api("/api/sidecar");
    const b = $("#sidecar");
    if (s.ready) { b.textContent = "Julia ready"; b.className = "badge ok"; }
    else if (s.running) { b.textContent = "loading MeanFieldHom…"; b.className = "badge"; }
    else if (s.error) { b.textContent = "Julia unavailable"; b.className = "badge bad"; b.title = s.error; }
    else { b.textContent = "starting…"; b.className = "badge"; }
    if (s.ready && !S.catalog) await boot();
  } catch { /* the server will answer eventually */ }
}

/* ── boot ───────────────────────────────────────────────────────── */

async function boot() {
  try {
    S.catalog = await api("/api/catalog");
  } catch (e) {
    return; // still loading; the poll will retry
  }
  const st = await api("/api/state");
  apply(st);
}

function wire() {
  document.querySelectorAll(".tabs button").forEach((b) =>
    b.addEventListener("click", () => {
      document.querySelectorAll(".tabs button").forEach((x) => x.classList.remove("on"));
      document.querySelectorAll(".tab").forEach((x) => x.classList.remove("on"));
      b.classList.add("on");
      $("#tab-" + b.dataset.tab).classList.add("on");
    }));

  $("#cell-select").addEventListener("change", (e) => {
    S.cellId = e.target.value; S.phaseIdx = 0; render(); draw3d();
  });
  $("#cell-add").addEventListener("click", () => {
    const n = S.model.cells.length + 1;
    // Stagger new boxes so they never land on top of each other; the user can
    // drag them wherever afterwards and the position is remembered.
    S.model.cells.push({
      id: "", name: `scale${n}`, matrix_name: "MATRIX", params: [],
      builder_name: null, rve_options: {},
      ui: { x: 40 + ((n - 1) % 3) * 210, y: 30 + Math.floor((n - 1) / 3) * 150 },
      phases: [{
        name: "MATRIX", is_matrix: true, amount_kind: "fraction", amount: 0, symmetrize: "none",
        geometry: { kind: "spheroid", args: { omega: 1.0 }, euler_angles: [], layers: [] },
        properties: [{ key: ":C", source: "builder", builder: "iso_stiffness", form: "iso_kmu", args: { k: 10, mu: 5 }, scheme_options: {} }],
      }],
    });
    push();
  });
  $("#cell-del").addEventListener("click", () => {
    const i = S.model.cells.findIndex((c) => c.id === S.cellId);
    if (i >= 0 && S.model.cells.length > 1) { S.model.cells.splice(i, 1); S.cellId = null; push(); }
  });
  $("#cell-name").addEventListener("change", (e) => { cell().name = e.target.value; push(); });
  $("#cell-matrix").addEventListener("change", (e) => { cell().matrix_name = e.target.value; push(); });
  $("#phase-add").addEventListener("click", () => {
    cell().phases.push({
      name: "PHASE" + cell().phases.length, is_matrix: false,
      amount_kind: "fraction", amount: 0.1, symmetrize: "none",
      geometry: { kind: "spheroid", args: { omega: 1.0 }, euler_angles: [], layers: [] },
      properties: [{ key: ":C", source: "builder", builder: "iso_stiffness", form: "iso_kmu", args: { k: 10, mu: 5 }, scheme_options: {} }],
    });
    S.phaseIdx = cell().phases.length - 1;
    push();
  });

  $("#cutaway").addEventListener("change", () => { lastExpr = null; draw3d(); });
  $("#run").addEventListener("click", run);
  $("#open").addEventListener("click", openFile);
  $("#save").addEventListener("click", saveFile);
}

wire();
pollSidecar();
setInterval(pollSidecar, 2000);
