# MFH Studio

A local web interface for building MeanFieldHom scripts: draw the
microstructure, run the model, and read existing scripts back without damaging
them.

```bash
cd tools/mfhstudio
python3 -m mfhstudio                 # start and open a browser
python3 -m mfhstudio --port 9000     # pick the port
python3 -m mfhstudio --no-browser    # stay in the terminal
python3 -m mfhstudio --check         # verify the Julia side and exit
```

Requires Python 3.10+ (standard library only — no `pip install`) and a `julia`
on `PATH` able to load MeanFieldHom. Set `JULIA` to point at a specific
executable.

The user-facing guide, with screenshots, is
`docs/src/manual/mfhstudio.md`. This file covers how the thing is put together.

## Shape

```
browser (HTML/JS, Plotly, draggable scale graph)
    │  REST
    ▼
Python  mfhstudio/            model state, Julia generation, script read-back
    │  JSON lines over stdio
    ▼
Julia   sidecar               3-D traces, Meta.parse, execution
```

The model lives in Python. That is what makes the sidecar disposable: it can be
restarted after a wedge without losing work.

| file | role |
| :--- | :--- |
| `mfhstudio/model.py` | the authoring model — a *graph of cells*, not one RVE |
| `mfhstudio/codegen.py` | model → Julia |
| `mfhstudio/readback.py` | script → model, with verbatim preservation |
| `mfhstudio/juliabridge.py` | the sidecar process: start, call, restart |
| `mfhstudio/server.py` | HTTP endpoints and session state |
| `julia/sidecar.jl` | the JSON-lines loop |
| `julia/introspect.jl` | the feature catalogue, read from the live package |
| `julia/geometry.jl` | 3-D traces |
| `julia/parse_script.jl` | `Meta.parse` → nodes that tile the file exactly |
| `web/graph.js` | the draggable scale graph |

## Three design decisions worth knowing

**The catalogue is introspected, never hard-coded.** Schemes come from
`subtypes(HomogenizationScheme)`, and each scheme's solver options are read
from the constant it declares for the purpose (`_SC_SOLVER_KWARGS`,
`_DIFF_RESERVED_OPTIONS`). Probing the constructor would not do: those schemes
take a `kwargs...` bag that accepts anything, so `SelfConsistent(; nsteps = 3)`
succeeds while `nsteps` is meaningless there. Reading the declared keys is what
keeps the interface exactly in step with the schemes.

**3-D reuses the parametrizations, not the trace builders.**
`scripts/common/docviz.jl` — the code that draws the documentation figures —
returns JavaScript object literals from its `*_trace` functions, which are not
JSON and would have to be injected into the page as executable script. The
numeric layer under them (`ellipsoid_surface`, `cylinder_surface`,
`disc_surface`, `param_surface`) returns plain arrays, so real JSON is built
from those. Shapes therefore look the same in the interface as in the manual,
without the injection.

**Read-back only claims what it can prove.** Beyond the structural checks,
before accepting a construct the reader renders it exactly as it would be saved
and compares against the source it came from. If the two differ, the original
text is kept and the reason is reported. Parsing something is not the same as
understanding it, and the difference is exactly where silent damage would come
from.

## The protocol

One JSON object per line, both directions:

```json
{"id": 7, "op": "traces", "payload": {"expr": "Spheroid(0.4)", "cutaway": true}}
{"id": 7, "ok": true, "result": {"data": [...], "layout": {...}}}
```

Ops: `ping`, `catalogue`, `traces`, `parse`, `run`. The first line the sidecar
ever writes is `{"event": "ready", ...}`.

`run` executes in a fresh anonymous module with stdout redirected to a
temporary file. The redirect has to be undone however the script ends —
otherwise the protocol's own replies would vanish into the capture — so a
timeout *interrupts* the task rather than abandoning it, and the reply says
`wedged` if the task refused to unwind, which makes the client restart.

## Tests

```bash
python3 tests/test_studio.py          # no Julia needed
python3 tests/test_studio.py --julia  # adds the sidecar-backed tests
```

The load-bearing one is preservation: every script under `scripts/` is opened
and written back, and no line may be lost.

## Known limits

- Sensitivities, laminates and custom/FE/neural inclusions are not modelled.
  Scripts using them open and are preserved, but those parts are not editable.
- Viscoelastic properties are entered as Julia expressions rather than through
  a form.
- An inner `Homogenized` cannot sit inside an ALV chain (MeanFieldHom cannot
  re-express a homogenized result as a `ViscoLaw`); the interface blocks the
  combination.
- `web/vendor/plotly.min.js` is the official *plotly-gl3d* partial bundle,
  vendored so the interface works with no outbound network.
