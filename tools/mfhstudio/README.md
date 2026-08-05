# MFH Studio

A local web interface for building MeanFieldHom scripts: draw the shape of each
phase, run the model, and read existing scripts back without damaging them.

```bash
cd tools/mfhstudio
python3 -m mfhstudio                 # start and open a browser
python3 -m mfhstudio --port 9000     # pick the port
python3 -m mfhstudio --no-browser    # stay in the terminal
python3 -m mfhstudio --check         # verify the Julia side and exit
python3 -m mfhstudio --project @env  # use another Julia environment
```

On Windows the command is `python -m mfhstudio`.

Requires Python 3.10+ (standard library only — no `pip install`) and a `julia`
on `PATH` able to load MeanFieldHom. Set `JULIA` to point at a specific
executable.

### Getting the Julia side to start

`--check` diagnoses this without paying the ten-second load, and is the right
first command on a new machine.

Two things go wrong in practice.

**The environment has never been instantiated.** The sidecar's `using JSON3`
then dies with a stack trace whose only advice is to run `Pkg.instantiate()`.
The sidecar now runs it itself on first start, which takes a few minutes once.
By hand:

    julia --project=<MeanFieldHom.jl> -e 'using Pkg; Pkg.instantiate()'

**The committed `Manifest.toml` pins dependencies to sibling checkouts.** It
records, for instance:

    [[deps.TensND]]
    path = "../TensND.jl"

which only resolves when `TensND.jl` sits next to `MeanFieldHom.jl`. That is a
development override: MeanFieldHom is not registered, but its dependencies are,
so on a machine with only the MeanFieldHom clone the way out is a separate
environment that takes them from the registry:

    julia -e 'using Pkg; Pkg.activate("mfhstudio", shared=true); \
              Pkg.add("TensND"); Pkg.develop(path=raw"<MeanFieldHom.jl>")'

then `python -m mfhstudio --project @mfhstudio`. Adding the dependencies before
the develop is what stops the resolver from turning them back into path
entries. `--check` reads the manifest and names the missing checkout itself, so
it will tell you which packages to add.

If Julia is unavailable for any reason the interface still comes up: you can
build and save a script, and a banner says what is off. Only the 3-D view,
reading a script back, and Run need the sidecar.

The user-facing guide, with screenshots, is
`docs/src/tools/mfhstudio.md`. This file covers how the thing is put together.

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
| `mfhstudio/catalog.py` | the form definitions — available without Julia |
| `mfhstudio/server.py` | HTTP endpoints, session state, the file browser |
| `julia/sidecar.jl` | the JSON-lines loop |
| `julia/introspect.jl` | the feature catalog, read from the live package |
| `julia/geometry.jl` | 3-D traces |
| `julia/parse_script.jl` | `Meta.parse` → nodes that tile the file exactly |
| `web/graph.js` | the draggable scale graph |
| `web/picker.js` | the file dialog |

## Three design decisions worth knowing

**The catalog has two halves.** Form definitions — which fields a spheroid
needs, which lenses exist — are interface concerns and live in Python, so they
are there immediately. Only the scheme list and each scheme's solver options
come from the sidecar, and they replace the fallback wholesale rather than
merging, so a scheme MeanFieldHom drops disappears. Putting the forms behind
the sidecar is what once made every control dead on a machine where Julia
failed to start.

**The scheme list is introspected, never hard-coded.** Schemes come from
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

Ops: `ping`, `catalog`, `traces`, `parse`, `run`. The first line the sidecar
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

- Sensitivities, laminates and custom/FE/neural inclusions are not modeled.
  Scripts using them open and are preserved, but those parts are not editable.
- Viscoelastic properties are entered as Julia expressions rather than through
  a form.

- An inner `Homogenized` cannot sit inside an ALV chain (MeanFieldHom cannot
  re-express a homogenized result as a `ViscoLaw`); the interface blocks the
  combination.
- `web/vendor/plotly.min.js` is the official *plotly-gl3d* partial bundle,
  vendored so the interface works with no outbound network.
