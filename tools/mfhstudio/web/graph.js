/* The scale graph — draggable boxes joined by connectors.
 *
 * A multiscale model in MeanFieldHomogenization is a directed graph: a phase property may
 * hold `Homogenized(inner_cell, scheme)`, and the outer scheme resolves the
 * inner scale when it reads that key. A form full of dropdowns hides that
 * shape; a picture is the shape.
 *
 * One box per scale. Each box has an output port on its right (its effective
 * property) and one input port per phase property slot. Dragging from an
 * output to an input creates the seam; the generated Julia updates as the
 * connector lands.
 *
 * Everything here manipulates `S.model` and calls `push()`, exactly like the
 * form panels — the graph is a second view on the same model, never a
 * separate source of truth.
 */

(function () {
  const NS = "http://www.w3.org/2000/svg";
  const svgEl = (tag, attrs = {}) => {
    const n = document.createElementNS(NS, tag);
    for (const [k, v] of Object.entries(attrs)) n.setAttribute(k, v);
    return n;
  };

  /** Drag state: either a box being moved or a connector being pulled. */
  let drag = null;

  function host() {
    return document.querySelector("#graph");
  }

  /* ── geometry helpers ───────────────────────────────────────────── */

  const portId = (cellId, phaseName, key) => `${cellId}|${phaseName}|${key}`;

  function centerOf(node) {
    const h = host().getBoundingClientRect();
    const r = node.getBoundingClientRect();
    return {
      x: r.left - h.left + r.width / 2 + host().scrollLeft,
      y: r.top - h.top + r.height / 2 + host().scrollTop,
    };
  }

  /** A cubic that leaves horizontally, so connectors read left-to-right. */
  function curve(a, b) {
    const dx = Math.max(40, Math.abs(b.x - a.x) * 0.5);
    return `M ${a.x} ${a.y} C ${a.x + dx} ${a.y}, ${b.x - dx} ${b.y}, ${b.x} ${b.y}`;
  }

  /* ── drawing ────────────────────────────────────────────────────── */

  function drawGraph() {
    const g = host();
    if (!g) return;
    g.replaceChildren();

    const svg = svgEl("svg", { class: "edges" });
    const defs = svgEl("defs");
    const marker = svgEl("marker", {
      id: "arrow", viewBox: "0 0 10 10", refX: "9", refY: "5",
      markerWidth: "7", markerHeight: "7", orient: "auto-start-reverse",
    });
    marker.append(svgEl("path", { d: "M 0 0 L 10 5 L 0 10 z", class: "arrowhead" }));
    defs.append(marker);
    svg.append(defs);
    g.append(svg);

    for (const c of S.model.cells) g.append(cellBox(c));

    // Edges are drawn after the boxes so the ports have real positions.
    requestAnimationFrame(() => drawEdges(svg));
  }

  function drawEdges(svg) {
    const g = host();
    if (!g || !svg.isConnected) return;
    [...svg.querySelectorAll("path.edge, path.edge-hit")].forEach((p) => p.remove());

    const byId = Object.fromEntries(S.model.cells.map((c) => [c.id, c]));
    for (const c of S.model.cells) {
      for (const ph of c.phases) {
        for (const pr of ph.properties || []) {
          if (pr.source !== "cell" || !byId[pr.cell]) continue;
          const from = g.querySelector(`[data-out="${pr.cell}"]`);
          const to = g.querySelector(`[data-in="${CSS.escape(portId(c.id, ph.name, pr.key))}"]`);
          if (!from || !to) continue;
          const d = curve(centerOf(from), centerOf(to));
          svg.append(svgEl("path", { class: "edge", d, "marker-end": "url(#arrow)" }));
          // A fat invisible copy makes the thin curve clickable.
          const hit = svgEl("path", { class: "edge-hit", d });
          hit.addEventListener("click", () => {
            pr.source = "builder";
            pr.cell = null;
            push();
          });
          hit.addEventListener("mouseenter", () => hit.classList.add("hot"));
          hit.addEventListener("mouseleave", () => hit.classList.remove("hot"));
          svg.append(hit);
        }
      }
    }
    sizeCanvas(svg);
  }

  function sizeCanvas(svg) {
    const g = host();
    let w = 0, h = 0;
    for (const box of g.querySelectorAll(".node")) {
      w = Math.max(w, box.offsetLeft + box.offsetWidth + 40);
      h = Math.max(h, box.offsetTop + box.offsetHeight + 40);
    }
    g.style.minHeight = Math.max(h, 220) + "px";
    svg.setAttribute("width", Math.max(w, g.clientWidth));
    svg.setAttribute("height", Math.max(h, 220));
  }

  function cellBox(c) {
    const box = document.createElement("div");
    box.className = "node" + (c.id === S.cellId ? " on" : "");
    box.style.left = (c.ui?.x ?? 40) + "px";
    box.style.top = (c.ui?.y ?? 40) + "px";
    box.dataset.cell = c.id;

    const head = document.createElement("header");
    head.append(Object.assign(document.createElement("b"), { textContent: c.name }));
    const pick = document.createElement("button");
    pick.className = "small";
    pick.textContent = "edit";
    pick.title = "Edit this scale in the panel";
    pick.addEventListener("click", (e) => {
      e.stopPropagation();
      S.cellId = c.id;
      S.phaseIdx = 0;
      render();
      draw3d();
    });
    head.append(pick);
    box.append(head);

    // Dragging the header moves the box; the position rides in the model.
    head.addEventListener("pointerdown", (e) => {
      if (e.target.closest("button")) return;
      const h = host().getBoundingClientRect();
      drag = {
        kind: "move", cell: c, box,
        dx: e.clientX - h.left - box.offsetLeft,
        dy: e.clientY - h.top - box.offsetTop,
      };
      head.setPointerCapture(e.pointerId);
      e.preventDefault();
    });

    for (const ph of c.phases) {
      const row = document.createElement("div");
      row.className = "phase";
      const tag = document.createElement("span");
      tag.className = "pill" + (ph.is_matrix ? " matrix" : "");
      tag.textContent = ph.name;
      row.append(tag);
      for (const pr of ph.properties || []) {
        row.append(inputPort(c, ph, pr));
      }
      box.append(row);
    }

    // The output port: this scale's effective property.
    const out = document.createElement("div");
    out.className = "port out";
    out.dataset.out = c.id;
    out.title = "Drag from here onto a property of another scale";
    out.addEventListener("pointerdown", (e) => {
      e.stopPropagation();
      const svg = host().querySelector("svg.edges");
      const path = svgEl("path", { class: "edge pending" });
      svg.append(path);
      drag = { kind: "link", from: c, path, start: centerOf(out) };
      out.setPointerCapture(e.pointerId);
      e.preventDefault();
    });
    box.append(out);

    return box;
  }

  function inputPort(c, ph, pr) {
    const p = document.createElement("span");
    const seam = pr.source === "cell";
    p.className = "port in" + (seam ? " linked" : "");
    p.dataset.in = portId(c.id, ph.name, pr.key);
    p.textContent = pr.key;
    p.title = seam
      ? "Fed by another scale — click the connector to detach"
      : "Drop a scale's output here to feed this property";
    return p;
  }

  /* ── pointer handling ───────────────────────────────────────────── */

  document.addEventListener("pointermove", (e) => {
    if (!drag) return;
    const g = host();
    const h = g.getBoundingClientRect();

    if (drag.kind === "move") {
      const x = Math.max(0, e.clientX - h.left - drag.dx);
      const y = Math.max(0, e.clientY - h.top - drag.dy);
      drag.box.style.left = x + "px";
      drag.box.style.top = y + "px";
      drag.cell.ui = { x: Math.round(x), y: Math.round(y) };
      drawEdges(g.querySelector("svg.edges"));
      drag.moved = true;
      return;
    }

    if (drag.kind === "link") {
      const pt = { x: e.clientX - h.left + g.scrollLeft, y: e.clientY - h.top + g.scrollTop };
      drag.path.setAttribute("d", curve(drag.start, pt));
      const over = document.elementFromPoint(e.clientX, e.clientY);
      g.querySelectorAll(".port.in.target").forEach((n) => n.classList.remove("target"));
      if (over && over.classList.contains("in")) over.classList.add("target");
    }
  });

  document.addEventListener("pointerup", (e) => {
    if (!drag) return;
    const g = host();

    if (drag.kind === "move") {
      const moved = drag.moved;
      drag = null;
      // Only round-trip to the server when something actually changed.
      if (moved) push();
      return;
    }

    if (drag.kind === "link") {
      const from = drag.from;
      drag.path.remove();
      g.querySelectorAll(".port.in.target").forEach((n) => n.classList.remove("target"));
      const over = document.elementFromPoint(e.clientX, e.clientY);
      drag = null;
      if (!over || !over.dataset.in) return;

      const [cellId, phaseName, key] = over.dataset.in.split("|");
      if (cellId === from.id) {
        toast("A scale cannot feed itself.", true);
        return;
      }
      const target = S.model.cells.find((c) => c.id === cellId);
      const ph = target && target.phases.find((p) => p.name === phaseName);
      const pr = ph && (ph.properties || []).find((x) => x.key === key);
      if (!pr) return;

      pr.source = "cell";
      pr.cell = from.id;
      pr.scheme = pr.scheme || "MoriTanaka";
      pr.scheme_options = pr.scheme_options || {};
      // The server validates the whole graph and reports a cycle if this
      // connection closed one; nothing is silently accepted here.
      push();
    }
  });

  window.drawGraph = drawGraph;
})();
