"""The HTTP server.

Deliberately built on `http.server` from the standard library rather than a
framework: the interface must start with `python3 -m mfhstudio` on any machine
that already runs MeanFieldHom, and adding a dependency tree to achieve that
would be a poor trade. The API is small — a handful of JSON endpoints and the
static files.
"""

from __future__ import annotations

import json
import os
import threading
import traceback
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from typing import Any, Optional
from urllib.parse import urlparse

from .codegen import extract_embedded, generate
from .juliabridge import Bridge, SidecarError, SidecarUnavailable
from .model import Model, default_model
from .readback import model_from_script

WEB_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "web"))

_MIME = {
    ".html": "text/html; charset=utf-8",
    ".js": "text/javascript; charset=utf-8",
    ".css": "text/css; charset=utf-8",
    ".json": "application/json",
    ".svg": "image/svg+xml",
    ".map": "application/json",
}


class Session:
    """Everything one running studio holds.

    The model lives here, in Python, which is what makes the sidecar
    disposable: it can be restarted after a wedge without losing work.
    """

    def __init__(self) -> None:
        self.model: Model = default_model()
        self.path: Optional[str] = None
        self.bridge = Bridge()
        self.lock = threading.Lock()
        self._catalogue: Optional[dict] = None

    def catalogue(self) -> dict:
        if self._catalogue is None:
            self._catalogue = self.bridge.catalogue()
        return self._catalogue

    def script(self) -> str:
        return generate(self.model)


class Handler(BaseHTTPRequestHandler):
    session: Session  # injected by `serve`

    server_version = "MFHStudio"

    def log_message(self, fmt: str, *args) -> None:  # quieter console
        return

    # -- helpers ----------------------------------------------------------

    def _send(self, code: int, body: bytes, ctype: str) -> None:
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        try:
            self.wfile.write(body)
        except BrokenPipeError:
            pass

    def _json(self, obj: Any, code: int = 200) -> None:
        self._send(code, json.dumps(obj).encode("utf-8"), "application/json")

    def _error(self, exc: Exception, code: int = 500) -> None:
        # The interface shows these verbatim: a precise message beats a
        # generic failure the user cannot act on.
        payload = {"error": str(exc), "kind": type(exc).__name__}
        if not isinstance(exc, (SidecarError, SidecarUnavailable, ValueError)):
            payload["traceback"] = traceback.format_exc()
        self._json(payload, code)

    def _body(self) -> dict:
        n = int(self.headers.get("Content-Length") or 0)
        if n <= 0:
            return {}
        return json.loads(self.rfile.read(n).decode("utf-8"))

    # -- routing ----------------------------------------------------------

    def do_GET(self) -> None:  # noqa: N802
        path = urlparse(self.path).path
        try:
            if path == "/api/state":
                return self._json(self._state())
            if path == "/api/catalogue":
                return self._json(self.session.catalogue())
            if path == "/api/script":
                return self._json({"source": self.session.script()})
            if path == "/api/sidecar":
                return self._json(self.session.bridge.status)
            return self._static(path)
        except Exception as exc:  # noqa: BLE001
            return self._error(exc)

    def do_POST(self) -> None:  # noqa: N802
        path = urlparse(self.path).path
        try:
            body = self._body()
            if path == "/api/model":
                return self._set_model(body)
            if path == "/api/traces":
                return self._traces(body)
            if path == "/api/run":
                return self._run(body)
            if path == "/api/open":
                return self._open(body)
            if path == "/api/save":
                return self._save(body)
            if path == "/api/sidecar/restart":
                self.session.bridge.restart()
                self.session._catalogue = None
                return self._json(self.session.bridge.status)
            return self._json({"error": f"unknown endpoint {path}"}, 404)
        except Exception as exc:  # noqa: BLE001
            return self._error(exc)

    # -- endpoints --------------------------------------------------------

    def _state(self) -> dict:
        s = self.session
        return {
            "model": s.model.to_dict(),
            "script": s.script(),
            "problems": s.model.validate(),
            "path": s.path,
            "multiscale": s.model.uses_multiscale(),
        }

    def _set_model(self, body: dict) -> None:
        with self.session.lock:
            self.session.model = Model.from_dict(body.get("model", {}))
        return self._json(self._state())

    def _traces(self, body: dict) -> None:
        expr = body.get("expr") or ""
        if not expr:
            return self._json({"data": [], "layout": {}})
        kw = {}
        if "cutaway" in body:
            kw["cutaway"] = bool(body["cutaway"])
        return self._json(self.session.bridge.traces(expr, **kw))

    def _run(self, body: dict) -> None:
        source = body.get("source") or self.session.script()
        timeout = float(body.get("timeout", 300.0))
        result = self.session.bridge.run(source, timeout=timeout)
        # A wedged sidecar cannot be trusted for the next request.
        if result.get("wedged"):
            try:
                self.session.bridge.restart()
                self.session._catalogue = None
            except Exception:
                pass
        return self._json(result)

    def _open(self, body: dict) -> None:
        path = body.get("path") or ""
        source = body.get("source")
        if source is None:
            if not os.path.isfile(path):
                raise ValueError(f"no such file: {path}")
            with open(path, encoding="utf-8") as fh:
                source = fh.read()

        model, report = model_from_script(source, self.session.bridge)
        with self.session.lock:
            self.session.model = model
            self.session.path = path or None
        out = self._state()
        out["read_report"] = report
        return self._json(out)

    def _save(self, body: dict) -> None:
        path = body.get("path") or self.session.path
        if not path:
            raise ValueError("no path given")
        source = body.get("source") or self.session.script()
        os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(source)
        with self.session.lock:
            self.session.path = path
        return self._json({"path": path, "bytes": len(source)})

    # -- static files -----------------------------------------------------

    def _static(self, path: str) -> None:
        rel = "index.html" if path in ("/", "") else path.lstrip("/")
        full = os.path.abspath(os.path.join(WEB_DIR, rel))
        # Never serve outside the web directory.
        if not full.startswith(WEB_DIR) or not os.path.isfile(full):
            return self._send(404, b"not found", "text/plain; charset=utf-8")
        ext = os.path.splitext(full)[1]
        with open(full, "rb") as fh:
            data = fh.read()
        return self._send(200, data, _MIME.get(ext, "application/octet-stream"))


def serve(host: str = "127.0.0.1", port: int = 8765, open_browser: bool = True) -> None:
    session = Session()
    handler = type("BoundHandler", (Handler,), {"session": session})
    httpd = ThreadingHTTPServer((host, port), handler)

    url = f"http://{host}:{port}/"
    print(f"MFH Studio → {url}")
    print("Starting the Julia sidecar (MeanFieldHom takes ~10 s to load)…")

    def warm() -> None:
        try:
            session.bridge.start()
            session.catalogue()
            print("Sidecar ready.")
        except Exception as exc:  # noqa: BLE001
            print(f"Sidecar unavailable: {exc}")
            print("The interface still runs; 3-D, read-back and execution are off.")

    threading.Thread(target=warm, daemon=True).start()

    if open_browser:
        import webbrowser

        threading.Timer(0.5, lambda: webbrowser.open(url)).start()

    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print("\nStopping.")
    finally:
        session.bridge.stop()
        httpd.server_close()
