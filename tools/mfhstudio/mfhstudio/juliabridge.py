"""The Python side of the Julia sidecar.

Loading MeanFieldHomogenization costs about ten seconds, so one long-lived Julia process
serves every request over JSON lines. The bridge owns that process: it starts
it lazily, serializes requests, and restarts it if it dies or wedges — the
model itself lives in Python, so a restart loses nothing but time.
"""

from __future__ import annotations

import json
import os
import queue
import shutil
import subprocess
import threading
import time
from dataclasses import dataclass, field
from typing import Any, Optional


class SidecarError(RuntimeError):
    """The sidecar answered, and the answer was a failure."""


class SidecarUnavailable(RuntimeError):
    """The sidecar could not be started or has stopped responding."""


def _diagnose(log: str) -> str:
    """Turn a Julia start-up failure into something the user can act on.

    A raw stack trace tells you what happened but not what to do, and the two
    common failures here have a one-line remedy each.
    """
    hint = None
    if "does not seem to be installed" in log or "Pkg.instantiate()" in log:
        hint = (
            "the MeanFieldHomogenization project has not been instantiated. The sidecar "
            "now does this itself on start-up; if you are seeing this, run it "
            "by hand:\n"
            "    julia --project=<MeanFieldHomogenization.jl> -e 'using Pkg; Pkg.instantiate()'"
        )
    elif "MeanFieldHomogenization" in log and "not found in current path" in log:
        hint = (
            "Julia started but could not find MeanFieldHomogenization. Check that "
            "tools/mfhstudio sits inside the package checkout, or point the "
            "sidecar at it explicitly."
        )
    elif "UndefVarError" in log or "LoadError" in log:
        hint = "the sidecar scripts failed to load; the Julia error is below."

    head = "Julia exited during start-up"
    if hint:
        head += ": " + hint
    return head + "\n\n" + log


HERE = os.path.dirname(os.path.abspath(__file__))
SIDECAR_JL = os.path.join(HERE, "..", "julia", "sidecar.jl")
PROJECT_ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))


def find_julia() -> Optional[str]:
    """The Julia executable, honoring `JULIA` if the user set one."""
    return os.environ.get("JULIA") or shutil.which("julia")


@dataclass
class Bridge:
    julia: Optional[str] = None
    project: str = PROJECT_ROOT
    script: str = os.path.abspath(SIDECAR_JL)
    # A first run may have to instantiate the project, which can take minutes.
    # Waiting is safe: a sidecar that *dies* is detected at once, so this only
    # bounds the "alive but slow" case.
    boot_timeout: float = 900.0

    _proc: Optional[subprocess.Popen] = field(default=None, init=False, repr=False)
    # `_lock` serializes round-trips (the sidecar handles one request at a
    # time); `_boot_lock` serializes start-up.  They are distinct on purpose:
    # start-up can take minutes, and holding the round-trip lock for it would
    # make every queued request inherit that wait as if it were its own.
    _lock: threading.Lock = field(default_factory=threading.Lock, init=False, repr=False)
    _boot_lock: threading.Lock = field(
        default_factory=threading.Lock, init=False, repr=False
    )
    # One reply queue PER pending request id.  A single shared queue lets
    # concurrent waiters steal each other's messages: whoever calls `get`
    # first receives whatever arrives and drops it if the id does not match,
    # which silently loses the readiness event to an in-flight request.
    _pending: dict = field(default_factory=dict, init=False, repr=False)
    _ready_evt: threading.Event = field(
        default_factory=threading.Event, init=False, repr=False
    )
    _reader: Optional[threading.Thread] = field(default=None, init=False, repr=False)
    _next_id: int = field(default=0, init=False, repr=False)
    _ready: bool = field(default=False, init=False, repr=False)
    _boot_error: Optional[str] = field(default=None, init=False, repr=False)
    _stderr: list = field(default_factory=list, init=False, repr=False)

    # -- lifecycle --------------------------------------------------------

    @property
    def status(self) -> dict:
        alive = self._proc is not None and self._proc.poll() is None
        return {
            "running": alive,
            "ready": self._ready,
            "error": self._boot_error,
            "julia": self.julia or find_julia(),
        }

    def start(self) -> None:
        """Launch the sidecar and wait for its readiness line.

        Idempotent and safe to call from several threads: the first caller
        boots, the others wait for the same readiness event instead of
        starting a second Julia or — worse — proceeding as if the first one
        were already usable.
        """
        with self._boot_lock:
            if self._ready and self._proc is not None and self._proc.poll() is None:
                return
            self._start_locked()

    def _start_locked(self) -> None:
        if self._proc is not None and self._proc.poll() is None and self._ready:
            return
        exe = self.julia or find_julia()
        if exe is None:
            raise SidecarUnavailable(
                "no `julia` on PATH. Install Julia, or point the JULIA "
                "environment variable at the executable."
            )
        if not os.path.isfile(self.script):
            raise SidecarUnavailable(f"sidecar script not found: {self.script}")

        self._ready = False
        self._boot_error = None
        self._stderr = []
        self._ready_evt.clear()
        self._pending = {}
        self._proc = subprocess.Popen(
            [exe, f"--project={self.project}", "--startup-file=no", self.script],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            # Julia speaks UTF-8; `text=True` alone decodes with the locale
            # encoding, which on Windows is cp1252 — enough to turn a `φ` in
            # the script's output into `Ï†` on its way back through the pipe.
            encoding="utf-8",
            errors="replace",
            bufsize=1,
        )
        self._reader = threading.Thread(target=self._read_loop, daemon=True)
        self._reader.start()
        threading.Thread(target=self._drain_stderr, daemon=True).start()

        # The readiness event is routed to us by the reader thread, so an
        # in-flight request can no longer consume it.
        deadline = time.time() + self.boot_timeout
        while time.time() < deadline:
            if self._ready_evt.wait(0.2):
                if not self._ready:
                    raise SidecarUnavailable(
                        "MeanFieldHomogenization failed to load in the sidecar:\n"
                        + str(self._boot_error)
                    )
                return
            if self._proc.poll() is not None:
                raise SidecarUnavailable(_diagnose("".join(self._stderr[-60:])))
        raise SidecarUnavailable(
            f"sidecar did not become ready within {self.boot_timeout:.0f} s. "
            "A first run on a fresh checkout has to instantiate and precompile "
            "the whole dependency tree; re-running is cheaper because that "
            "work is cached. Run `--check` for a diagnosis."
        )

    def stop(self) -> None:
        proc, self._proc = self._proc, None
        self._ready = False
        self._ready_evt.clear()
        self._pending = {}
        if proc is None:
            return
        try:
            proc.stdin and proc.stdin.close()
            proc.wait(timeout=5)
        except Exception:
            proc.kill()

    def restart(self) -> None:
        self.stop()
        self.start()

    # -- plumbing ---------------------------------------------------------

    def _read_loop(self) -> None:
        proc = self._proc
        assert proc is not None and proc.stdout is not None
        for line in proc.stdout:
            line = line.strip()
            if not line:
                continue
            try:
                msg = json.loads(line)
            except json.JSONDecodeError:
                # Anything the sidecar prints outside the protocol is kept for
                # diagnostics rather than silently dropped.
                self._stderr.append(line + "\n")
                continue
            if msg.get("event") == "ready":
                self._ready = bool(msg.get("ok"))
                self._boot_error = msg.get("error")
                self._ready_evt.set()
                continue
            # Route to the waiter that asked for it.  A reply with no pending
            # waiter is stale (its caller timed out, or a restart intervened)
            # and is dropped rather than handed to whoever asks next.
            q = self._pending.get(msg.get("id"))
            if q is not None:
                q.put(msg)

    def _drain_stderr(self) -> None:
        proc = self._proc
        if proc is None or proc.stderr is None:
            return
        for line in proc.stderr:
            self._stderr.append(line)
            del self._stderr[:-200]

    def call(self, op: str, payload: Optional[dict] = None, timeout: float = 300.0) -> Any:
        """One request, one reply.

        `timeout` bounds the time the sidecar spends *answering*, so the wait
        for it to finish booting is taken first and separately.  Sending on a
        process that is alive but still loading MeanFieldHomogenization would spend the
        whole budget waiting for a reply that cannot come yet, and report a
        wedged sidecar when nothing is wrong but a cold precompilation.
        """
        # Not under `_lock`: booting can take minutes and is shared work.
        self.start()

        with self._lock:
            if self._proc is None or self._proc.poll() is not None or not self._ready:
                self.start()
            assert self._proc is not None and self._proc.stdin is not None

            self._next_id += 1
            rid = self._next_id
            inbox: "queue.Queue[dict]" = queue.Queue()
            self._pending[rid] = inbox
            req = {"id": rid, "op": op, "payload": payload or {}}
            try:
                try:
                    self._proc.stdin.write(json.dumps(req) + "\n")
                    self._proc.stdin.flush()
                except (BrokenPipeError, ValueError) as exc:
                    raise SidecarUnavailable(f"sidecar pipe closed: {exc}") from exc

                deadline = time.time() + timeout
                while time.time() < deadline:
                    try:
                        msg = inbox.get(timeout=0.2)
                    except queue.Empty:
                        if self._proc.poll() is not None:
                            raise SidecarUnavailable(
                                "sidecar died:\n" + "".join(self._stderr[-40:])
                            )
                        continue
                    if not msg.get("ok"):
                        raise SidecarError(msg.get("error") or "unknown sidecar error")
                    return msg.get("result")
            finally:
                self._pending.pop(rid, None)

        raise SidecarUnavailable(
            f"no reply to `{op}` within {timeout:.0f} s; the sidecar may be wedged"
        )

    # -- the operations ---------------------------------------------------

    def ping(self) -> dict:
        return self.call("ping", timeout=30)

    def catalog(self) -> dict:
        return self.call("catalog", timeout=120)

    def traces(self, expr: str, **kw) -> dict:
        return self.call("traces", {"expr": expr, **kw}, timeout=60)

    def parse(self, source: str) -> dict:
        return self.call("parse", {"source": source}, timeout=60)

    def run(self, source: str, timeout: float = 300.0) -> dict:
        # The sidecar enforces its own limit; the outer one is slightly longer
        # so a clean "timed out" answer wins over a transport error.
        return self.call(
            "run", {"source": source, "timeout": timeout}, timeout=timeout + 30
        )
