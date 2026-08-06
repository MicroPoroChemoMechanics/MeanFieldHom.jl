"""Opening an Echoes script converts it.

`tools/echoes2mfh` already turns an Echoes Python script into idiomatic Julia,
and MFH Studio already reads Julia back into a model. Putting a "convert"
button next to "open" would make the user decide which of the two they wanted;
the file extension already says it. So **Open** takes a `.py`, translates it,
asks where to put the `.jl`, and opens that — one action, one mental step.

The translator is imported rather than shelled out to: it is a sibling package
under `tools/`, and importing it keeps its findings as objects instead of text
that would have to be parsed back.
"""

from __future__ import annotations

import os
import sys
from dataclasses import dataclass, field
from typing import Optional

# `tools/echoes2mfh` sits beside `tools/mfhstudio`; neither is installed.
_TOOLS = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
_ECHOES2MFH = os.path.join(_TOOLS, "echoes2mfh")
if _ECHOES2MFH not in sys.path:
    sys.path.insert(0, _ECHOES2MFH)


class ConversionUnavailable(RuntimeError):
    """The translator is not importable from here."""


@dataclass
class Conversion:
    source: str
    julia: str = ""
    findings: list = field(default_factory=list)
    blocking: int = 0

    @property
    def clean(self) -> bool:
        return not self.findings

    def summary(self) -> str:
        if self.clean:
            return "Translated with no findings."
        n = len(self.findings)
        return (
            f"Translated with {n} finding{'s' if n != 1 else ''} "
            f"({self.blocking} blocking). Each is marked UNTRANSLATED in the "
            f"script and raises at run time, so a half-translated model cannot "
            f"quietly produce numbers."
        )


def available() -> bool:
    try:
        import echoes2mfh.emit  # noqa: F401
        import echoes2mfh.extract  # noqa: F401
    except Exception:  # noqa: BLE001
        return False
    return True


def is_echoes_script(path: str) -> bool:
    return os.path.splitext(path)[1].lower() == ".py"


def suggest_output(path: str) -> str:
    """Where the translation would go if the user does not say otherwise."""
    base = os.path.splitext(os.path.basename(path))[0]
    return os.path.join(os.path.dirname(os.path.abspath(path)), base + ".jl")


def convert(path: str, source: Optional[str] = None) -> Conversion:
    """Translate one Echoes script. Raises with an actionable message."""
    try:
        from echoes2mfh.emit import emit
        from echoes2mfh.extract import extract
    except Exception as exc:  # noqa: BLE001
        raise ConversionUnavailable(
            "the echoes2mfh translator could not be imported from "
            f"{_ECHOES2MFH}: {exc}"
        ) from exc

    if source is None:
        with open(path, encoding="utf-8", errors="replace") as fh:
            source = fh.read()

    import ast

    try:
        ast.parse(source)
    except SyntaxError as exc:
        # Python 2 is the one failure with a one-line remedy.
        raise ValueError(
            f"{os.path.basename(path)} does not parse as Python 3: {exc.msg} "
            f"(line {exc.lineno}). If it is Python 2, run `2to3 -w` on it first."
        ) from exc

    name = os.path.splitext(os.path.basename(path))[0]
    script = extract(source, path)
    julia = emit(script, name)

    return Conversion(
        source=path,
        julia=julia,
        findings=[
            {
                "line": f.lineno,
                "severity": f.severity,
                "reason": f.reason,
                "suggestion": f.suggestion,
            }
            for f in script.findings
        ],
        blocking=script.finding_count("blocking"),
    )
