"""MFH Studio — a graphical builder for MeanFieldHom scripts.

    python3 -m mfhstudio                 start and open a browser
    python3 -m mfhstudio --port 9000     pick the port
    python3 -m mfhstudio --no-browser    stay in the terminal
    python3 -m mfhstudio --check         verify the sidecar and exit
"""

from __future__ import annotations

import argparse
import sys


def main(argv: list | None = None) -> int:
    p = argparse.ArgumentParser(prog="mfhstudio", description=__doc__)
    p.add_argument("--host", default="127.0.0.1")
    p.add_argument("--port", type=int, default=8765)
    p.add_argument("--no-browser", action="store_true")
    p.add_argument(
        "--check", action="store_true",
        help="start the Julia sidecar, report what it found, and exit",
    )
    args = p.parse_args(argv)

    if args.check:
        return _check()

    from .server import serve

    serve(args.host, args.port, open_browser=not args.no_browser)
    return 0


def _check() -> int:
    from .juliabridge import Bridge, SidecarUnavailable

    b = Bridge()
    try:
        b.start()
    except SidecarUnavailable as exc:
        print(f"sidecar unavailable: {exc}", file=sys.stderr)
        return 1
    try:
        cat = b.catalogue()
        print(f"MeanFieldHom {cat['mfh_version']} on Julia {cat['julia_version']}")
        print(f"  schemes    : {len(cat['schemes'])}")
        print(f"  geometries : {len(cat['geometries'])}")
        print(f"  interfaces : {len(cat['interfaces'])}")
        print(f"  lenses     : {len(cat['lenses'])}")
        return 0
    finally:
        b.stop()


if __name__ == "__main__":
    sys.exit(main())
