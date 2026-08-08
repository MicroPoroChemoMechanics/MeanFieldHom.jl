"""Check the Julia environment before paying ten seconds to find out.

MeanFieldHomogenization's committed `Manifest.toml` records some dependencies by
*relative path* to a sibling checkout:

    [[deps.TensND]]
    path = "../TensND.jl"

which resolves only when those siblings are actually there. On a machine where
someone cloned MeanFieldHomogenization.jl alone, `Pkg.instantiate()` fails and the sidecar
dies with a message about a package that "does not seem to be installed" —
true, but not the useful part.

Reading the manifest here costs nothing and names the missing checkout and the
directory it is expected in, which is what the user needs.
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass, field

try:  # Python 3.11+
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - 3.10
    tomllib = None


@dataclass
class Report:
    project: str
    ok: bool = True
    problems: list = field(default_factory=list)
    hints: list = field(default_factory=list)
    missing_paths: list = field(default_factory=list)

    def fail(self, problem: str, hint: str = "") -> None:
        self.ok = False
        self.problems.append(problem)
        if hint:
            self.hints.append(hint)

    def text(self) -> str:
        if self.ok:
            return f"Julia environment looks usable: {self.project}"
        out = [f"The Julia environment at {self.project} cannot be used:"]
        out += [f"  · {p}" for p in self.problems]
        if self.hints:
            out.append("")
            out += [f"  {h}" for h in self.hints]
        return "\n".join(out)


def _path_deps(manifest: str) -> list:
    """Every `path = "..."` dependency, as (name, path) pairs."""
    if tomllib is not None:
        try:
            with open(manifest, "rb") as fh:
                data = tomllib.load(fh)
        except Exception:  # noqa: BLE001
            return []
        deps = data.get("deps", data)
        out = []
        for name, entries in deps.items():
            if not isinstance(entries, list):
                continue
            for e in entries:
                if isinstance(e, dict) and "path" in e:
                    out.append((name, e["path"]))
        return out

    # 3.10 fallback: the format is regular enough to read without a parser.
    out = []
    name = None
    for line in open(manifest, encoding="utf-8", errors="replace"):
        m = re.match(r"\[\[deps\.([^\]]+)\]\]", line.strip())
        if m:
            name = m.group(1)
            continue
        m = re.match(r'path\s*=\s*"([^"]+)"', line.strip())
        if m and name:
            out.append((name, m.group(1)))
    return out


def check_project(project: str) -> Report:
    project = os.path.abspath(project)
    rep = Report(project=project)

    proj_toml = os.path.join(project, "Project.toml")
    if not os.path.isfile(proj_toml):
        rep.fail(
            f"no Project.toml in {project}",
            "Point the studio at the MeanFieldHomogenization checkout with "
            "`--project <path>`.",
        )
        return rep

    manifest = os.path.join(project, "Manifest.toml")
    if not os.path.isfile(manifest):
        # Nothing to verify: Pkg will resolve from the registries, and the
        # sidecar instantiates on first start.
        return rep

    for name, rel in _path_deps(manifest):
        if rel in (".", ""):
            continue  # the project itself
        full = os.path.normpath(os.path.join(project, rel))
        if os.path.isdir(full):
            continue
        rep.missing_paths.append((name, full))
        rep.fail(
            f"`{name}` is recorded in Manifest.toml as a local checkout at "
            f"{rel}, and there is nothing at {full}",
        )

    if rep.missing_paths:
        names = ", ".join(n for n, _ in rep.missing_paths)
        parent = os.path.dirname(project)
        # These entries are development overrides, not requirements: the
        # packages they point at are published, so the quickest way out is to
        # stop using this manifest rather than to reproduce its layout.
        adds = "; ".join(f'Pkg.add("{n}")' for n, _ in rep.missing_paths)
        rep.hints.append(
            "These are development overrides pinning *published* packages to "
            "local checkouts. If you are developing MeanFieldHomogenization itself, "
            "the fix is a separate environment that takes the siblings from the "
            "registry and develops only MeanFieldHomogenization by path:\n"
            '      julia -e \'using Pkg; Pkg.activate("mfhstudio", shared=true); '
            f'{adds}; Pkg.develop(path=raw"{project}")\'\n'
            "    then start the studio with `--project @mfhstudio`.\n"
            "    (Adding them before the develop is what stops the resolver "
            "from turning them back into path entries.)"
        )
        rep.hints.append(
            f"Alternatively, clone {names} next to the MeanFieldHomogenization checkout "
            f"(inside {parent}), which is the layout the committed manifest "
            f"assumes."
        )
    return rep
