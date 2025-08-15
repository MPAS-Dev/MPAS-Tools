#! /usr/bin/env python
import argparse
import re
import sys
from datetime import date
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Bump version across files (CITATION.cff, __init__.py, meta.yaml)."
        )
    )
    parser.add_argument(
        "-v",
        "--version",
        required=True,
        help="New version in the form X.Y.Z (e.g., 1.3.0)",
    )
    return parser.parse_args()


def ensure_semver(version: str):
    m = re.fullmatch(r"(\d+)\.(\d+)\.(\d+)", version)
    if not m:
        sys.exit("Error: version must be in the form X.Y.Z (e.g., 1.3.0).")
    return tuple(int(p) for p in m.groups())


def read_text(path: Path):
    try:
        return path.read_text(encoding="utf-8")
    except FileNotFoundError:
        return None


def write_text_if_changed(path: Path, content: str) -> bool:
    old = read_text(path)
    if old is None:
        return False
    if old != content:
        path.write_text(content, encoding="utf-8")
        return True
    return False


def bump_citation(root: Path, version: str, today: str) -> bool:
    path = root / "CITATION.cff"
    content = read_text(path)
    if content is None:
        print(f"Skip (not found): {path}")
        return False

    # version: 1.2.1  -> version: 1.2.2
    content_new = re.sub(
        r"(?m)^(version:\s*)(.+)\s*$", rf"\g<1>{version}", content, count=1
    )
    # date-released: '2025-06-12' -> date-released: 'YYYY-MM-DD'
    content_new = re.sub(
        r"""(?m)^(date-released:\s*)['"]?\d{4}-\d{2}-\d{2}['"]?\s*$""",
        rf"\g<1>'{today}'",
        content_new,
        count=1,
    )
    # Ensure a final newline at EOF
    if not content_new.endswith("\n"):
        content_new += "\n"

    changed = write_text_if_changed(path, content_new)
    if changed:
        print(f"Updated: {path}")
    else:
        print(f"No changes: {path}")
    return changed


def bump_init(pkg_dir: Path, version_tuple) -> bool:
    path = pkg_dir / "mpas_tools" / "__init__.py"
    content = read_text(path)
    if content is None:
        print(f"Skip (not found): {path}")
        return False

    major, minor, patch = version_tuple
    content_new = re.sub(
        r"(?m)^__version_info__\s*=\s*\(.+\)\s*$",
        f"__version_info__ = ({major}, {minor}, {patch})",
        content,
        count=1,
    )

    changed = write_text_if_changed(path, content_new)
    if changed:
        print(f"Updated: {path}")
    else:
        print(f"No changes: {path}")
    return changed


def bump_meta(recipe_dir: Path, version: str) -> bool:
    path = recipe_dir / "meta.yaml"
    content = read_text(path)
    if content is None:
        print(f"Skip (not found): {path}")
        return False

    # {% set version = "1.2.1" %} -> {% set version = "1.2.2" %}
    content_new = re.sub(
        r'(?m)^\{\%\s*set\s+version\s*=\s*["\']([^"\']+)["\']\s*\%\}',
        f'{{% set version = "{version}" %}}',
        content,
        count=1,
    )

    changed = write_text_if_changed(path, content_new)
    if changed:
        print(f"Updated: {path}")
    else:
        print(f"No changes: {path}")
    return changed


def main():
    args = parse_args()
    version = args.version.strip()
    version_tuple = ensure_semver(version)
    today = date.today().isoformat()

    # Resolve repo layout relative to this script
    script_dir = Path(__file__).resolve().parent
    repo_root = script_dir
    pkg_dir = script_dir / "conda_package"
    recipe_dir = pkg_dir / "recipe"

    changed = []
    changed.append(bump_citation(repo_root, version, today))
    changed.append(bump_init(pkg_dir, version_tuple))
    changed.append(bump_meta(recipe_dir, version))

    if not any(changed):
        print("No files modified.")
        return 0
    print(f"Done. Bumped to {version} (date-released: {today}).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
