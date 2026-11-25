"""
Tests that the ``[tool.setuptools.package-data]`` patterns in
``pyproject.toml`` match every non-Python file under ``mpas_tools`` that git
tracks, so no data file is silently left out of the package.
"""

import shutil
import subprocess
import tomllib
from pathlib import Path

import pytest

CONDA_PACKAGE_DIR = Path(__file__).resolve().parents[1]
PYPROJECT = CONDA_PACKAGE_DIR / 'pyproject.toml'
PACKAGE_DIR = CONDA_PACKAGE_DIR / 'mpas_tools'


def _git_tracked_package_files():
    """
    Get the files under ``mpas_tools`` that git tracks, relative to the
    package directory, or ``None`` if this is not a git checkout
    """
    if shutil.which('git') is None or not PACKAGE_DIR.is_dir():
        return None
    result = subprocess.run(
        ['git', 'ls-files', '-z', '.'],
        cwd=PACKAGE_DIR,
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        return None
    files = {Path(name) for name in result.stdout.split('\0') if name}
    return files or None


TRACKED_FILES = _git_tracked_package_files()

# skip these tests when they are run from the conda package rather than from
# a source checkout
pytestmark = pytest.mark.skipif(
    not PYPROJECT.exists() or TRACKED_FILES is None,
    reason='package data can only be checked in a git source checkout',
)


def test_package_data_patterns_match_tracked_files():
    with open(PYPROJECT, 'rb') as f:
        pyproject = tomllib.load(f)
    patterns = pyproject['tool']['setuptools']['package-data']['mpas_tools']

    matched = set()
    for pattern in patterns:
        matched.update(
            path.relative_to(PACKAGE_DIR)
            for path in PACKAGE_DIR.glob(pattern)
            if path.is_file()
        )

    data_files = {path for path in TRACKED_FILES if path.suffix != '.py'}
    missing = sorted(str(path) for path in data_files - matched)
    assert not missing, (
        'These files tracked by git are not matched by the package-data '
        'patterns in pyproject.toml, so they would be left out of the '
        f'package: {missing}'
    )
