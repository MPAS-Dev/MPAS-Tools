"""
Tests that the dependency specs in the various files that describe the
``mpas_tools`` package stay consistent with one another.

``recipe/recipe.yaml`` is the source of truth: every package in its
``requirements: run:`` section must show up with the same version constraints
in ``pyproject.toml`` (for those available on PyPI), ``dev-spec.txt`` and
``pixi.toml``.  The latter two are also allowed to list extra packages needed
for development, testing and building the documentation.
"""

import re
import tomllib
from pathlib import Path

import pytest
import yaml

CONDA_PACKAGE_DIR = Path(__file__).resolve().parents[1]
RECIPE = CONDA_PACKAGE_DIR / 'recipe' / 'recipe.yaml'
PYPROJECT = CONDA_PACKAGE_DIR / 'pyproject.toml'
DEV_SPEC = CONDA_PACKAGE_DIR / 'dev-spec.txt'
PIXI = CONDA_PACKAGE_DIR / 'pixi.toml'

# conda-forge packages in the recipe's run requirements that are not available
# on PyPI, so they cannot be listed in pyproject.toml
CONDA_ONLY = frozenset(
    {
        'geometric-features',
        'nco',
    }
)

# conda-forge package names that differ from their PyPI equivalents
CONDA_TO_PYPI = {
    'dask-core': 'dask',
    'matplotlib-base': 'matplotlib',
    'python-igraph': 'igraph',
}

PYPI_TO_CONDA = {value: key for key, value in CONDA_TO_PYPI.items()}

# python is expressed as "requires-python" in pyproject.toml rather than as a
# dependency, so it is checked separately from the other packages
PYTHON = 'python'

# skip these tests when the spec files are not available, e.g. when the tests
# are run from the conda package rather than from a source checkout
pytestmark = pytest.mark.skipif(
    not all(path.exists() for path in (RECIPE, PYPROJECT, DEV_SPEC, PIXI)),
    reason='dependency specs are only available in a source checkout',
)


def _normalize_name(name):
    """Normalize a package name following PEP 503"""
    return re.sub(r'[-_.]+', '-', name.strip().lower())


def _normalize_constraints(constraints):
    """
    Convert a version constraint such as ``>=2.0,<3.0`` into a set of
    individual constraints so that the order and the spacing do not matter
    """
    constraints = constraints.strip()
    if constraints in ('', '*'):
        return frozenset()
    parts = constraints.split(',')
    return frozenset(part.replace(' ', '') for part in parts if part.strip())


def _parse_spec(spec):
    """Split a requirement into a normalized name and its constraints"""
    match = re.match(r'^\s*([A-Za-z0-9_.\-]+)\s*(.*)$', spec)
    if match is None:
        raise ValueError(f'Could not parse the requirement "{spec}"')
    return _normalize_name(match.group(1)), _normalize_constraints(
        match.group(2)
    )


def _parse_specs(specs):
    """Convert a list of requirements into a dict of names to constraints"""
    return dict(_parse_spec(spec) for spec in specs)


def _recipe_run_requirements():
    """The run requirements from the conda recipe"""
    with open(RECIPE) as recipe_file:
        recipe = yaml.safe_load(recipe_file)
    return _parse_specs(recipe['requirements']['run'])


def _pyproject_dependencies():
    """The (non-optional) dependencies from pyproject.toml"""
    with open(PYPROJECT, 'rb') as pyproject_file:
        pyproject = tomllib.load(pyproject_file)
    return _parse_specs(pyproject['project']['dependencies'])


def _pyproject_requires_python():
    """The python versions supported according to pyproject.toml"""
    with open(PYPROJECT, 'rb') as pyproject_file:
        pyproject = tomllib.load(pyproject_file)
    return _normalize_constraints(pyproject['project']['requires-python'])


def _dev_spec_dependencies():
    """The dependencies from dev-spec.txt"""
    specs = []
    for line in DEV_SPEC.read_text().splitlines():
        line = line.split('#')[0].strip()
        if line:
            specs.append(line)
    return _parse_specs(specs)


def _pixi_dependencies():
    """The dependencies from the pixi environment"""
    with open(PIXI, 'rb') as pixi_file:
        pixi = tomllib.load(pixi_file)
    return {
        _normalize_name(name): _normalize_constraints(constraints)
        for name, constraints in pixi['dependencies'].items()
    }


def _compare_with_recipe(dependencies, filename, rename=None, skip=()):
    """
    Check that ``dependencies`` includes each run requirement of the conda
    recipe with the same version constraints, returning a list of problems
    """
    rename = {} if rename is None else rename
    problems = []
    for name, constraints in _recipe_run_requirements().items():
        if name == PYTHON or name in skip:
            continue
        name = rename.get(name, name)
        if name not in dependencies:
            problems.append(
                f'{filename} is missing "{name}", a run requirement of '
                f'recipe/recipe.yaml'
            )
        elif dependencies[name] != constraints:
            problems.append(
                f'{filename} constrains "{name}" to '
                f'"{",".join(sorted(dependencies[name]))}" but '
                f'recipe/recipe.yaml uses "{",".join(sorted(constraints))}"'
            )
    return problems


def test_pyproject_matches_recipe():
    """pyproject.toml matches the run requirements available on PyPI"""
    problems = _compare_with_recipe(
        _pyproject_dependencies(),
        'pyproject.toml',
        rename=CONDA_TO_PYPI,
        skip=CONDA_ONLY,
    )
    assert not problems, '\n'.join(problems)


def test_pyproject_has_no_extra_dependencies():
    """pyproject.toml does not list anything the conda recipe is missing"""
    run_requirements = _recipe_run_requirements()
    problems = []
    for name in _pyproject_dependencies():
        conda_name = PYPI_TO_CONDA.get(name, name)
        if conda_name not in run_requirements:
            problems.append(
                f'pyproject.toml lists "{name}", which is not a run '
                f'requirement of recipe/recipe.yaml'
            )
    assert not problems, '\n'.join(problems)


def test_dev_spec_matches_recipe():
    """dev-spec.txt includes all the run requirements of the conda recipe"""
    problems = _compare_with_recipe(_dev_spec_dependencies(), 'dev-spec.txt')
    assert not problems, '\n'.join(problems)


def test_pixi_matches_recipe():
    """pixi.toml includes all the run requirements of the conda recipe"""
    problems = _compare_with_recipe(_pixi_dependencies(), 'pixi.toml')
    assert not problems, '\n'.join(problems)


def test_python_version_matches_pyproject():
    """
    The python versions in the other specs match "requires-python" from
    pyproject.toml.  The conda recipe is allowed to leave python
    unconstrained, since conda pins it to the version being built.
    """
    requires_python = _pyproject_requires_python()
    others = {
        'recipe/recipe.yaml': _recipe_run_requirements(),
        'dev-spec.txt': _dev_spec_dependencies(),
        'pixi.toml': _pixi_dependencies(),
    }
    problems = []
    for filename, dependencies in others.items():
        constraints = dependencies.get(PYTHON)
        if constraints is None:
            problems.append(f'{filename} does not require python')
        elif constraints and constraints != requires_python:
            problems.append(
                f'{filename} constrains python to '
                f'"{",".join(sorted(constraints))}" but the '
                f'"requires-python" in pyproject.toml is '
                f'"{",".join(sorted(requires_python))}"'
            )
    assert not problems, '\n'.join(problems)
