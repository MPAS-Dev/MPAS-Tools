"""
Constants from the Physical Constants Dictionary (PCD)

The PCD (https://icams-us.github.io/PhysicalConstantsDictionary/) is the
source of the physical constants in E3SM.  ``pcd.yaml`` in this package is a
copy of ``share/pcd.yaml`` from E3SM.
"""

import ast
import operator
from collections.abc import Callable
from importlib.resources import files

from ruamel.yaml import YAML

# the arithmetic that a PCD ``formula`` may use, as in the PCD's own
# ``tools/generate_source.py``
_BINARY_OPERATORS: dict[type, Callable[[float, float], float]] = {
    ast.Add: operator.add,
    ast.Sub: operator.sub,
    ast.Mult: operator.mul,
    ast.Div: operator.truediv,
    ast.Pow: operator.pow,
}

_UNARY_OPERATORS: dict[type, Callable[[float], float]] = {
    ast.UAdd: operator.pos,
    ast.USub: operator.neg,
}


def get_constant(name):
    """
    Get a constant from the Physical Constants Dictionary (PCD).

    Parameters
    ----------
    name : str
        The name of the constant to retrieve (e.g. ``'mean_radius'``).

    Returns
    -------
    value : float
        The value of the constant.

    Raises
    ------
    KeyError
        If the constant is not in the PCD.
    """
    if name not in PCD_CONSTANTS:
        raise KeyError(
            f'Constant {name} not found in the Physical Constants Dictionary'
        )
    return PCD_CONSTANTS[name]


def get_pcd_version():
    """
    Get the version of the PCD that MPAS-Tools packages.

    Returns
    -------
    version : str
        The PCD version.
    """
    return str(PCD_FULL['physical_constants_dictionary']['version_number'])


def parse_pcd(text):
    """
    Parse the constants from the text of a PCD YAML file.

    Parameters
    ----------
    text : str
        The contents of a ``pcd.yaml`` file.

    Returns
    -------
    pcd_full : dict
        The full contents of the PCD.

    pcd_constants : dict
        The value of each constant, keyed by name, with any constants defined
        by a ``formula`` evaluated in double precision.
    """
    yaml_data = YAML(typ='rt')
    pcd_full = yaml_data.load(text)

    entries = {}
    # constants are within the physical_constants_dictionary and set sections
    sets = pcd_full['physical_constants_dictionary']['set']
    # they may be within any of the available subsections
    for set_dict in sets:
        for entry in next(iter(set_dict.values()))['entries']:
            entries[entry['name']] = entry

    return pcd_full, _evaluate_constants(entries)


def _evaluate_constants(entries):
    """
    Get the value of each PCD constant, evaluating any that are defined by a
    ``formula`` rather than a ``value``

    Parameters
    ----------
    entries : dict
        PCD entries, keyed by constant name

    Returns
    -------
    pcd_constants : dict
        The value of each constant, keyed by name
    """
    pcd_constants: dict = {}
    for name in entries:
        _evaluate_constant(name, entries, pcd_constants, stack=[])
    return pcd_constants


def _evaluate_constant(name, entries, pcd_constants, stack):
    """
    Get the value of a PCD constant, first evaluating any constants its
    formula references
    """
    if name in pcd_constants:
        return pcd_constants[name]

    if name in stack:
        cycle = ' -> '.join(stack + [name])
        raise ValueError(f'Circular dependency between PCD formulas: {cycle}')

    entry = entries[name]
    if 'formula' in entry:
        formula = str(entry['formula'])
        tree = ast.parse(formula, mode='eval')
        value = _evaluate_node(
            tree.body, entries, pcd_constants, stack=stack + [name]
        )
    else:
        value = float(entry['value'])

    pcd_constants[name] = value
    return value


def _evaluate_node(node, entries, pcd_constants, stack):
    """
    Evaluate a node of a PCD formula in double precision, allowing only
    arithmetic on numbers and other PCD constants
    """
    name = stack[-1]
    if isinstance(node, ast.BinOp) and type(node.op) in _BINARY_OPERATORS:
        left = _evaluate_node(node.left, entries, pcd_constants, stack)
        right = _evaluate_node(node.right, entries, pcd_constants, stack)
        return _BINARY_OPERATORS[type(node.op)](left, right)

    if isinstance(node, ast.UnaryOp) and type(node.op) in _UNARY_OPERATORS:
        operand = _evaluate_node(node.operand, entries, pcd_constants, stack)
        return _UNARY_OPERATORS[type(node.op)](operand)

    if (
        isinstance(node, ast.Constant)
        and isinstance(node.value, (int, float))
        and not isinstance(node.value, bool)
    ):
        return float(node.value)

    if isinstance(node, ast.Name):
        if node.id not in entries:
            raise ValueError(
                f'The PCD formula for {name} references unknown constant '
                f'{node.id}.'
            )
        return _evaluate_constant(node.id, entries, pcd_constants, stack)

    raise ValueError(
        f'The PCD formula for {name} uses unsupported syntax: '
        f'{ast.unparse(node)}'
    )


# parsed full dictionary and constants for quick lookup
PCD_FULL, PCD_CONSTANTS = parse_pcd(
    files('mpas_tools.constants')
    .joinpath('pcd.yaml')
    .read_text(encoding='utf-8')
)
