import ast
import operator
import re

import requests

from mpas_tools.cime.constants import constants
from mpas_tools.constants.pcd import PCD_CONSTANTS, get_pcd_version, parse_pcd

E3SM_TAG = 'master'

_BINARY_OPERATORS = {
    ast.Add: operator.add,
    ast.Sub: operator.sub,
    ast.Mult: operator.mul,
    ast.Div: operator.truediv,
}

_UNARY_OPERATORS = {
    ast.UAdd: operator.pos,
    ast.USub: operator.neg,
}


def test_pcd_matches_e3sm(e3sm_tag=E3SM_TAG):
    """
    Check that the packaged PCD matches ``share/pcd.yaml`` in E3SM

    Parameters
    ----------
    e3sm_tag : str, optional
        The E3SM tag or branch to download the PCD from
    """
    text = _get_e3sm_file('share/pcd.yaml', e3sm_tag)
    pcd_full, pcd_constants = parse_pcd(text)

    e3sm_version = str(
        pcd_full['physical_constants_dictionary']['version_number']
    )
    assert get_pcd_version() == e3sm_version, (
        f'MPAS-Tools has PCD version {get_pcd_version()} but E3SM has '
        f'{e3sm_version}.  Please update mpas_tools/constants/pcd.yaml.'
    )
    assert PCD_CONSTANTS == pcd_constants


def test_cime_constants(e3sm_tag=E3SM_TAG):
    """
    Check that each constant in :py:mod:`mpas_tools.cime.constants` has the
    value that E3SM's ``shr_const_mod`` gives it, evaluating each expression
    with constants from the packaged PCD

    Parameters
    ----------
    e3sm_tag : str, optional
        The E3SM tag or branch to download ``shr_const_mod.F90`` from
    """
    text = _get_e3sm_file('share/util/shr_const_mod.F90', e3sm_tag)
    statements = _get_statements(text)

    names = {}
    for statement in statements:
        names.update(_parse_pcd_aliases(statement))

    assert len(names) > 0, 'Found no constants used from the pcd module'

    for statement in statements:
        constant, expression = _parse_parameter(statement)
        if constant is None:
            continue
        value = _evaluate(ast.parse(expression, mode='eval').body, names)
        print(f'{constant} = {expression} = {value}')
        names[constant] = value

    missing = [constant for constant in constants if constant not in names]
    assert not missing, f'Not found in shr_const_mod: {missing}'

    mismatched = {
        constant: (constants[constant], names[constant])
        for constant in constants
        if constants[constant] != names[constant]
    }
    assert not mismatched, (
        f'Constants that differ from E3SM (MPAS-Tools, E3SM): {mismatched}'
    )


def _get_e3sm_file(path, e3sm_tag):
    resp = requests.get(
        f'https://raw.githubusercontent.com/E3SM-Project/E3SM/{e3sm_tag}/'
        f'{path}',
        timeout=60,
    )
    resp.raise_for_status()
    return resp.text


def _get_statements(text):
    """Remove comments and join continuation lines into Fortran statements"""
    statements = []
    current = ''
    for line in text.split('\n'):
        line = line.split('!', 1)[0].rstrip()
        if line.endswith('&'):
            current = f'{current} {line[:-1].strip()}'
        else:
            statements.append(f'{current} {line.strip()}'.strip())
            current = ''
    return statements


def _parse_pcd_aliases(statement):
    """Get the PCD value of each local name in a ``use pcd`` statement"""
    match = re.match(r'use\s+pcd\s*,\s*only\s*:(.*)', statement, re.IGNORECASE)
    if match is None:
        return {}

    aliases = {}
    for item in match.group(1).split(','):
        local, pcd_name = (part.strip() for part in item.split('=>'))
        aliases[local] = PCD_CONSTANTS[pcd_name]
    return aliases


def _parse_parameter(statement):
    """Get the name and expression of a ``SHR_CONST_*`` parameter"""
    match = re.match(
        r'real\s*\(\s*R8\s*\)\s*,\s*parameter\s*::\s*(SHR_CONST_\w+)\s*=(.*)',
        statement,
        re.IGNORECASE,
    )
    if match is None:
        return None, None

    # remove the kind from Fortran literals like 1.0_R8
    expression = re.sub(r'_[rR]8\b', '', match.group(2).strip())
    return match.group(1), expression


def _evaluate(node, names):
    """Evaluate an arithmetic expression from shr_const_mod"""
    if isinstance(node, ast.BinOp) and type(node.op) in _BINARY_OPERATORS:
        left = _evaluate(node.left, names)
        right = _evaluate(node.right, names)
        return _BINARY_OPERATORS[type(node.op)](left, right)

    if isinstance(node, ast.UnaryOp) and type(node.op) in _UNARY_OPERATORS:
        return _UNARY_OPERATORS[type(node.op)](_evaluate(node.operand, names))

    if isinstance(node, ast.Constant) and isinstance(node.value, (int, float)):
        return float(node.value)

    if isinstance(node, ast.Name):
        return names[node.id]

    raise ValueError(f'Unsupported syntax: {ast.unparse(node)}')


if __name__ == '__main__':
    test_pcd_matches_e3sm()
    test_cime_constants()
