import textwrap

import pytest

from mpas_tools.constants import get_constant, get_pcd_version
from mpas_tools.constants.pcd import _evaluate_constants, parse_pcd


@pytest.mark.parametrize(
    'name, expected',
    [
        ('pi', 3.141592653589793238462643),
        ('speed_of_light_in_vacuum', 299792458.0),
        ('standard_acceleration_of_gravity', 9.80665),
        ('angular_velocity', 7.292115e-05),
        ('seconds_per_day', 86400.0),
        ('ocean_reference_salinity_teos10', 35.16504),
    ],
)
def test_get_constant_expected_values(name, expected):
    value = get_constant(name)
    assert isinstance(value, float)
    assert value == pytest.approx(expected)


def test_get_constant_missing_raises():
    with pytest.raises(KeyError, match='dummy_constant_that_does_not_exist'):
        get_constant('dummy_constant_that_does_not_exist')


@pytest.mark.parametrize(
    'name, expected',
    [
        ('mean_radius', (2.0 * 6378137.0 + 6356752.3141) / 3.0),
        (
            'latent_heat_of_sublimation_reference',
            333700.0 + 2501000.0,
        ),
    ],
)
def test_get_constant_formula(name, expected):
    # formulas are evaluated in double precision exactly as written, so the
    # values are bit-for-bit those that the compiled models use
    assert get_constant(name) == expected


def test_get_pcd_version():
    major, minor, patch = (int(part) for part in get_pcd_version().split('.'))
    assert (major, minor, patch) >= (0, 4, 0)


def test_parse_pcd():
    text = textwrap.dedent(
        """
        physical_constants_dictionary:
          version_number: 1.2.3
          set:
            - first:
                entries:
                  - name: a
                    value: 7292115E-11
            - second:
                entries:
                  - name: b
                    formula: 2 * a
        """
    )
    pcd_full, pcd_constants = parse_pcd(text)
    assert pcd_full['physical_constants_dictionary']['version_number'] == (
        '1.2.3'
    )
    assert pcd_constants == {'a': 7.292115e-05, 'b': 2 * 7.292115e-05}


def test_evaluate_constants_out_of_order():
    entries = {
        'half_b': {'name': 'half_b', 'formula': 'b / 2'},
        'b': {'name': 'b', 'formula': '-(a + 1) ** 2'},
        'a': {'name': 'a', 'value': 3},
    }
    constants = _evaluate_constants(entries)
    assert constants == {'a': 3.0, 'b': -16.0, 'half_b': -8.0}


@pytest.mark.parametrize(
    'entries, match',
    [
        (
            {
                'a': {'name': 'a', 'formula': 'b'},
                'b': {'name': 'b', 'formula': '2 * a'},
            },
            'Circular dependency between PCD formulas: a -> b -> a',
        ),
        (
            {'a': {'name': 'a', 'formula': 'b + 1'}},
            'references unknown constant b',
        ),
        (
            {
                'a': {'name': 'a', 'value': 4.0},
                'b': {'name': 'b', 'formula': 'sqrt(a)'},
            },
            'unsupported syntax: sqrt\\(a\\)',
        ),
        (
            {'a': {'name': 'a', 'formula': "'1'"}},
            'unsupported syntax',
        ),
    ],
)
def test_evaluate_constants_errors(entries, match):
    with pytest.raises(ValueError, match=match):
        _evaluate_constants(entries)
