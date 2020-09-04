from mpas_tools.cime.constants import constants
import requests


def test_cime_constants(cime_tag='master'):
    """
    Parse relevant constants from CIME

    Parameters
    ----------
    cime_tag : str, optional
        The CIME tag to download constants from
    """

    resp = requests.get(
        'https://raw.githubusercontent.com/ESMCI/cime/{}/src/share'
        '/util/shr_const_mod.F90'.format(cime_tag))

    text = resp.text

    text = text.split('\n')

    found = {}
    for constant in constants:
        found[constant] = False

    for line in text:
        for constant in constants:
            if constant in line:
                print('verifying {}'.format(constant))
                value = _parse_value(line, constant)
                assert value == constants[constant]
                found[constant] = True

    allFound = True
    for constant in found:
        if not found[constant]:
            print('{} was not found!'.format(constant))
            allFound = False

    assert allFound


def _parse_value(line, key):
    line, _ = line.split('!', 1)
    _, line = line.split('=')
    if '&' in line:
        raise ValueError('This parser is too dumb to handle multi-line Fortran')

    line, _ = line.split('_R8')

    try:
        value = float(line)
    except ValueError:
        value = line

    return value

if __name__ == '__main__':
    test_cime_constants()
