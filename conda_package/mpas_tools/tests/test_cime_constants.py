from mpas_tools.cime.constants import constants
import requests


def test_cime_constants(e3sm_tag='master'):
    """
    Parse relevant constants from CIME

    Parameters
    ----------
    e3sm_tag : str, optional
        The E3SM tag to download CIME constants from
    """

    resp = requests.get(
        'https://raw.githubusercontent.com/E3SM-Project/E3SM/{}/cime/src/share'
        '/util/shr_const_mod.F90'.format(e3sm_tag))

    text = resp.text

    text = text.split('\n')

    for line in text:
        for constant in constants:
            if constant in line:
                print('verifying {}'.format(constant))
                value = _parse_value(line, constant)
                assert value == constants[constant]


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
