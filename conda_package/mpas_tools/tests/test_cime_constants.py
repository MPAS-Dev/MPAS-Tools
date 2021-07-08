from mpas_tools.cime.constants import constants
import requests


def test_cime_constants(e3sm_tag='master'):
    """
    Parse relevant constants from CIME

    Parameters
    ----------
    e3sm_tag : str, optional
        The E3SM tag to download constants from
    """

    resp = requests.get(
        'https://raw.githubusercontent.com/E3SM-Project/E3SM/{}/share/util/'
        'shr_const_mod.F90'.format(e3sm_tag))

    text = resp.text

    text = text.split('\n')

    found = {}
    for constant in constants:
        found[constant] = False

    for line in text:
        constant, value = _parse_value(line)
        if constant is None:
            continue
        print(line)
        print('parsed: {} = {}'.format(constant, value))
        if constant in constants:
            print('verifying {}'.format(constant))
            assert value == constants[constant]
            found[constant] = True

    allFound = True
    for constant in found:
        if not found[constant]:
            print('{} was not found!'.format(constant))
            allFound = False

    assert allFound


def _parse_value(line):
    if '::' not in line or '=' not in line:
        return None, None

    start = line.find('::') + 2
    end = line.find('=')

    key = line[start:end]
    line = line[end+1:]

    if '!' in line:
        line, _ = line.split('!', 1)

    if '_R8' in line:
        line, _ = line.split('_R8')

    if '_r8' in line:
        line, _ = line.split('_r8')

    try:
        value = float(line)
    except ValueError:
        value = line.strip()

    return key.strip(), value


if __name__ == '__main__':
    test_cime_constants()
