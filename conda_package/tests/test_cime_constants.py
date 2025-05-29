import requests

from mpas_tools.cime.constants import constants


def test_cime_constants(e3sm_tag='master'):
    """
    Parse relevant constants from CIME

    Parameters
    ----------
    e3sm_tag : str, optional
        The E3SM tag to download constants from
    """

    resp = requests.get(
        f'https://raw.githubusercontent.com/E3SM-Project/E3SM/{e3sm_tag}/'
        f'share/util/shr_const_mod.F90'
    )

    text = resp.text

    text = text.split('\n')

    found = {}
    for constant in constants:
        found[constant] = False

    for line in text:
        constant, value = _parse_value(line)
        if constant is None:
            continue
        print(f'line: {line}')
        print(f'parsed: {constant} = {value}')
        if constant in constants:
            if isinstance(value, float):
                print(f'verifying {constant}')
                assert value == constants[constant]
            else:
                print(f'skipping verification for {constant}')

            found[constant] = True
        else:
            print('not in constants')

        print('')

    all_found = True
    for constant in found:
        if not found[constant]:
            print(f'{constant} was not found!')
            all_found = False

    assert all_found


def _parse_value(line):
    if '::' not in line or '=' not in line:
        return None, None

    start = line.find('::') + 2
    end = line.find('=')

    key = line[start:end]
    line = line[end + 1 :]

    if '!' in line:
        line, _ = line.split('!', 1)

    line = line.replace('_R8', '').replace('_r8', '')

    try:
        value = float(line)
    except ValueError:
        value = line.strip()

    return key.strip(), value


if __name__ == '__main__':
    test_cime_constants()
