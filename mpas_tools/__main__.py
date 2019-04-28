"""
MPAS mesh tools
"""

from __future__ import absolute_import, division, print_function, \
    unicode_literals

import mpas_tools

import argparse


def main():
    """
    Entry point for the main script ``mpas_tools``
    """

    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-v', '--version',
                        action='version',
                        version='mpas_tools {}'.format(
                                mpas_tools.__version__),
                        help="Show version number and exit")

    args = parser.parse_args()


if __name__ == "__main__":
    main()

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
