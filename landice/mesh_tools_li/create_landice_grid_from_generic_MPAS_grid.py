#!/usr/bin/env python
"""
Script to create a grid with land ice variables from an MPAS grid.
Currently variable attributes are not copied.
This script could be modified to copy them (looping over dir(var), skipping
over variable function names "assignValue", "getValue", "typecode").
"""

from mpas_tools.landice.create import create_from_generic_mpas_grid


if __name__ == '__main__':
    create_from_generic_mpas_grid()
