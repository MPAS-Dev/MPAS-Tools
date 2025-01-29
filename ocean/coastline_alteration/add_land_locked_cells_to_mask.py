#!/usr/bin/env python
"""
Find ocean cells that are land-locked, and alter the cell
mask so that they are counted as land cells.
"""

from mpas_tools.ocean.coastline_alteration import (
    main_add_land_locked_cells_to_mask
)

if __name__ == '__main__':
    main_add_land_locked_cells_to_mask()
