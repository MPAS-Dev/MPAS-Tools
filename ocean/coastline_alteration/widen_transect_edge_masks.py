#!/usr/bin/env python
"""
Alter transects to be at least two cells wide.  This is used for critical
passages, to avoid sea ice blockage.  Specifically, mark cells on both sides
of each transect edge mask as a water cell.
"""
from mpas_tools.ocean.coastline_alteration import (
    main_widen_transect_edge_masks
)

if __name__ == '__main__':
    main_widen_transect_edge_masks()
