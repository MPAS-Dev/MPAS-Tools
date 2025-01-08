#!/usr/bin/env python
"""
Script for adding a field named cullMask to an MPAS land ice grid for use with
the MpasCellCuller tool that actually culls the unwanted cells.
"""

from mpas_tools.landice.cull import define_cull_mask


if __name__ == '__main__':
    define_cull_mask()
