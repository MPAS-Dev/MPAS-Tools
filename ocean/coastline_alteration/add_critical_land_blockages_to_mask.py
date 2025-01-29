#!/usr/bin/env python
"""
Add transects that identify critical regions where narrow strips of land block
ocean flow.  These are, essentially, the opposite of critical passages, which
must remain open for ocean flow.
"""

from mpas_tools.ocean.coastline_alteration import (
    main_add_critical_land_blockages
)

if __name__ == '__main__':
    main_add_critical_land_blockages()
