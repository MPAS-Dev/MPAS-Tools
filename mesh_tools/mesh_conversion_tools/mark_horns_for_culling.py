#!/usr/bin/env python
'''
This script identifies "horns" on a mesh (cells with two or fewer neighbors),
and marks them for culling.  In some cores/configurations, these weakly
connected cells can be dynamically inactive, and, therefore, undesirable to
keep in a mesh.

The method used will work on both planar and spherical meshes.
It adds the new masked cell to an existing 'cullCell' field if it exists,
otherwise it creates a new field.
'''

from mpas_tools.mesh.mark_horns_for_culling import main


if __name__ == '__main__':
    main()
