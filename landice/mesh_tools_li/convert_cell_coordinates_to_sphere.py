#!/usr/bin/env python
'''
Script to allow plotting MALI output on a sphere in Paraview.
This works by replacing x/y/zCell fields and modifying global attributes so
Paraview will plot as a spherical mesh.  Paraview only uses those 3 fields for
plotting.
The script modifies the file in place, so be sure to make a copy first!
If you screw that up, you can always copy the planar cell coordinate fields
from another file and set on_a_sphere back to NO.
Note that the modified coordinates will make the file no longer compatible
with MALI or many postprocessing tools, so this is meant as a modification
for visuatlization in Paraview only.
'''

import argparse
from netCDF4 import Dataset
import numpy as np

parser = argparse.ArgumentParser(
    prog='convert_cell_coordinates_to_sphere.py',
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-f", dest="filename", default="output.nc",
                    help="file to convert to a spherical mesh convention")
parser.add_argument("-r", dest="radius", type=float, default=6371220.,
                    help="sphere radius to use (m)")
args = parser.parse_args()


def compute_xyz(latCell_deg, lonCell_deg, radius):
    '''
    Copied from lonlat2xyz in mesh/creation/util.py
    '''
    lat_rad = np.radians(latCell_deg)
    lon_rad = np.radians(lonCell_deg)

    xCell = radius * np.cos(lat_rad) * np.cos(lon_rad)
    yCell = radius * np.cos(lat_rad) * np.sin(lon_rad)
    zCell = radius * np.sin(lat_rad)

    return xCell, yCell, zCell


ds = Dataset(args.filename, "r+")
latCell = np.degrees(ds.variables["latCell"][:])
lonCell = np.degrees(ds.variables["lonCell"][:])

xCell, yCell, zCell = compute_xyz(latCell, lonCell, args.radius)

# Add to NetCDF file (assumes variables already exist)
ds.variables["xCell"][:] = xCell
ds.variables["yCell"][:] = yCell
ds.variables["zCell"][:] = zCell
# Update global attributes
ds.setncattr("on_a_sphere", "YES")
ds.setncattr("sphere_radius", args.radius)
ds.close()
