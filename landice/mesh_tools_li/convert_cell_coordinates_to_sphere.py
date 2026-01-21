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
parser.add_argument("-w", dest="warp_surface", action='store_true',
                    help="warp surface by upperSurface variable. " \
                         "Uses scale factor from -s" \
                         "If upperSurface is not in file, will attempt " \
                         "to calculate it from thickness and bedTopography. " \
                         "Uses geometry from first time level only; it is not " \
                         "possible to have scaling evolve in time with this method.")
parser.add_argument("-s", dest="scale", type=float, default=100.,
                    help="scale factor for warping surface")
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

if args.warp_surface:
    if 'upperSurface' in ds.variables:
        sfc = ds.variables['upperSurface'][0,:]
    else:
        thickness = ds.variables['thickness'][0,:]
        bedTopography = ds.variables['bedTopography'][0,:]
        sfc = np.maximum(bedTopography + thickness,
                         (1.0 - 910.0 / 1028.0) * thickness)
        #sfc = thickness + np.minimum(bedTopography, 0.0) * 1028.0 / 910.0
    radius = args.radius + np.maximum(sfc, 0.0) * args.scale
else:
    radius = args.radius
xCell, yCell, zCell = compute_xyz(latCell, lonCell, radius)

# Add to NetCDF file (assumes variables already exist)
ds.variables["xCell"][:] = xCell
ds.variables["yCell"][:] = yCell
ds.variables["zCell"][:] = zCell
# Update global attributes
ds.setncattr("on_a_sphere", "YES")
ds.setncattr("sphere_radius", args.radius)
ds.close()
