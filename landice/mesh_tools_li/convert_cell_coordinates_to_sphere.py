#!/usr/bin/env python
'''
Script to allow plotting MALI output on a sphere in Paraview.
This works by replacing x/y/zCell fields and modifying global attributes so
Paraview will plot as a spherical mesh.  Paraview only uses those 3 fields for
plotting.
The script writes to a separate output file and leaves the input file unchanged.
Note that the modified coordinates will make the file no longer compatible
with MALI or many postprocessing tools, so this is meant as a modification
for visualization in Paraview only.
'''

import argparse
import shutil
from pathlib import Path
from netCDF4 import Dataset
import numpy as np
import sys

parser = argparse.ArgumentParser(
    prog='convert_cell_coordinates_to_sphere.py',
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-i", "--infile", dest="infile", default="output.nc",
                    help="input file to convert to spherical mesh convention")
parser.add_argument("-o", "--outfile", dest="outfile", default=None,
                    help="output file to write (default: <infile stem>_sphere<suffix>)")
parser.add_argument("-r", dest="radius", type=float, default=6371220.,
                    help="sphere radius to use (m)")
parser.add_argument("-w", dest="warp_surface", action='store_true',
                    help="warp surface by upperSurface variable. " \
                         "Uses scale factor from -s. " \
                         "If upperSurface is not in file, will attempt " \
                         "to calculate it from thickness and bedTopography. " \
                         "Uses geometry from first time level only; it is not " \
                         "possible to have scaling evolve in time with this method.")
parser.add_argument("-s", dest="scale", type=float, default=100.,
                    help="scale factor for warping surface")
args = parser.parse_args()

infile = Path(args.infile)
if args.outfile is None:
    outfile = infile.with_name(f'{infile.stem}_sphere{infile.suffix}')
else:
    outfile = Path(args.outfile)

if infile.resolve() == outfile.resolve():
    raise ValueError('Output file must be different from input file.')

shutil.copyfile(infile, outfile)


def compute_xyz(latCell_deg, lonCell_deg, radius):
    '''
    Copied from mpas_tools.mesh.creation.util.lonlat2xyz in
    conda_package/mpas_tools/mesh/creation/util.py
    '''
    lat_rad = np.radians(latCell_deg)
    lon_rad = np.radians(lonCell_deg)

    xCell = radius * np.cos(lat_rad) * np.cos(lon_rad)
    yCell = radius * np.cos(lat_rad) * np.sin(lon_rad)
    zCell = radius * np.sin(lat_rad)

    return xCell, yCell, zCell


with Dataset(outfile, "r+") as ds:
    latCell = np.degrees(ds.variables["latCell"][:])
    lonCell = np.degrees(ds.variables["lonCell"][:])

    if args.warp_surface:
        has_upper_surface = 'upperSurface' in ds.variables
        has_thickness = 'thickness' in ds.variables
        has_bed_topography = 'bedTopography' in ds.variables

        if not (has_upper_surface or (has_thickness and has_bed_topography)):
            raise ValueError(
                "Cannot warp surface: provide 'upperSurface' or both "
                "'thickness' and 'bedTopography'.")

        if has_upper_surface:
            sfc = ds.variables['upperSurface'][0,:]
        else:
            thickness = ds.variables['thickness'][0,:]
            bedTopography = ds.variables['bedTopography'][0,:]
            sfc = np.maximum(bedTopography + thickness,
                             (1.0 - 910.0 / 1028.0) * thickness)
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

    # Record this modification in the global history attribute
    prev_history = getattr(ds, "history", "")
    cmd = Path(__file__).name + " " + " ".join(sys.argv[1:])
    cmd = cmd.strip()
    if prev_history:
        new_history = prev_history + "\n" + cmd
    else:
        new_history = cmd
    ds.setncattr("history", new_history)
