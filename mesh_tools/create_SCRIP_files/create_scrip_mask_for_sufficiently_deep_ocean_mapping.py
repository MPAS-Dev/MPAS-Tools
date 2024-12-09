#!/usr/bin/env python
"""
Script for defining mask of ocean cells that are deeper than a threshold.
To be used for updating a scrip file so that it only remaps
from ocean cells deeper than that threshold.
The resulting mapping file can be used as an OCN2GLC_TF_SMAPNAME
mapping file in E3SM.
Note that config_2d_thermal_forcing_depth in the MPAS-Ocean namelist
needs to match the depth used here!
"""

import argparse
import sys
from datetime import datetime
from netCDF4 import Dataset

parser = argparse.ArgumentParser(
                    description=__doc__)
parser.add_argument("-s", dest="scrip",
                    required=True,
                    help="scrip file to which to add mask")
parser.add_argument("-m", dest="mpas",
                    required=True,
                    help="MPAS-Ocean mesh file")
parser.add_argument("-d", dest="depth",
                    required=True,
                    help="depth threshold (m), should be a positive value",
                    type=float)
args = parser.parse_args()

assert args.depth > 0.0, "'depth' argument should be a positive value"

# open mesh file
fmesh = Dataset(args.mpas, 'r')
nCells = len(fmesh.dimensions['nCells'])
bottomDepth = fmesh.variables['bottomDepth'][:]

# identify cells shallower than target depth
mask = bottomDepth > args.depth

# insert mask into scrip file
fscrip = Dataset(args.scrip, 'r+')
if 'grid_imask' not in fscrip.variables:
    fscrip.createVariable('grid_imask', 'i', ('nCells',))
fscrip.variables['grid_imask'][:] = mask

# Update history attribute of scrip netCDF file
thiscommand = datetime.now().strftime("%a %b %d %H:%M:%S %Y") + ": " + " ".join(sys.argv[:])
if hasattr(fscrip, 'history'):
   newhist = '\n'.join([thiscommand, getattr(fscrip, 'history')])
else:
   newhist = thiscommand
setattr(fscrip, 'history', newhist )

fmesh.close()
fscrip.close()
