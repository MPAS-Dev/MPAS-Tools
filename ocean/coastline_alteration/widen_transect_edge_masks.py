#!/usr/bin/env python
"""
Name: widen_transect_edge_masks.py
Author: Mark Petersen

Alter transects to be at least two cells wide.  This is used for critical
passages, to avoid sea ice blockage.  Specifically, mark cells on both sides
of each transect edge mask as a water cell.
"""
import numpy as np
from netCDF4 import Dataset
import argparse

parser = \
    argparse.ArgumentParser(description=__doc__,
                            formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-f", "--mask_file", dest="mask_filename",
                    help="Mask file with cell and edge transect masks.",
                    metavar="MASKFILE",
                    required=True)
parser.add_argument("-m", "--mesh_file", dest="mesh_filename",
                    help="MPAS Mesh filename.", metavar="MESHFILE",
                    required=True)
parser.add_argument("-l", "--latitude_threshold", dest="latitude_threshold",
                    help="Minimum latitude, degrees, for transect widening.",
                    required=False, type=float, default=43.0)
args = parser.parse_args()

latitude_threshold_radians = args.latitude_threshold*3.1415/180.

# Obtain mesh variables
meshFile = Dataset(args.mesh_filename, "r")
nEdges = len(meshFile.dimensions["nEdges"])
cellsOnEdge = meshFile.variables["cellsOnEdge"][:, :]
latEdge = meshFile.variables["latEdge"][:]
meshFile.close()

# Obtain transect mask variables
maskFile = Dataset(args.mask_filename, "a")
nTransects = len(maskFile.dimensions["nTransects"])
transectCellMasks = maskFile.variables["transectCellMasks"][:, :]
transectEdgeMasks = maskFile.variables["transectEdgeMasks"][:, :]

print("widen_transect_edge_masks.py: Widening transects to two cells wide")
for iEdge in range(nEdges):
    if abs(latEdge[iEdge]) > latitude_threshold_radians:
        for iTransect in range(nTransects):
            if transectEdgeMasks[iEdge, iTransect] == 1:
                maskFile['transectCellMasks'][cellsOnEdge[iEdge, 0]-1, iTransect] = 1
                maskFile['transectCellMasks'][cellsOnEdge[iEdge, 1]-1, iTransect] = 1

maskFile.close()
