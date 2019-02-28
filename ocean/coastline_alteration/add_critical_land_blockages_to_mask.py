#!/usr/bin/env python
"""
Name: add_critical_land_blockages_to_mask.py
Author: Xylar Asay-Davis

Add transects that identify critical regions where narrow strips of land block
ocean flow.  These are, essentially, the opposite of critical passages, which
must remain open for ocean flow.
"""

from __future__ import absolute_import, division, print_function, \
    unicode_literals

import os
import shutil
from netCDF4 import Dataset
import numpy as np
import argparse


def removeFile(fileName):
    try:
        os.remove(fileName)
    except OSError:
        pass


parser = \
    argparse.ArgumentParser(description=__doc__,
                            formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-f", "--input_mask_file", dest="input_mask_filename",
                    help="Mask file that includes cell and edge masks.",
                    metavar="INPUTMASKFILE", required=True)
parser.add_argument("-o", "--output_mask_file", dest="output_mask_filename",
                    help="Mask file that includes cell and edge masks.",
                    metavar="OUTPUTMASKFILE", required=True)
parser.add_argument("-b", "--blockage_file", dest="blockage_file",
                    help="Masks for each transect identifying critical land"
                         "blockage.", metavar="BLOCKFILE",
                    required=True)
args = parser.parse_args()

removeFile(args.output_mask_filename)
shutil.copyfile(args.input_mask_filename, args.output_mask_filename)

outMaskFile = Dataset(args.output_mask_filename, "r+")
nRegions = len(outMaskFile.dimensions["nRegions"])
regionCellMasks = outMaskFile.variables["regionCellMasks"]

blockageFile = Dataset(args.blockage_file, "r+")
nTransects = len(blockageFile.dimensions["nTransects"])
transectCellMasks = blockageFile.variables["transectCellMasks"]
for transectIndex in range(nTransects):
    # make sure the regionCellMasks for the first region is 1 anywhere a
    # transectCellMask is 1
    regionCellMasks[:, 0] = np.maximum(transectCellMasks[:, transectIndex],
                                       regionCellMasks[:, 0])

blockageFile.close()
outMaskFile.close()
