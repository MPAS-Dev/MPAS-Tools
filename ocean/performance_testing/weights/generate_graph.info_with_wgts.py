#!/usr/bin/env python

"""
Name: generate_graph.info_with_wgts.py
Author: Divya Jaganathan
Date: 17 July, 2018

Assigns vertex weight to each horizontal cell in graph.info (in gpmetis format)
Reads: <data_file>, <unweighed_graph_file>
Writes: graph.info_with_wgts_<X>

Flags(s) in call-command:
 -x <floating-point value> or --vertex_weight=<floating-point value>, default=0.0
 -d <data_filename> or --data_file=<data_filename>, default=init.nc
 -g <graph_file> or --graph_file=<graph_filename>, default=graph.info

"""

import numpy as np
import netCDF4 as nc4
from netCDF4 import MFDataset
import argparse

parser = \
    argparse.ArgumentParser(description=__doc__,
                            formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument(
    "-x",
    "--vertex_weight",
    dest="vertex_weight",
    help="Exponent factor in the weighing function defining dependence on depth (maxLevelCell)",
    default=0.0)

parser.add_argument(
    "-d",
    "--data_file",
    dest="data_filename",
    help="File containing the maxLevelCell data (Default: init.nc)",
    default="init.nc")

parser.add_argument(
    "-g",
    "--graph_file",
    dest="graph_filename",
    help="Unweighed graph file (Default: graph.info)",
    default="graph.info")


args = parser.parse_args()

depth_dependence_factor_x = float(args.vertex_weight)
graph_filename = args.graph_filename
data_filename = args.data_filename

file = MFDataset(data_filename)

levels = file.variables['maxLevelCell'][:]

minimum = np.amin(levels)

ratio = np.divide(levels, minimum)
weights = np.ceil((np.float_power(ratio, depth_dependence_factor_x)))
weights = weights.astype(int)
file.close()

filename = "graph.info_with_wgts_" + str(depth_dependence_factor_x)
fr = open(graph_filename, 'r')
fw = open(filename, 'w')

counter = -1

for line in fr:
    if counter == -1:
        temp = line.split("\n", 1)[0]
        # 010 indicates that the graph.info file is formatted to include the
        # cell weights
        fw.write("%s 010 \n" % temp)
    else:
        temp = line.split("\n", 1)[0]
        fw.write("%d %s \n" % (weights[counter], temp))
    counter = counter + 1

fr.close()
fw.close()
