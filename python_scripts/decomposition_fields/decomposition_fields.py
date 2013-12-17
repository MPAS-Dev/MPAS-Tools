#!/usr/bin/python
import sys, os, glob, shutil, numpy, math

from netCDF4 import *
from netCDF4 import Dataset as NetCDFFile
from pylab import *

import matplotlib
import matplotlib.pyplot as plt

from optparse import OptionParser

parser = OptionParser()
parser.add_option("-f", "--file", dest="filename", help="Path to grid file", metavar="FILE")
parser.add_option("-d", "--decomp", dest="decomposition", help="Decomposition file", metavar="DECOMP")

options, args = parser.parse_args()

if not options.filename:
	parser.error("A grid file is required.")

if not options.decomposition:
	parser.error("A decomposition file is required.")

grid = NetCDFFile(options.filename, 'a')

nCells = len(grid.dimensions['nCells'])
nEdges = len(grid.dimensions['nEdges'])
nVertices = len(grid.dimensions['nVertices'])
vertexDegree = len(grid.dimensions['vertexDegree'])

cellsOnEdge_full = grid.variables['cellsOnEdge'][:] -1
cellsOnVertex_full = grid.variables['cellsOnVertex'][:] -1

try:
	cellDecomposition_full = grid.createVariable('cellDecompositon', 'i4', ( 'nCells' ,) ) 
except:
	cellDecomposition_full = grid.variables['cellDecompositon']

try:
	triangleDecomposition_full = grid.createVariable('triangleDecomposition', 'i4', ( 'nVertices' ,) ) 
except:
	triangleDecomposition_full = grid.variables['triangleDecomposition']

try:
	edgeDecomposition_full = grid.createVariable('edgeDecomposition', 'i4', ( 'nEdges', ) )
except:
	edgeDecomposition_full = grid.variables['edgeDecomposition']

print 'Reading decomposition file'
i = 0
decomp_file = open(options.decomposition, 'r')
for block in iter(lambda: decomp_file.readline(), ""):
	block_arr = block.split()
	cellDecomposition_full[i] = int(block_arr[0])
	i = i + 1

decomp_file.close()

print 'Computing vertex owner'
for i in arange(0, nVertices):
	owner = cellDecomposition_full[cellsOnVertex_full[i][0]]

	triangleDecomposition_full[i] = owner

print 'Computing edge owner'
for i in arange(0, nEdges):
	owner = cellDecomposition_full[cellsOnEdge_full[i][0]]

	edgeDecomposition_full[i] = owner

grid.close()
