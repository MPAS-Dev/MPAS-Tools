#!/usr/bin/env python
from __future__ import print_function, division

import os
import numpy as np
import subprocess
from optparse import OptionParser
from collections import defaultdict
from netCDF4 import Dataset as NetCDFFile

parser = OptionParser()
parser.add_option("-f", "--file", dest="filename", help="Path to grid file", metavar="FILE")
parser.add_option("-m", "--metis", dest="metis_path", help="Path or name of metis executable", metavar="METIS")
parser.add_option("-p", "--procs", dest="num_procs", help="Number of processors for decomposition", metavar="PROCS")
parser.add_option("-b", "--blocks", dest="num_blocks", help="Number of blocks for decomposition", metavar="BLOCKS")
parser.add_option("-w", "--weights", dest="weight_field", help="Field to weight block partition file on.", metavar="VAR")

options, args = parser.parse_args()

if not options.metis_path:
	parser.error("Path to metis is required.")

if not options.filename:
	parser.error("A grid file is required.")

if not options.num_procs:
	parser.error("Number of processors is required.")

if not options.num_blocks:
	parser.error("Number of blocks is required.")

if not options.weight_field:
	print("Weight field missing. Defaulting to unweighted graphs.")
	weighted_parts = False
else:
	weighted_parts = True

if int(options.num_procs) > 1:
	proc_decomp = True
else:
	proc_decomp = False

dev_null = open(os.devnull, 'w')
grid = NetCDFFile(options.filename, 'r')

nCells = len(grid.dimensions['nCells'])
nEdges = len(grid.dimensions['nEdges'])

nEdgesOnCell = grid.variables['nEdgesOnCell'][:]
cellsOnCell = grid.variables['cellsOnCell'][:] - 1
if weighted_parts:
	if options.weight_field not in grid.variables:
		raise ValueError('Weight field {} not found in file.'.format(
			options.weight_field))
	weights = grid.variables[options.weight_field][:]
else:
	weights = None

grid.close()

num_blocks = 0
owning_block = [0] * nCells
block_owner = [0] * int(options.num_blocks)

nEdges = 0
for i in np.arange(0, nCells):
	for j in np.arange(0,nEdgesOnCell[i]):
		if cellsOnCell[i][j] != -1:
			nEdges = nEdges + 1

nEdges = nEdges//2

graph = open('graph.info', 'w+')
graph.write('%s %s\n'%(nCells, nEdges))
if weighted_parts:
	wgraph = open('weighted.graph.info', 'w+')
	wgraph.write('%s %s 010\n'%(nCells, nEdges))
else:
	wgraph = None

for i in np.arange(0, nCells):
	if weighted_parts:
		wgraph.write('%s '%int(weights[i]))

	for j in np.arange(0,nEdgesOnCell[i]):
		if weighted_parts:
			if(cellsOnCell[i][j] >= 0):
				wgraph.write('%s '%(cellsOnCell[i][j]+1))

		if(cellsOnCell[i][j] >= 0):
			graph.write('%s '%(cellsOnCell[i][j]+1))
	graph.write('\n')

	if weighted_parts:
		wgraph.write('\n')
graph.close()

if weighted_parts:
	wgraph.close()

command = "%s"%(options.metis_path)
if weighted_parts:
	arg1 = "weighted.graph.info"
else:
	arg1 = "graph.info"
arg2 = "%s"%options.num_blocks
subprocess.call([command, arg1, arg2], stdout=dev_null, stderr=dev_null)

if weighted_parts:
	graph = open('weighted.graph.info.part.%s'%options.num_blocks, 'r')
else:
	graph = open('graph.info.part.%s'%options.num_blocks, 'r')
i = -1
for block in iter(lambda: graph.readline(), ""):
	if i >= 0:
		block_arr = block.split()
		owning_block[i] = int(block_arr[0])+1
		num_blocks = max(num_blocks, owning_block[i])

	i = i + 1
graph.close()

if proc_decomp:
	blocksOnBlock = defaultdict(list)
	nEdges = 0

	for i in np.arange(0, nCells):
		for j in np.arange(0, nEdgesOnCell[i]):
			iCell = cellsOnCell[i][j]
			try:
				can_add = True
				for block in blocksOnBlock[owning_block[i]]:
					if block == owning_block[iCell]:
						can_add = False
			except:
				can_add = True

			if iCell == -1:
				can_add = False

			if owning_block[iCell] == owning_block[i]:
				can_add = False

			if owning_block[iCell] <= 0:
				can_add = False

			if owning_block[i] <= 0:
				can_add = False

			if can_add:
				nEdges = nEdges + 1
				blocksOnBlock[owning_block[i]].append(owning_block[iCell])

	del blocksOnBlock[0]

	block_graph = open('block.graph.info', 'w+')
	block_graph.write('%s %s\n'%(int(num_blocks), int(nEdges//2)))
	for i in np.arange(1, num_blocks+1):
		for block in blocksOnBlock[i]:
			block_graph.write('%s '%int(block))
		block_graph.write('\n')

	block_graph.close()

	command = "%s"%(options.metis_path)
	arg1 = "block.graph.info"
	arg2 = "%s"%options.num_procs
	subprocess.call([command, arg1, arg2], stdout=dev_null, stderr=dev_null)

	block_graph = open('block.graph.info.part.%s'%options.num_procs, 'r')
	iBlock = 0
	for block in iter(lambda: block_graph.readline(), ""):
		block_arr = block.split()
		block_owner[iBlock] = int(block_arr[0])
		iBlock = iBlock + 1

	block_graph.close()

	block_location = open('int_ext_blocks.dat', 'w+')
	interior_blocks = 0
	exterior_blocks = 0
	for i in np.arange(1, num_blocks+1):
		owner = block_owner[i-1]
		interior = True
		for block in blocksOnBlock[i-1]:
			if block_owner[block-1] != owner:
				interior = False

		if interior:
			block_location.write('1\n')
			interior_blocks = interior_blocks + 1
		else:
			block_location.write('0\n')
			exterior_blocks = exterior_blocks + 1

	block_location.close()

	print('Interior blocks: {}'.format(interior_blocks))
	print('Exterior blocks: {}'.format(exterior_blocks))
