#!/usr/bin/env python
import sys, os, glob, shutil, numpy, math

from netCDF4 import *
from netCDF4 import Dataset as NetCDFFile
from pylab import *

import matplotlib
import matplotlib.pyplot as plt

from optparse import OptionParser

parser = OptionParser()
parser.add_option("-f", "--file", dest="filename", help="Path to grid file", metavar="FILE")
parser.add_option("-b", "--bins", dest="num_bins", help="Number of bins for histogram", metavar="NUMBINS")
parser.add_option("-q", "--quiet", dest="quiet", action="store_true", default=False, help="Turns off display of plots, though they are still saved in quality.png")
parser.add_option("-r", "--radius", dest="radius", help="Radius of sphere to expand for grid spacing", metavar="RAD")

options, args = parser.parse_args()

if not options.filename:
	parser.error("A grid file is required.")

if not options.num_bins:
	print "Bin number not supplied. Using default of 10."
	num_bins = 10
else:
	num_bins = int(options.num_bins)

if not options.radius:
	radius = 1.0
else:
	radius = float(options.radius)

grid = NetCDFFile(options.filename, 'a')

nCells = len(grid.dimensions['nCells'])
nEdges = len(grid.dimensions['nEdges'])
nVertices = len(grid.dimensions['nVertices'])

latCell_full = grid.variables['latCell'][:]
dcEdge_full = grid.variables['dcEdge'][:]
nEdgesOnCell_full = grid.variables['nEdgesOnCell'][:]
edgesOnCell_full = grid.variables['edgesOnCell'][:] - 1
edgesOnVertex_full = grid.variables['edgesOnVertex'][:] - 1

global_min_dc = min(dcEdge_full)
global_max_dc = max(dcEdge_full)

try:
	cellQuality_full = grid.createVariable('cellQuality', 'f8', ( 'nCells' ,) ) 
except:
	cellQuality_full = grid.variables['cellQuality']

try:
	gridSpacing_full = grid.createVariable('gridSpacing', 'f8', ( 'nCells' ,) ) 
except:
	gridSpacing_full = grid.variables['gridSpacing']

try:
	triangleQuality_full = grid.createVariable('triangleQuality', 'f8', ( 'nVertices' ,) ) 
except:
	triangleQuality_full = grid.variables['triangleQuality']

try:
	triangleAngleQuality_full = grid.createVariable('triangleAngleQuality', 'f8', ( 'nVertices' ,) ) 
except:
	triangleAngleQuality_full = grid.variables['triangleAngleQuality']

try:
	obtuseTriangle_full = grid.createVariable('obtuseTriangle', 'i4', ( 'nVertices', ) )
except:
	obtuseTriangle_full = grid.variables['obtuseTriangle']

grid_spacing_coords = np.zeros(nCells)
grid_spacing_bins = np.zeros(nCells)

bin_coords = np.arange(0.0, 1.0, 1.0/num_bins)
vor_quality_bins = np.zeros(num_bins)
del_quality_bins = np.zeros(num_bins)

edge_number_coords = np.arange(0.0, 10.0)
edge_number_bins = np.zeros(10)

vor_quality_min = 10.0
vor_quality_ave = 0.0
del_quality_min = 10.0
del_quality_ave = 0.0
del_angle_quality_min = 10.0
del_angle_quality_ave = 0.0
obtuse_triangles = 0

obtuse_angle = math.pi/2.0

for i in arange(0, nCells):
	#nEdgesOnCell = grid.variables['nEdgesOnCell'][i]
	#cellsOnCell = grid.variables['cellsOnCell'][i][:] - 1
	#edgesOnCell = grid.variables['edgesOnCell'][i][:] - 1

	nEdgesOnCell = nEdgesOnCell_full[i]
	edgesOnCell = edgesOnCell_full[i][:]

	min_dc = 10000000000.000
	max_dc = 0.0
	ave_dc = 0.0

	for j in arange(0, nEdgesOnCell):
		#dcEdge = grid.variables['dcEdge'][edgesOnCell[j]]
		dcEdge = dcEdge_full[edgesOnCell[j]]
		min_dc = min(min_dc, dcEdge)
		max_dc = max(max_dc, dcEdge)
		ave_dc = ave_dc + dcEdge

	ave_dc = ave_dc / nEdgesOnCell
	grid_spacing_coords[i] = latCell_full[i] * 180.0 / math.pi
	grid_spacing_bins[i] = ave_dc * radius
	gridSpacing_full[i] = ave_dc * radius

	sigma = min_dc / max_dc
	cellQuality_full[i] = sigma
	bin_num = min(int(sigma * num_bins), num_bins-1)
	vor_quality_bins[bin_num] = vor_quality_bins[bin_num] + 1
	vor_quality_min = min(vor_quality_min, sigma)
	vor_quality_ave = vor_quality_ave + sigma

	edge_number_bins[nEdgesOnCell-1] = edge_number_bins[nEdgesOnCell-1] + 1

for i in arange(0, nVertices):
	#edgesOnVertex = grid.variables['edgesOnVertex'][i][:] - 1
	edgesOnVertex = edgesOnVertex_full[i][:]

	#a_len = grid.variables['dcEdge'][edgesOnVertex[0]]
	#b_len = grid.variables['dcEdge'][edgesOnVertex[1]]
	#c_len = grid.variables['dcEdge'][edgesOnVertex[2]]

	if edgesOnVertex[0] < 0 or edgesOnVertex[1] < 0 or edgesOnVertex[2] < 0:
		q = 1.0
		obtuseTriangle_full[i] = 0
		angle_quality = 1.0
	else:
		a_len = dcEdge_full[edgesOnVertex[0]]
		b_len = dcEdge_full[edgesOnVertex[1]]
		c_len = dcEdge_full[edgesOnVertex[2]]

		q = (b_len + c_len - a_len) * (c_len + a_len - b_len) * (a_len + b_len - c_len) / (a_len * b_len * c_len)

		angle1 = math.acos( max(-1.0, min(1.0, (b_len**2 + c_len**2 - a_len**2) / (2 * b_len * c_len))))
		angle2 = math.acos( max(-1.0, min(1.0, (a_len**2 + c_len**2 - b_len**2) / (2 * a_len * c_len))))
		angle3 = math.acos( max(-1.0, min(1.0, (a_len**2 + b_len**2 - c_len**2) / (2 * a_len * b_len))))

		angle_quality = min(min(angle1, angle2), angle3) / max(max(angle1, angle2), angle3)

		if (angle1 > obtuse_angle) or (angle2 > obtuse_angle) or (angle3 > obtuse_angle):
			obtuse_triangles = obtuse_triangles + 1
			obtuseTriangle_full[i] = 1
		else:
			obtuseTriangle_full[i] = 0


	triangleAngleQuality_full[i] = angle_quality
	triangleQuality_full[i] = q
	bin_num = min(int(q * num_bins),num_bins-1)
	del_quality_bins[bin_num] = del_quality_bins[bin_num] + 1
	del_quality_min = min(del_quality_min, q)
	del_quality_ave = del_quality_ave + q

	if angle_quality < del_angle_quality_min:
		del_angle_quality_min = angle_quality
	del_angle_quality_ave = del_angle_quality_ave + angle_quality

vor_quality_ave = vor_quality_ave / nCells
del_quality_ave = del_quality_ave / nVertices
del_angle_quality_ave = del_angle_quality_ave / nVertices

print 'Average Cell Quality: ', vor_quality_ave
print 'Minimum Cell Quality: ', vor_quality_min

print
print 'Average Triangle Quality: ', del_quality_ave
print 'Minimum Triangle Quality: ', del_quality_min

print
print 'Average Angle Based Triangle Quality: ', del_angle_quality_ave
print 'Minimum Angle Based Triangle Quality: ', del_angle_quality_min

print
print 'Number of obtuse triangles: ', obtuse_triangles, ' out of ', nVertices

fig = plt.figure(1, figsize=(15,10))
ax = fig.add_subplot(221)

ax.set_title('Cell Quality Histogram')
ax.set_xlabel('Cell quality bins')
ax.set_ylabel('Number of cells in bin')
ax.bar(bin_coords, vor_quality_bins, 1.0/num_bins, color='r')
ax.set_ylim([0,nCells])
ax.set_xlim([0.0,1.0])

ax = fig.add_subplot(222)
ax.set_title('Triangle Quality Histogram')
ax.set_xlabel('Triangle quality bins')
ax.set_ylabel('Number of triangles in bin')
ax.bar(bin_coords, del_quality_bins, 1.0/num_bins, color='r')
ax.set_ylim([0,nVertices])
ax.set_xlim([0.0,1.0])

ax = fig.add_subplot(223)
ax.set_title('Number of cells with edge counts')
ax.set_xlabel('Number of edges')
ax.set_ylabel('Number of cells with edge count')
ax.bar(edge_number_coords, edge_number_bins, 1.0, color='r')
ax.set_ylim([0,nCells])
ax.set_xlim([0,10])

ax = fig.add_subplot(224)
ax.set_title('Grid spacing histogram')
ax.set_xlabel('Latitude (deg)')
ax.set_ylabel('Grid spacing (m)')
ax.scatter(grid_spacing_coords, grid_spacing_bins, color='r')

plt.tight_layout()
plt.draw()
plt.savefig('quality.png', bbox_inches=0)

if not options.quiet:
	plt.show()

grid.close()
