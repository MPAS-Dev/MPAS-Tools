#!/usr/bin/env python
"""
Name: paraview_cell_plots.py
Author: Doug Jacobsen
Date: 03/01/2016

This script is used to extract a field from a time series of NetCDF files as
VTK files for plotting in paraview.

It can extract a field across multiple files by passing in a regular expression
for the filename patter. As an example, one can run the script using:

`./paraview_cell_plots.py -v areaCell -f "hist.comp.*.nc"`

To extract a time series of areaCell that spans multiple files.

Requirements:
This script requires access to the following non standard modules:
vtk
netCDF4
progressbar
numpy
"""
import vtk
 
import sys, os, glob
import numpy as np

from netCDF4 import *
from netCDF4 import Dataset as NetCDFFile
import argparse

from progressbar import *

def setup_time_indices(fn_pattern, local_indices, global_indices, file_names):#{{{
    # Build file list and time indices
    file_list = sorted(glob.glob(fn_pattern))
    #local_indices = []
    #file_names = []
    #global_indices = []

    widgets = ['Build time indices: ', Percentage(), ' ', Bar(), ' ', ETA()]
    time_bar = ProgressBar(widgets=widgets, maxval=len(file_list)).start()
    i_global = 0
    i_file = 0
    for file in file_list:
        nc_file = NetCDFFile(file, 'r')

        times = len(nc_file.dimensions['Time'])

        for time_idx in np.arange(0, times):
            local_indices.append(time_idx)
            file_names.append(file)
            global_indices.append(i_global)

            i_global = i_global + 1

        i_file = i_file + 1
        nc_file.close()
        time_bar.update(i_file)
    time_bar.finish()
#}}}

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-f", "--file_pattern", dest="filename_pattern", help="MPAS Filename pattern.", metavar="FILE", required=True)
    parser.add_argument("-b", "--blocking", dest="blocking", help="Size of blocks when reading MPAS file", metavar="BLK")
    parser.add_argument("-v", "--variable", dest="variable", help="Variable to plot.", metavar="VAR", required=True)

    args = parser.parse_args()

    if not args.blocking:
        args.blocking = int(10000)
    else:
        args.blocking = int(args.blocking)

    Points = vtk.vtkPoints()
    Cells = vtk.vtkCellArray()

    time_indices = []
    time_file = []
    time_global_indices = []
    setup_time_indices(args.filename_pattern, time_indices, time_global_indices, time_file)

    nc_file = NetCDFFile(time_file[0], 'r')

    # Setting dimension values:
    first_dim_set = True
    field_dims = nc_file.variables[args.variable].dimensions
    excluded_dims = {'Time', 'nCells', 'nEdges', 'nVertices'}
    dim_vals = []
    for dim in field_dims:
        if dim == 'Time':
            dim_vals.append(-1)

        if dim == 'nCells' or dim == 'nEdges' or dim == 'nVertices':
            dim_vals.append(-2)

        if not dim in excluded_dims:
            if first_dim_set:
                print ""
                print "Need to define addition dimension values."
                first_dim_set = False

            valid = False
            while not valid:
                print "Valid range for dimension %s between 0 and %d"%(dim, len(nc_file.dimensions[dim])-1)
                val = raw_input("Enter a value for dimension %s: "%(dim))
                int_val = int(val)
                if int_val >= 0 and int_val <= len(nc_file.dimensions[dim])-1:
                    dim_vals.append( int(val) )
                    valid = True
                else:
                    print " -- Invalid value, please re-enter --"

    nCells = len(nc_file.dimensions['nCells'])
    nVertices = len(nc_file.dimensions['nVertices'])

    xVertex_var = nc_file.variables['xVertex']
    yVertex_var = nc_file.variables['yVertex']
    zVertex_var = nc_file.variables['zVertex']

    nEdgesOnCell_var = nc_file.variables['nEdgesOnCell']
    verticesOnCell_var = nc_file.variables['verticesOnCell']

    # Build vertex list
    nBlocks = 1 + nVertices / args.blocking
    widgets = ['Build vertex list: ', Percentage(), ' ', Bar(), ' ', ETA()]
    vertex_bar = ProgressBar(widgets=widgets, maxval=nBlocks).start()
    for iBlock in np.arange(0, nBlocks):
        blockStart = iBlock * args.blocking
        blockEnd = min( (iBlock + 1) * args.blocking, nVertices )
        blockCount = blockEnd - blockStart

        #print "  On vertex block %d out of %d (vertices %d through %d out of %d, contains %d vertices)"%(iBlock, nBlocks, blockStart, blockEnd, nVertices, blockCount)

        xVertex = xVertex_var[blockStart:blockEnd]
        yVertex = yVertex_var[blockStart:blockEnd]
        zVertex = zVertex_var[blockStart:blockEnd]

        for idx in np.arange(0, blockCount):
            id = Points.InsertNextPoint(xVertex[idx], yVertex[idx], zVertex[idx])
        
        del xVertex
        del yVertex
        del zVertex
        vertex_bar.update(iBlock)

    vertex_bar.finish()

    # Build cells
    nBlocks = 1 + nCells / args.blocking
    widgets = ['Build cell connectivity: ', Percentage(), ' ', Bar(), ' ', ETA()]
    cell_bar = ProgressBar(widgets=widgets, maxval=nBlocks).start()
    for iBlock in np.arange(0, nBlocks):
        blockStart = iBlock * args.blocking
        blockEnd = min( (iBlock + 1) * args.blocking, nCells )
        blockCount = blockEnd - blockStart

        nEdgesOnCell = nEdgesOnCell_var[blockStart:blockEnd]
        verticesOnCell = verticesOnCell_var[blockStart:blockEnd][:]

        for idx in np.arange(0, blockCount):
            polygon = vtk.vtkPolygon()
            polygon.GetPointIds().SetNumberOfIds(nEdgesOnCell[idx])
            for iVertex in np.arange(0, nEdgesOnCell[idx]):
                vert = verticesOnCell[idx][iVertex] - 1
                polygon.GetPointIds().SetId(iVertex, vert)

            Cells.InsertNextCell(polygon)

        del nEdgesOnCell
        del verticesOnCell
        cell_bar.update(iBlock)

    cell_bar.finish()

    nc_file.close()

    polydata = vtk.vtkPolyData()
    polydata.SetPoints(Points)
    polydata.SetPolys(Cells)

    if not os.path.exists('vtk_files/time_series'):
        os.makedirs('vtk_files/time_series')
        
    # Output time series
    widgets = ['Writing time series: ', Percentage(), ' ', Bar(), ' ', ETA()]
    field_bar = ProgressBar(widgets=widgets, maxval=len(time_global_indices)).start()
    prev_file = ""
    for time_index in time_global_indices:

        if prev_file != time_file[time_index]:
            if prev_file != "":
                nc_file.close()
            nc_file = NetCDFFile(time_file[time_index], 'r')
            prev_file = time_file[time_index]

        Colors = vtk.vtkTypeFloat64Array()
        Colors.SetNumberOfComponents(1);
        Colors.SetName(args.variable);

        field_var = nc_file.variables[args.variable]
        field_ndims = len(field_var.dimensions)

        # Build data
        for iBlock in np.arange(0, nBlocks):
            blockStart = iBlock * args.blocking
            blockEnd = min( (iBlock + 1) * args.blocking, nCells )
            blockCount = blockEnd - blockStart

            if field_ndims == 1:
                field = field_var[blockStart:blockEnd]
            if field_ndims == 2:
                field = field_var[time_indices[time_index], blockStart:blockEnd]
            elif field_ndims == 3:
                field = field_var[time_indices[time_index], blockStart:blockEnd, 0]
            elif field_ndims == 4:
                field = field_var[time_indices[time_index], blockStart:blockEnd, dim_vals[2], dim_vals[3]]
            elif field_ndims == 5:
                field = field_var[time_indices[time_index], blockStart:blockEnd, dim_vals[2], dim_vals[3], dim_vals[4]]

            for idx in np.arange(0, blockCount):
                Colors.InsertNextTuple1(field[idx]);

        polydata.GetCellData().SetScalars(Colors)
        polydata.Modified()

        out_file = "vtk_files/time_series/%s.%d.vtp"%(args.variable, time_index)

        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName( out_file )
        if vtk.VTK_MAJOR_VERSION <= 5:
            writer.SetInput(polydata)
        else:
            writer.SetInputData(polydata)
        writer.Write()

        del Colors
        del writer
        field_bar.update(time_index)

    nc_file.close()
    field_bar.finish()

    # Write pvd file, based on: http://www.paraview.org/Wiki/ParaView/Data_formats
    pvd_file = open('vtk_files/%s.pvd'%(args.variable), 'w')
    pvd_file.write('<?xml version="1.0"?>\n')
    pvd_file.write('<VTKFile type="Collection" version="0.1"\n')
    pvd_file.write('\tbyte_order="LittleEndian"\n')
    pvd_file.write('\tcompressor="vtkZLibDataCompressor">\n')
    pvd_file.write('<Collection>\n')
    for time_index in time_global_indices:
        pvd_file.write('<DataSet timestep="%d" group="" part="0"\n'%(time_index))
        pvd_file.write('\tfile="time_series/%s.%d.vtp"/>\n'%(args.variable, time_index))
    pvd_file.write('</Collection>\n')
    pvd_file.write('</VTKFile>\n')

# vim: set expandtab:
