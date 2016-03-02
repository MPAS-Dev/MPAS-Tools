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

def setup_dimension_values_and_sort_vars( nc_file, variable_list, all_dim_vals, cellVarTime, cellVarNoTime, vertexVarTime, vertexVarNoTime, edgeVarTime, edgeVarNoTime):#{{{
    for var in variable_list.split(','):
        captured_input = False

        dim_vals = []
        field_var = nc_file.variables[var]

        # Setting dimension values:
        first_dim_set = True
        field_dims = field_var.dimensions
        excluded_dims = {'Time', 'nCells', 'nEdges', 'nVertices'}
        for dim in field_dims:
            if dim == 'Time':
                dim_vals.append(-1)

            if dim == 'nCells' or dim == 'nEdges' or dim == 'nVertices':
                dim_vals.append(-2)

            if not dim in excluded_dims:
                captured_input = True
                if first_dim_set:
                    print ""
                    print "Need to define additional dimension values for field %s"%(var)
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

        if "Time" in field_dims:
            if "nCells" in field_dims:
                cellVarTime.append(var)
            elif "nVertices" in field_dims:
                vertexVarTime.append(var)
            elif "nEdges" in field_dims:
                edgeVarTime.append(var)
        else:
            if "nCells" in field_dims:
                cellVarNoTime.append(var)
            elif "nVertices" in field_dims:
                vertexVarNoTime.append(var)
            elif "nEdges" in field_dims:
                edgeVarNoTime.append(var)

        all_dim_vals[var] = dim_vals
        del dim_vals

        if captured_input:
            print ""
#}}}

def summarize_extraction(all_dim_vals, cellVarTime, cellVarNoTime, vertexVarTime, vertexVarNoTime, edgeVarTime, edgeVarNoTime):#{{{
    print ""
    print "The following variables will be extracted from the input file(s)."
    print ""

    if len(cellVarTime) > 0:
        print "   Variables with 'Time' and 'nCells' as dimensions:"
        for var in cellVarTime:
            print "      name: %s"%(var)

    if len(cellVarNoTime) > 0:
        print "   Variables with 'nCells' and without 'Time' as dimensions:"
        for var in cellVarNoTime:
            print "      name: %s"%(var)

    if len(vertexVarTime) > 0:
        print "   Variables with 'Time' and 'nVertices' as dimensions:"
        for var in vertexVarTime:
            print "      name: %s"%(var)

    if len(vertexVarNoTime) > 0:
        print "   Variables with 'nVertices' and without 'Time' as dimensions:"
        for var in vertexVarNoTime:
            print "      name: %s"%(var)

    if len(edgeVarTime) > 0:
        print " Edge fields are currently not supported, any fields listed are ignored."
        print "   Variables with 'Time' and 'nEdges' as dimensions:"
        for var in edgeVarTime:
            print "      name: %s"%(var)

    if len(edgeVarNoTime) > 0:
        print " Edge fields are currently not supported, any fields listed are ignored."
        print "   Variables with 'nEdges' and without 'Time' as dimensions:"
        for var in edgeVarNoTime:
            print "      name: %s"%(var)

    print ""
#}}}

def build_location_lists_xyz( nc_file, blocking, dimName, xName, yName, zName, point_list ):#{{{
    nLocs = len(nc_file.dimensions[dimName])

    xLoc_var = nc_file.variables[xName]
    yLoc_var = nc_file.variables[yName]
    zLoc_var = nc_file.variables[zName]

    # Build vertex list
    nBlocks = 1 + nLocs / blocking
    widgets = ['Build location list: ', Percentage(), ' ', Bar(), ' ', ETA()]
    loc_bar = ProgressBar(widgets=widgets, maxval=nBlocks).start()
    for iBlock in np.arange(0, nBlocks):
        blockStart = iBlock * blocking
        blockEnd = min( (iBlock + 1) * blocking, nLocs )
        blockCount = blockEnd - blockStart

        #print "  On vertex block %d out of %d (vertices %d through %d out of %d, contains %d vertices)"%(iBlock, nBlocks, blockStart, blockEnd, nVertices, blockCount)

        xLoc = xLoc_var[blockStart:blockEnd]
        yLoc = yLoc_var[blockStart:blockEnd]
        zLoc = zLoc_var[blockStart:blockEnd]

        for idx in np.arange(0, blockCount):
            id = point_list.InsertNextPoint(xLoc[idx], yLoc[idx], zLoc[idx])
        
        del xLoc
        del yLoc
        del zLoc
        loc_bar.update(iBlock)

    loc_bar.finish()
#}}}

def build_cell_lists( nc_file, blocking, cell_list ):#{{{
    nCells = len(nc_file.dimensions['nCells'])

    nEdgesOnCell_var = nc_file.variables['nEdgesOnCell']
    verticesOnCell_var = nc_file.variables['verticesOnCell']

    # Build cells
    nBlocks = 1 + nCells / blocking
    widgets = ['Build cell connectivity: ', Percentage(), ' ', Bar(), ' ', ETA()]
    cell_bar = ProgressBar(widgets=widgets, maxval=nBlocks).start()
    for iBlock in np.arange(0, nBlocks):
        blockStart = iBlock * blocking
        blockEnd = min( (iBlock + 1) * blocking, nCells )
        blockCount = blockEnd - blockStart

        nEdgesOnCell = nEdgesOnCell_var[blockStart:blockEnd]
        verticesOnCell = verticesOnCell_var[blockStart:blockEnd][:]

        for idx in np.arange(0, blockCount):
            polygon = vtk.vtkPolygon()
            polygon.GetPointIds().SetNumberOfIds(nEdgesOnCell[idx])
            for iVertex in np.arange(0, nEdgesOnCell[idx]):
                vert = verticesOnCell[idx][iVertex] - 1
                polygon.GetPointIds().SetId(iVertex, vert)

            cell_list.InsertNextCell(polygon)

        del nEdgesOnCell
        del verticesOnCell
        cell_bar.update(iBlock)

    cell_bar.finish()

#}}}

def build_dual_cell_lists( nc_file, blocking, cell_list ):#{{{
    nVertices = len(nc_file.dimensions['nVertices'])
    vertexDegree = len(nc_file.dimensions['vertexDegree'])

    cellsOnVertex_var = nc_file.variables['cellsOnVertex']

    # Build cells
    nBlocks = 1 + nVertices / blocking
    widgets = ['Build dual connectivity: ', Percentage(), ' ', Bar(), ' ', ETA()]
    dual_bar = ProgressBar(widgets=widgets, maxval=nBlocks).start()
    for iBlock in np.arange(0, nBlocks):
        blockStart = iBlock * blocking
        blockEnd = min( (iBlock + 1) * blocking, nVertices )
        blockCount = blockEnd - blockStart

        cellsOnVertex = cellsOnVertex_var[blockStart:blockEnd][:]

        keep_cell = True

        for idx in np.arange(0, blockCount):
            polygon = vtk.vtkPolygon()
            polygon.GetPointIds().SetNumberOfIds(vertexDegree)

            for iCell in np.arange(0, vertexDegree):
                cell = cellsOnVertex[idx][iCell] - 1
                polygon.GetPointIds().SetId(iVertex, cell)
                keep_cell = keep_cell and ( cell >= 0 )

            if keep_cell:
                cell_list.InsertNextCell(polygon)

        del cellsOnVertex
        dual_bar.update(iBlock)

    dual_bar.finish()
#}}}

def build_field_time_series( local_indices, global_indices, file_names, blocking, all_dim_vals, variable_list, point_list, cell_list ):#{{{
    polydata = vtk.vtkPolyData()
    polydata.SetPoints(Points)
    polydata.SetPolys(Cells)

    written_fields = []

    # Output time series
    widgets = ['Writing time series: ', Percentage(), ' ', Bar(), ' ', ETA()]
    field_bar = ProgressBar(widgets=widgets, maxval=len(global_indices) ).start()
    prev_file = ""
    for time_index in global_indices:

        if prev_file != file_names[time_index]:
            if prev_file != "":
                nc_file.close()
            nc_file = NetCDFFile(file_names[time_index], 'r')
            prev_file = file_names[time_index]

        nCells = len(nc_file.dimensions['nCells'])
        nBlocks = 1 + nCells / blocking

        for var_name in variable_list:
            Colors = vtk.vtkTypeFloat64Array()
            Colors.SetNumberOfComponents(1);
            Colors.SetName(var_name);

            dim_vals = all_dim_vals[var_name]
            field_var = nc_file.variables[var_name]
            field_ndims = len(dim_vals)
            field_dims = field_var.dimensions

            # Build data
            for iBlock in np.arange(0, nBlocks):
                blockStart = iBlock * blocking
                blockEnd = min( (iBlock + 1) * blocking, nCells )
                blockCount = blockEnd - blockStart

                if field_ndims == 1:
                    field = field_var[blockStart:blockEnd]
                if field_ndims == 2:
                    field = field_var[local_indices[time_index], blockStart:blockEnd]
                elif field_ndims == 3:
                    field = field_var[local_indices[time_index], blockStart:blockEnd, dim_vals[2]]
                elif field_ndims == 4:
                    field = field_var[local_indices[time_index], blockStart:blockEnd, dim_vals[2], dim_vals[3]]
                elif field_ndims == 5:
                    field = field_var[local_indices[time_index], blockStart:blockEnd, dim_vals[2], dim_vals[3], dim_vals[4]]

                for idx in np.arange(0, blockCount):
                    Colors.InsertNextTuple1(field[idx]);

            polydata.GetCellData().SetScalars(Colors)
            polydata.Modified()

            out_file = "vtk_files/time_series/%s.%d.vtp"%(var_name, time_index)

            writer = vtk.vtkXMLPolyDataWriter()
            writer.SetFileName( out_file )
            if vtk.VTK_MAJOR_VERSION <= 5:
                writer.SetInput(polydata)
            else:
                writer.SetInputData(polydata)
            writer.Write()

            del Colors
            del writer
            del field
            del field_ndims
            del field_var
            del dim_vals

        field_bar.update(time_index)

    nc_file.close()
    field_bar.finish()

    # Write pvd file, based on: http://www.paraview.org/Wiki/ParaView/Data_formats
    for var_name in variable_list:
        pvd_file = open('vtk_files/%s.pvd'%(var_name), 'w')
        pvd_file.write('<?xml version="1.0"?>\n')
        pvd_file.write('<VTKFile type="Collection" version="0.1"\n')
        pvd_file.write('\tbyte_order="LittleEndian"\n')
        pvd_file.write('\tcompressor="vtkZLibDataCompressor">\n')
        pvd_file.write('<Collection>\n')
        for time_index in time_global_indices:
            pvd_file.write('<DataSet timestep="%d" group="" part="0"\n'%(time_index))
            pvd_file.write('\tfile="time_series/%s.%d.vtp"/>\n'%(var_name, time_index))
        pvd_file.write('</Collection>\n')
        pvd_file.write('</VTKFile>\n')

    del written_fields
#}}}

def build_field_single_time_field( local_indices, global_indices, file_names, blocking, all_dim_vals, variable_list, point_list, cell_list ):#{{{
    polydata = vtk.vtkPolyData()
    polydata.SetPoints(Points)
    polydata.SetPolys(Cells)

    # Output time series
    widgets = ['Writing single fields: ', Percentage(), ' ', Bar(), ' ', ETA()]
    field_bar = ProgressBar(widgets=widgets, maxval=len(variable_list) ).start()
    prev_file = ""
    nc_file = NetCDFFile( file_names[0], 'r' )

    nCells = len(nc_file.dimensions['nCells'])
    nBlocks = 1 + nCells / blocking

    var_index = 0
    for var_name in variable_list:
        Colors = vtk.vtkTypeFloat64Array()
        Colors.SetNumberOfComponents(1);
        Colors.SetName(var_name);

        dim_vals = all_dim_vals[var_name]
        field_var = nc_file.variables[var_name]
        field_ndims = len(dim_vals)

        # Build data
        for iBlock in np.arange(0, nBlocks):
            blockStart = iBlock * blocking
            blockEnd = min( (iBlock + 1) * blocking, nCells )
            blockCount = blockEnd - blockStart

            if field_ndims == 1:
                field = field_var[blockStart:blockEnd]
            if field_ndims == 2:
                field = field_var[blockStart:blockEnd, dim_vals[1]]
            elif field_ndims == 3:
                field = field_var[blockStart:blockEnd, dim_vals[1]]
            elif field_ndims == 4:
                field = field_var[blockStart:blockEnd, dim_vals[1], dim_vals[2]]
            elif field_ndims == 5:
                field = field_var[blockStart:blockEnd, dim_vals[1], dim_vals[2], dim_vals[3]]

            for idx in np.arange(0, blockCount):
                Colors.InsertNextTuple1(field[idx]);

        polydata.GetCellData().SetScalars(Colors)
        polydata.Modified()

        out_file = "vtk_files/%s.vtp"%(var_name)

        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName( out_file )
        if vtk.VTK_MAJOR_VERSION <= 5:
            writer.SetInput(polydata)
        else:
            writer.SetInputData(polydata)
        writer.Write()

        del Colors
        del writer
        del field
        del field_ndims
        del field_var
        del dim_vals

        var_index = var_index + 1
        field_bar.update(var_index)

    nc_file.close()
    field_bar.finish()
#}}}

def build_location_list_xyz( nc_file, blocking, dimName, xName, yName, zName, point_list ):#{{{
    nLocs = len(nc_file.dimensions[dimName])

    xLoc_var = nc_file.variables[xName]
    yLoc_var = nc_file.variables[yName]
    zLoc_var = nc_file.variables[zName]

    # Build vertex list
    nBlocks = 1 + nLocs / blocking
    widgets = ['Build location list: ', Percentage(), ' ', Bar(), ' ', ETA()]
    loc_bar = ProgressBar(widgets=widgets, maxval=nBlocks).start()
    for iBlock in np.arange(0, nBlocks):
        blockStart = iBlock * blocking
        blockEnd = min( (iBlock + 1) * blocking, nLocs )
        blockCount = blockEnd - blockStart

        #print "  On vertex block %d out of %d (vertices %d through %d out of %d, contains %d vertices)"%(iBlock, nBlocks, blockStart, blockEnd, nVertices, blockCount)

        xLoc = xLoc_var[blockStart:blockEnd]
        yLoc = yLoc_var[blockStart:blockEnd]
        zLoc = zLoc_var[blockStart:blockEnd]

        for idx in np.arange(0, blockCount):
            id = point_list.InsertNextPoint(xLoc[idx], yLoc[idx], zLoc[idx])
        
        del xLoc
        del yLoc
        del zLoc
        loc_bar.update(iBlock)

    id = point_list.InsertNextPoint(0.0, 0.0, 0.0)

    loc_bar.finish()
#}}}

def build_cell_lists( nc_file, blocking, cell_list ):#{{{
    nCells = len(nc_file.dimensions['nCells'])

    nEdgesOnCell_var = nc_file.variables['nEdgesOnCell']
    verticesOnCell_var = nc_file.variables['verticesOnCell']

    # Build cells
    nBlocks = 1 + nCells / blocking
    widgets = ['Build cell connectivity: ', Percentage(), ' ', Bar(), ' ', ETA()]
    cell_bar = ProgressBar(widgets=widgets, maxval=nBlocks).start()
    for iBlock in np.arange(0, nBlocks):
        blockStart = iBlock * blocking
        blockEnd = min( (iBlock + 1) * blocking, nCells )
        blockCount = blockEnd - blockStart

        nEdgesOnCell = nEdgesOnCell_var[blockStart:blockEnd]
        verticesOnCell = verticesOnCell_var[blockStart:blockEnd][:]

        for idx in np.arange(0, blockCount):
            polygon = vtk.vtkPolygon()
            polygon.GetPointIds().SetNumberOfIds(nEdgesOnCell[idx])
            keep_cell = True

            for iVertex in np.arange(0, nEdgesOnCell[idx]):
                vert = verticesOnCell[idx][iVertex] - 1
                if vert == -1:
                    keep_cell = False

                polygon.GetPointIds().SetId(iVertex, vert)

            if not keep_cell:
                for iVertex in np.arange(0, nEdgesOnCell[idx]):
                    polygon.GetPointIds().SetId(iVertex, nCells+1)

            cell_list.InsertNextCell(polygon)

        del nEdgesOnCell
        del verticesOnCell
        cell_bar.update(iBlock)

    cell_bar.finish()

#}}}

def build_dual_cell_lists( nc_file, blocking, cell_list ):#{{{
    nVertices = len(nc_file.dimensions['nVertices'])
    vertexDegree = len(nc_file.dimensions['vertexDegree'])

    cellsOnVertex_var = nc_file.variables['cellsOnVertex']

    # Build cells
    nBlocks = 1 + nVertices / blocking
    widgets = ['Build dual connectivity: ', Percentage(), ' ', Bar(), ' ', ETA()]
    dual_bar = ProgressBar(widgets=widgets, maxval=nBlocks).start()
    for iBlock in np.arange(0, nBlocks):
        blockStart = iBlock * blocking
        blockEnd = min( (iBlock + 1) * blocking, nVertices )
        blockCount = blockEnd - blockStart

        cellsOnVertex = cellsOnVertex_var[blockStart:blockEnd][:]

        keep_cell = True

        for idx in np.arange(0, blockCount):
            polygon = vtk.vtkPolygon()

            polygon.GetPointIds().SetNumberOfIds(vertexDegree)

            keep_cell = True
            for iCell in np.arange(0, vertexDegree):
                cell = cellsOnVertex[idx][iCell] - 1
                if cell == -1:
                    keep_cell = False
                else:
                    polygon.GetPointIds().SetId(iCell, cell)


            if not keep_cell:
                for iCell in np.arange(0, vertexDegree):
                    polygon.GetPointIds().SetId(iCell, nVertices+1)

            cell_list.InsertNextCell(polygon)

        del cellsOnVertex
        dual_bar.update(iBlock)

    dual_bar.finish()
#}}}

def build_field_time_series( local_indices, global_indices, file_names, blocking, all_dim_vals, blockDimName, variable_list, point_list, cell_list ):#{{{
    if len(variable_list) > 0:
        polydata = vtk.vtkPolyData()
        polydata.SetPoints(point_list)
        polydata.SetPolys(cell_list)

        written_fields = []

        # Output time series
        widgets = ['Writing time series: ', Percentage(), ' ', Bar(), ' ', ETA()]
        field_bar = ProgressBar(widgets=widgets, maxval=len(global_indices) ).start()
        prev_file = ""
        for time_index in global_indices:

            if prev_file != file_names[time_index]:
                if prev_file != "":
                    nc_file.close()
                nc_file = NetCDFFile(file_names[time_index], 'r')
                prev_file = file_names[time_index]

            blockDim = len(nc_file.dimensions[blockDimName])
            nBlocks = 1 + blockDim / blocking

            for var_name in variable_list:
                Colors = vtk.vtkTypeFloat64Array()
                Colors.SetNumberOfComponents(1);
                Colors.SetName(var_name);

                dim_vals = all_dim_vals[var_name]
                field_var = nc_file.variables[var_name]
                field_ndims = len(dim_vals)
                field_dims = field_var.dimensions

                if "Time" in field_dims:
                    written_fields.append(var_name)

                    # Build data
                    for iBlock in np.arange(0, nBlocks):
                        blockStart = iBlock * blocking
                        blockEnd = min( (iBlock + 1) * blocking, blockDim )
                        blockCount = blockEnd - blockStart

                        if field_ndims == 1:
                            field = field_var[blockStart:blockEnd]
                        if field_ndims == 2:
                            field = field_var[local_indices[time_index], blockStart:blockEnd]
                        elif field_ndims == 3:
                            field = field_var[local_indices[time_index], blockStart:blockEnd, dim_vals[2]]
                        elif field_ndims == 4:
                            field = field_var[local_indices[time_index], blockStart:blockEnd, dim_vals[2], dim_vals[3]]
                        elif field_ndims == 5:
                            field = field_var[local_indices[time_index], blockStart:blockEnd, dim_vals[2], dim_vals[3], dim_vals[4]]

                        for idx in np.arange(0, blockCount):
                            Colors.InsertNextTuple1(field[idx]);

                    polydata.GetCellData().SetScalars(Colors)
                    polydata.Modified()

                    out_file = "vtk_files/time_series/%s.%d.vtp"%(var_name, time_index)

                    writer = vtk.vtkXMLPolyDataWriter()
                    writer.SetFileName( out_file )
                    if vtk.VTK_MAJOR_VERSION <= 5:
                        writer.SetInput(polydata)
                    else:
                        writer.SetInputData(polydata)
                    writer.Write()

                    del Colors
                    del writer
                    del field
                    del field_ndims
                    del field_var
                    del dim_vals

            field_bar.update(time_index)

        nc_file.close()
        field_bar.finish()

        # Write pvd file, based on: http://www.paraview.org/Wiki/ParaView/Data_formats
        for var_name in written_fields:
            pvd_file = open('vtk_files/%s.pvd'%(var_name), 'w')
            pvd_file.write('<?xml version="1.0"?>\n')
            pvd_file.write('<VTKFile type="Collection" version="0.1"\n')
            pvd_file.write('\tbyte_order="LittleEndian"\n')
            pvd_file.write('\tcompressor="vtkZLibDataCompressor">\n')
            pvd_file.write('<Collection>\n')
            for time_index in time_global_indices:
                pvd_file.write('<DataSet timestep="%d" group="" part="0"\n'%(time_index))
                pvd_file.write('\tfile="time_series/%s.%d.vtp"/>\n'%(var_name, time_index))
            pvd_file.write('</Collection>\n')
            pvd_file.write('</VTKFile>\n')

        del written_fields
#}}}

def build_field_single_time_field( local_indices, global_indices, file_names, blocking, all_dim_vals, blockDimName, variable_list, point_list, cell_list ):#{{{
    if len(variable_list) > 0:
        polydata = vtk.vtkPolyData()
        polydata.SetPoints(point_list)
        polydata.SetPolys(cell_list)

        # Output time series
        widgets = ['Writing single fields: ', Percentage(), ' ', Bar(), ' ', ETA()]
        field_bar = ProgressBar(widgets=widgets, maxval=len(variable_list) ).start()
        prev_file = ""
        nc_file = NetCDFFile( file_names[0], 'r' )

        blockDim = len(nc_file.dimensions[blockDimName])
        nBlocks = 1 + blockDim / blocking

        var_index = 0
        for var_name in variable_list:
            Colors = vtk.vtkTypeFloat64Array()
            Colors.SetNumberOfComponents(1);
            Colors.SetName(var_name);

            dim_vals = all_dim_vals[var_name]
            field_var = nc_file.variables[var_name]
            field_dims = field_var.dimensions
            field_ndims = len(dim_vals)

            if not "Time" in field_dims:
                # Build data
                for iBlock in np.arange(0, nBlocks):
                    blockStart = iBlock * blocking
                    blockEnd = min( (iBlock + 1) * blocking, blockDim )
                    blockCount = blockEnd - blockStart

                    if field_ndims == 1:
                        field = field_var[blockStart:blockEnd]
                    if field_ndims == 2:
                        field = field_var[blockStart:blockEnd, dim_vals[1]]
                    elif field_ndims == 3:
                        field = field_var[blockStart:blockEnd, dim_vals[1]]
                    elif field_ndims == 4:
                        field = field_var[blockStart:blockEnd, dim_vals[1], dim_vals[2]]
                    elif field_ndims == 5:
                        field = field_var[blockStart:blockEnd, dim_vals[1], dim_vals[2], dim_vals[3]]

                    for idx in np.arange(0, blockCount):
                        Colors.InsertNextTuple1(field[idx]);

                polydata.GetCellData().SetScalars(Colors)
                polydata.Modified()

                out_file = "vtk_files/%s.vtp"%(var_name)

                writer = vtk.vtkXMLPolyDataWriter()
                writer.SetFileName( out_file )
                if vtk.VTK_MAJOR_VERSION <= 5:
                    writer.SetInput(polydata)
                else:
                    writer.SetInputData(polydata)
                writer.Write()

                del Colors
                del writer
                del field
                del field_ndims
                del field_var
                del dim_vals

            var_index = var_index + 1
            field_bar.update(var_index)

        nc_file.close()
        field_bar.finish()
#}}}

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-f", "--file_pattern", dest="filename_pattern", help="MPAS Filename pattern.", metavar="FILE", required=True)
    parser.add_argument("-b", "--blocking", dest="blocking", help="Size of blocks when reading MPAS file", metavar="BLK")
    parser.add_argument("-v", "--variable_list", dest="variable_list", help="List of variables to extract", metavar="VAR", required=True)

    args = parser.parse_args()

    if not args.blocking:
        args.blocking = int(10000)
    else:
        args.blocking = int(args.blocking)

    time_indices = []
    time_file_names = []
    time_global_indices = []
    setup_time_indices(args.filename_pattern, time_indices, time_global_indices, time_file_names)

    if not os.path.exists('vtk_files/time_series'):
        os.makedirs('vtk_files/time_series')

    # Setting dimension values:
    all_dim_vals = {}
    cellVarTime = []
    cellVarNoTime = []
    vertexVarTime = []
    vertexVarNoTime = []
    edgeVarTime = []
    edgeVarNoTime = []
    nc_file = NetCDFFile(time_file_names[0], 'r')
    setup_dimension_values_and_sort_vars( nc_file, args.variable_list, all_dim_vals, cellVarTime, cellVarNoTime, vertexVarTime, vertexVarNoTime, edgeVarTime, edgeVarNoTime)
    nc_file.close()

    summarize_extraction(all_dim_vals, cellVarTime, cellVarNoTime, vertexVarTime, vertexVarNoTime, edgeVarTime, edgeVarNoTime)

    # Handle cell variables
    if len(cellVarTime) > 0 or len(cellVarNoTime) > 0:
        print " -- Extracting cell fields --"
        CellVertices = vtk.vtkPoints()
        Cells = vtk.vtkCellArray()

        nc_file = NetCDFFile(time_file_names[0], 'r')
        # Build vertex list
        build_location_list_xyz( nc_file, args.blocking, 'nVertices', 'xVertex', 'yVertex', 'zVertex', CellVertices )

        # Build cell list
        build_cell_lists( nc_file, args.blocking, Cells )

        nc_file.close()

        build_field_single_time_field( time_indices, time_global_indices, time_file_names, args.blocking,
                                       all_dim_vals, 'nCells', cellVarNoTime, CellVertices, Cells )
        build_field_time_series( time_indices, time_global_indices, time_file_names, args.blocking,
                                 all_dim_vals, 'nCells', cellVarTime, CellVertices, Cells )

        print ""
        del Cells
        del CellVertices

    if len(vertexVarTime) > 0 or len(vertexVarNoTime) > 0:
        print " -- Extracting vertex fields --"
        VertexVertices = vtk.vtkPoints()
        Vertices = vtk.vtkCellArray()

        nc_file = NetCDFFile(time_file_names[0], 'r')
        # Build vertex list
        build_location_list_xyz( nc_file, args.blocking, 'nCells', 'xCell', 'yCell', 'zCell', VertexVertices )

        # Build cell list
        build_dual_cell_lists( nc_file, args.blocking, Vertices )

        nc_file.close()

        build_field_single_time_field( time_indices, time_global_indices, time_file_names, args.blocking,
                                       all_dim_vals, 'nVertices', vertexVarNoTime, VertexVertices, Vertices )
        build_field_time_series( time_indices, time_global_indices, time_file_names, args.blocking,
                                 all_dim_vals, 'nVertices', vertexVarTime, VertexVertices, Vertices )

        print ""
        del Vertices
        del VertexVertices

# vim: set expandtab:
