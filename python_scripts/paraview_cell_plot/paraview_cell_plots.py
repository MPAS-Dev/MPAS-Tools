#!/usr/bin/env python
"""
Name: paraview_cell_plots.py
Author: Doug Jacobsen
Date: 03/01/2016

This script is used to extract a field from a time series of NetCDF files as
VTK files for plotting in paraview.

It can extract a field across multiple files by passing in a regular expression
for the filename patter. As an example, one can run the script using:

`./paraview_cell_plots.py -v areaCell,latVertex -f "hist.comp.*.nc"`

To extract a time series of areaCell,latVertex that spans multiple files.

Requirements:
This script requires access to the following non standard modules:
vtk
netCDF4
numpy

Optional modules:
progressbar
"""
import vtk
 
import sys, os, glob
import numpy as np

from netCDF4 import *
from netCDF4 import Dataset as NetCDFFile
import argparse

try:
    from progressbar import *
    use_progress_bar = True
except:
    use_progress_bar = False

def setup_time_indices(fn_pattern, local_indices, global_indices, file_names):#{{{
    # Build file list and time indices
    file_list = sorted(glob.glob(fn_pattern))

    if len(file_list) == 0:
        print "No files to process."
        print "Exiting..."
        sys.exit(0)

    if use_progress_bar:
        widgets = ['Build time indices: ', Percentage(), ' ', Bar(), ' ', ETA()]
        time_bar = ProgressBar(widgets=widgets, maxval=len(file_list)).start()
    else:
        print "Build time indices..."

    i_global = 0
    i_file = 0
    for file in file_list:
        nc_file = NetCDFFile(file, 'r')

        try:
            times = len(nc_file.dimensions['Time'])
        except:
            times = 1

        times = max(times, 1)

        for time_idx in np.arange(0, times):
            local_indices.append(time_idx)
            file_names.append(file)
            global_indices.append(i_global)

            i_global = i_global + 1

        i_file = i_file + 1
        nc_file.close()
        if use_progress_bar:
            time_bar.update(i_file)

    if use_progress_bar:
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

def summarize_extraction(all_dim_vals, time_global_indices, cellVarTime, cellVarNoTime, vertexVarTime, vertexVarNoTime, edgeVarTime, edgeVarNoTime):#{{{
    print ""
    print "Extracting a total of %d time levels."%(len(time_global_indices))
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

def build_cell_lists( nc_file, blocking, cell_list ):#{{{
    nCells = len(nc_file.dimensions['nCells'])

    nEdgesOnCell_var = nc_file.variables['nEdgesOnCell']
    verticesOnCell_var = nc_file.variables['verticesOnCell']

    # Build cells
    nBlocks = 1 + nCells / blocking
    if use_progress_bar:
        widgets = ['Build cell connectivity: ', Percentage(), ' ', Bar(), ' ', ETA()]
        cell_bar = ProgressBar(widgets=widgets, maxval=nBlocks).start()
    else:
        print "Build cell connectivity..."

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
        if use_progress_bar:
            cell_bar.update(iBlock)

    if use_progress_bar:
        cell_bar.finish()

#}}}

def build_dual_cell_lists( nc_file, blocking, cell_list ):#{{{
    nVertices = len(nc_file.dimensions['nVertices'])
    vertexDegree = len(nc_file.dimensions['vertexDegree'])

    cellsOnVertex_var = nc_file.variables['cellsOnVertex']

    # Build cells
    nBlocks = 1 + nVertices / blocking
    if use_progress_bar:
        widgets = ['Build dual connectivity: ', Percentage(), ' ', Bar(), ' ', ETA()]
        dual_bar = ProgressBar(widgets=widgets, maxval=nBlocks).start()
    else:
        print "Build dual connectivity...."

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
        if use_progress_bar:
            dual_bar.update(iBlock)

    if use_progress_bar:
        dual_bar.finish()
#}}}

def build_edge_cell_lists( nc_file, blocking, cell_list ):#{{{
    nCells = len(nc_file.dimensions['nCells'])
    nEdges = len(nc_file.dimensions['nEdges'])
    nVertices = len(nc_file.dimensions['nVertices'])

    cellsOnEdge_var = nc_file.variables['cellsOnEdge']
    verticesOnEdge_var = nc_file.variables['verticesOnEdge']

    # Build cells
    nBlocks = 1 + nEdges / blocking
    if use_progress_bar:
        widgets = ['Build edge connectivity: ', Percentage(), ' ', Bar(), ' ', ETA()]
        edge_bar = ProgressBar(widgets=widgets, maxval=nBlocks).start()
    else:
        print "Build edge connectivity...."

    for iBlock in np.arange(0, nBlocks):
        blockStart = iBlock * blocking
        blockEnd = min( (iBlock + 1) * blocking, nEdges )
        blockCount = blockEnd - blockStart

        verticesOnEdge = verticesOnEdge_var[blockStart:blockEnd][:]
        cellsOnEdge = cellsOnEdge_var[blockStart:blockEnd][:]

        for idx in np.arange(0, blockCount):
            polygon = vtk.vtkPolygon()
            polygon.GetPointIds().SetNumberOfIds(4)

            cell1 = cellsOnEdge[idx][0] - 1
            vertex1 = verticesOnEdge[idx][0] - 1
            cell2 = cellsOnEdge[idx][1] - 1
            vertex2 = verticesOnEdge[idx][1] - 1

            land_cell = False
            if ( cell1 < 0 ):
                land_cell = True

            if ( cell2 < 0 ):
                land_cell = True

            if ( vertex1 < 0 ):
                land_cell = True

            if ( vertex2 < 0 ):
                land_cell = True

            if land_cell:
                cell1 = nCells
                cell2 = nCells
                vertex1 = nVertices
                vertex2 = nVertices

            polygon.GetPointIds().SetId(0, cell1)
            polygon.GetPointIds().SetId(1, vertex1 + nCells + 1) # nCells has an extra garbage cell at the end
            polygon.GetPointIds().SetId(2, cell2)
            polygon.GetPointIds().SetId(3, vertex2 + nCells + 1) # nCells has an extra garbage cell at the end

            cell_list.InsertNextCell(polygon)

        del cellsOnEdge
        del verticesOnEdge
        if use_progress_bar:
            edge_bar.update(iBlock)

    if use_progress_bar:
        edge_bar.finish()
#}}}

def build_location_list_xyz( nc_file, blocking, dimName, xName, yName, zName, point_list ):#{{{
    nLocs = len(nc_file.dimensions[dimName])

    xLoc_var = nc_file.variables[xName]
    yLoc_var = nc_file.variables[yName]
    zLoc_var = nc_file.variables[zName]

    # Build vertex list
    nBlocks = 1 + nLocs / blocking
    if use_progress_bar:
        widgets = ['Build location list: ', Percentage(), ' ', Bar(), ' ', ETA()]
        loc_bar = ProgressBar(widgets=widgets, maxval=nBlocks).start()
    else:
        print "Build location list...."

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
        if use_progress_bar:
            loc_bar.update(iBlock)

    id = point_list.InsertNextPoint(0.0, 0.0, 0.0)

    if use_progress_bar:
        loc_bar.finish()
#}}}

def build_field_time_series( local_indices, global_indices, file_names, blocking, all_dim_vals, blockDimName, variable_list, point_list, cell_list, output_32bit ):#{{{
    if len(variable_list) > 0:
        polydata = vtk.vtkPolyData()
        polydata.SetPoints(point_list)
        polydata.SetPolys(cell_list)

        if output_32bit:
            Colors = vtk.vtkTypeFloat32Array()
        else:
            Colors = vtk.vtkTypeFloat64Array()

        Colors.SetNumberOfComponents(1);

        # Get dimension info to allocate the size of Colors
        nc_file = NetCDFFile(file_names[0], 'r')

        blockDim = len(nc_file.dimensions[blockDimName])

        # Pre-compute the number of blocks
        nBlocks = 1 + blockDim / blocking

        nc_file.close

        # Allocate the correct number of elements for the colors array
        Colors.Allocate( blockDim+1 )

        for i in np.arange(0, blockDim+1):
            Colors.InsertNextTuple1(0.0);

        # Output time series
        if use_progress_bar:
            widgets = ['Writing time series: ', Percentage(), ' ', Bar(), ' ', ETA()]
            field_bar = ProgressBar(widgets=widgets, maxval=len(global_indices) ).start()
        else:
            print "Writing time series...."

        prev_file = ""
        for time_index in global_indices:

            if prev_file != file_names[time_index]:
                if prev_file != "":
                    nc_file.close()
                nc_file = NetCDFFile(file_names[time_index], 'r')
                prev_file = file_names[time_index]

            for var_name in variable_list:
                Colors.SetName(var_name);
                dim_vals = all_dim_vals[var_name]
                field_var = nc_file.variables[var_name]
                field_ndims = len(dim_vals)
                field_dims = field_var.dimensions

                try:
                    missing_val = field_var.missing_value
                except:
                    missing_val = -9999999790214767953607394487959552.000000

                # Build data
                for iBlock in np.arange(0, nBlocks):
                    blockStart = iBlock * blocking
                    blockEnd = min( (iBlock + 1) * blocking, blockDim )
                    blockCount = blockEnd - blockStart

                    if field_ndims == 1:
                        # We should never encounter a 1D field, since it doesn't have time as a dimension.
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
                        if (field[idx] == missing_val):
                            Colors.SetTuple1(blockStart + idx, float('nan'));
                        else:
                            Colors.SetTuple1(blockStart + idx, field[idx]);

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

                del writer
                del field
                del field_ndims
                del field_var
                del dim_vals

            if use_progress_bar:
                field_bar.update(time_index)

        del Colors
        nc_file.close()
        if use_progress_bar:
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
#}}}

def build_field_single_time_field( local_indices, global_indices, file_names, blocking, all_dim_vals, blockDimName, variable_list, point_list, cell_list, output_32bit ):#{{{
    if len(variable_list) > 0:
        polydata = vtk.vtkPolyData()
        polydata.SetPoints(point_list)
        polydata.SetPolys(cell_list)

        # Output time series
        if use_progress_bar:
            widgets = ['Writing single fields: ', Percentage(), ' ', Bar(), ' ', ETA()]
            field_bar = ProgressBar(widgets=widgets, maxval=len(variable_list) ).start()
        else:
            print "Writing single fields...."

        prev_file = ""
        nc_file = NetCDFFile( file_names[0], 'r' )

        blockDim = len(nc_file.dimensions[blockDimName])
        nBlocks = 1 + blockDim / blocking

        if output_32bit:
            Colors = vtk.vtkTypeFloat32Array()
        else:
            Colors = vtk.vtkTypeFloat64Array()

        Colors.SetNumberOfComponents(1);

        # Get dimension info to allocate the size of Colors
        nc_file = NetCDFFile(file_names[0], 'r')

        blockDim = len(nc_file.dimensions[blockDimName])

        # Pre-compute the number of blocks
        nBlocks = 1 + blockDim / blocking

        nc_file.close

        # Allocate the correct number of elements for the colors array
        Colors.Allocate( blockDim )

        for i in np.arange(0, blockDim+1):
            Colors.InsertNextTuple1(0.0);

        var_index = 0
        for var_name in variable_list:
            Colors.SetName(var_name);

            dim_vals = all_dim_vals[var_name]
            field_var = nc_file.variables[var_name]
            field_ndims = len(dim_vals)

            try:
                missing_val = field_var.missing_value
            except:
                missing_val = -9999999790214767953607394487959552.000000

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
                    if ( field[idx] == missing_val  ):
                        Colors.SetTuple1(blockStart + idx, float('nan'))
                    else:
                        Colors.SetTuple1(blockStart + idx, field[idx])

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

            del writer
            del field
            del field_ndims
            del field_var
            del dim_vals

            var_index = var_index + 1
            if use_progress_bar:
                field_bar.update(var_index)

        del Colors
        nc_file.close()
        if use_progress_bar:
            field_bar.finish()
#}}}

if __name__ == "__main__":
    if use_progress_bar:
        print " -- Using progress bars --"
    else:
        print " -- Progress bars are not available--"
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-f", "--file_pattern", dest="filename_pattern", help="MPAS Filename pattern.", metavar="FILE", required=True)
    parser.add_argument("-b", "--blocking", dest="blocking", help="Size of blocks when reading MPAS file", metavar="BLK")
    parser.add_argument("-v", "--variable_list", dest="variable_list", help="List of variables to extract", metavar="VAR", required=True)
    parser.add_argument("-3", "--32bit", dest="output_32bit", help="If set, the vtk files will be written using 32bit floats.", action="store_true")

    args = parser.parse_args()

    if not args.output_32bit:
        use_32bit = False
    else:
        use_32bit = True

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

    summarize_extraction(all_dim_vals, time_global_indices, cellVarTime, cellVarNoTime, vertexVarTime, vertexVarNoTime, edgeVarTime, edgeVarNoTime)

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
                                       all_dim_vals, 'nCells', cellVarNoTime, CellVertices, Cells, use_32bit )
        build_field_time_series( time_indices, time_global_indices, time_file_names, args.blocking,
                                 all_dim_vals, 'nCells', cellVarTime, CellVertices, Cells, use_32bit )

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
                                       all_dim_vals, 'nVertices', vertexVarNoTime, VertexVertices, Vertices, use_32bit )
        build_field_time_series( time_indices, time_global_indices, time_file_names, args.blocking,
                                 all_dim_vals, 'nVertices', vertexVarTime, VertexVertices, Vertices, use_32bit )

        print ""
        del Vertices
        del VertexVertices

    if len(edgeVarTime) > 0 or len(edgeVarNoTime) > 0:
        print " -- Extracting edge fields --"
        EdgeVertices = vtk.vtkPoints()
        Edges = vtk.vtkCellArray()

        nc_file = NetCDFFile(time_file_names[0], 'r')

        # Build vertex list
        build_location_list_xyz( nc_file, args.blocking, 'nCells', 'xCell', 'yCell', 'zCell', EdgeVertices )
        build_location_list_xyz( nc_file, args.blocking, 'nVertices', 'xVertex', 'yVertex', 'zVertex', EdgeVertices )

        # Build cell list
        build_edge_cell_lists( nc_file, args.blocking, Edges )

        nc_file.close()

        build_field_single_time_field( time_indices, time_global_indices, time_file_names, args.blocking,
                                       all_dim_vals, 'nEdges', edgeVarNoTime, EdgeVertices, Edges, use_32bit )
        build_field_time_series( time_indices, time_global_indices, time_file_names, args.blocking,
                                 all_dim_vals, 'nEdges', edgeVarTime, EdgeVertices, Edges, use_32bit )


# vim: set expandtab:
