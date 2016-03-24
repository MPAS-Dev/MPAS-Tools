#!/usr/bin/env python
"""
Name: paraview_vtk_field_extractor.py
Author: Doug Jacobsen
Date: 03/01/2016

This script is used to extract a field from a time series of NetCDF files as
VTK files for plotting in paraview.

It can extract a field across multiple files by passing in a regular expression
for the filename patter. As an example, one can run the script using:

`./paraview_vtk_field_extractor.py -v areaCell,latVertex -f "hist.comp.*.nc"`

To extract a time series of areaCell,latVertex that spans multiple files.

Requirements:
This script requires access to the following non standard modules:
pyevtk
netCDF4
numpy

Optional modules:
progressbar
"""
from pyevtk.vtk import VtkFile, VtkPolyData
 
import sys, os, glob
import numpy as np

from netCDF4 import Dataset as NetCDFFile
import argparse

try:
    from progressbar import ProgressBar, Percentage, Bar, ETA
    use_progress_bar = True
except:
    use_progress_bar = False


def setup_time_indices(fn_pattern):#{{{
    # Build file list and time indices
    file_list = sorted(glob.glob(fn_pattern))
    
    local_indices = []
    global_indices = []
    file_names = []

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
        
    return (local_indices, global_indices, file_names)
#}}}

def setup_dimension_values_and_sort_vars( nc_file, variable_list):#{{{
    all_dim_vals = {}
    cellVars = []
    vertexVars = []
    edgeVars = []

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

        if "nCells" in field_dims:
            cellVars.append(var)
        elif "nVertices" in field_dims:
            vertexVars.append(var)
        elif "nEdges" in field_dims:
            edgeVars.append(var)

        all_dim_vals[var] = dim_vals
        del dim_vals

        if captured_input:
            print ""
    return (all_dim_vals, cellVars, vertexVars, edgeVars)
#}}}

def summarize_extraction(mesh_file, time_global_indices, cellVars, vertexVars, edgeVars):#{{{
    print ""
    print "Extracting a total of %d time levels."%(len(time_global_indices))
    print "Using file '%s' as the mesh file for this extraction."%(mesh_file)
    print ""
    print "The following variables will be extracted from the input file(s)."
    print ""

    if len(cellVars) > 0:
        print "   Variables with 'nCells' as a dimension:"
        for var in cellVars:
            print "      name: %s"%(var)

    if len(vertexVars) > 0:
        print "   Variables with 'nVertices' as a dimension:"
        for var in vertexVars:
            print "      name: %s"%(var)

    if len(edgeVars) > 0:
        print "   Variables with 'nEdges' as adimension:"
        for var in edgeVars:
            print "      name: %s"%(var)

    print ""
#}}}

def build_cell_lists( nc_file, blocking ):#{{{
    nCells = len(nc_file.dimensions['nCells'])

    nEdgesOnCell_var = nc_file.variables['nEdgesOnCell']
    verticesOnCell_var = nc_file.variables['verticesOnCell']
    
    offsets = np.cumsum(nEdgesOnCell_var[:], dtype=int)
    print offsets[0], offsets[-1]
    connectivity = np.zeros(offsets[-1], dtype=int)
    valid_mask = np.ones(nCells,bool)

    # Build cells
    nBlocks = 1 + nCells / blocking
    if use_progress_bar:
        widgets = ['Build cell connectivity: ', Percentage(), ' ', Bar(), ' ', ETA()]
        cell_bar = ProgressBar(widgets=widgets, maxval=nBlocks).start()
    else:
        print "Build cell connectivity..."

    outIndex = 0
    for iBlock in range(nBlocks):
        blockStart = iBlock * blocking
        blockEnd = min( (iBlock + 1) * blocking, nCells )
        blockCount = blockEnd - blockStart

        nEdgesOnCell = nEdgesOnCell_var[blockStart:blockEnd]
        verticesOnCell = verticesOnCell_var[blockStart:blockEnd,:] - 1
        
        for idx in range(blockCount):
            cellCount = nEdgesOnCell[idx]
            connectivity[outIndex:outIndex+cellCount] = verticesOnCell[idx,0:cellCount]
            outIndex += cellCount
 
        del nEdgesOnCell
        del verticesOnCell
        if use_progress_bar:
            cell_bar.update(iBlock)
    
    if use_progress_bar:
        cell_bar.finish()
    
    return (connectivity, offsets, valid_mask)

#}}}

def build_dual_cell_lists( nc_file, blocking ):#{{{
    print "Build dual connectivity...."
    vertexDegree = len(nc_file.dimensions['vertexDegree'])

    cellsOnVertex = nc_file.variables['cellsOnVertex'][:,:]-1
    
    valid_mask = np.all(cellsOnVertex >= 0, axis=1)
    connectivity = cellsOnVertex[valid_mask,:].ravel()
    validCount = np.count_nonzero(valid_mask)
    offsets = vertexDegree*np.arange(1,validCount+1)

    return (connectivity, offsets, valid_mask)

#}}}

def build_edge_cell_lists( nc_file, blocking ):#{{{
    nCells = len(nc_file.dimensions['nCells'])
    nEdges = len(nc_file.dimensions['nEdges'])

    cellsOnEdge_var = nc_file.variables['cellsOnEdge']
    verticesOnEdge_var = nc_file.variables['verticesOnEdge']

    valid_mask = np.zeros(nEdges,bool)
    
    connectivity = []
    offsets = []

    # Build cells
    nBlocks = 1 + nEdges / blocking
    if use_progress_bar:
        widgets = ['Build edge connectivity: ', Percentage(), ' ', Bar(), ' ', ETA()]
        edge_bar = ProgressBar(widgets=widgets, maxval=nBlocks).start()
    else:
        print "Build edge connectivity...."

    offset = 0
    for iBlock in np.arange(0, nBlocks):
        blockStart = iBlock * blocking
        blockEnd = min( (iBlock + 1) * blocking, nEdges )
        blockCount = blockEnd - blockStart

        verticesOnEdge = verticesOnEdge_var[blockStart:blockEnd,:] - 1
        cellsOnEdge = cellsOnEdge_var[blockStart:blockEnd,:] - 1
        
        vertices = np.zeros((blockCount,4))
        vertices[:,0] = cellsOnEdge[:,0]
        vertices[:,1] = verticesOnEdge[:,0]
        vertices[:,2] = cellsOnEdge[:,1]
        vertices[:,3] = verticesOnEdge[:,1]
        valid = vertices >= 0
        vertices[:,1] += nCells
        vertices[:,3] += nCells
        validCount = np.sum(np.array(valid,int),axis=1)
 
        for idx in range(blockCount):
            if(validCount[idx] < 3):
                continue
            valid_mask[blockStart + idx] = True
            verts = vertices[idx,valid[idx]]
            connectivity.extend(list(verts))
            offset += validCount[idx]
            offsets.append(offset)
       
        del cellsOnEdge
        del verticesOnEdge
        del vertices, valid, validCount
        
        if use_progress_bar:
            edge_bar.update(iBlock)

    if use_progress_bar:
        edge_bar.finish()

    connectivity = np.array(connectivity, dtype=int)
    offsets = np.array(offsets, dtype=int)

    return (connectivity, offsets, valid_mask)

#}}}

def build_location_list_xyz( nc_file, xName, yName, zName, output_32bit ):#{{{

    X = nc_file.variables[xName][:]
    Y = nc_file.variables[yName][:]
    Z = nc_file.variables[zName][:]
    if output_32bit:
        X = np.array(X,'f4')
        Y = np.array(Y,'f4')
        Z = np.array(Z,'f4')
    return (X,Y,Z)
    
#}}}
    
def build_field_time_series( local_indices, global_indices, file_names, blocking, all_dim_vals,
                             blockDimName, variable_list, vertices, connectivity, offsets,
                             valid_mask, output_32bit ):#{{{
    if len(variable_list) == 0:
        return
    
    if output_32bit:
        outType = 'f4'
    else:
        outType = 'f8'
        
                    
    # Get dimension info to allocate the size of Colors
    nc_file = NetCDFFile(file_names[0], 'r')

    blockDim = len(nc_file.dimensions[blockDimName])

    # Pre-compute the number of blocks
    nBlocks = 1 + blockDim / blocking

    nc_file.close
    
    var_has_time_dim = np.zeros(len(variable_list),bool)
    
    # Output time series
    if use_progress_bar:
        widgets = ['Writing time series: ', Percentage(), ' ', Bar(), ' ', ETA()]
        field_bar = ProgressBar(widgets=widgets, maxval=len(global_indices)*len(variable_list)*nBlocks ).start()
    else:
        print "Writing time series...."

    prev_file = ""
    for iTime in range(len(global_indices)):
        time_index = global_indices[iTime]

        if prev_file != file_names[time_index]:
            if prev_file != "":
                nc_file.close()
            nc_file = NetCDFFile(file_names[time_index], 'r')
            prev_file = file_names[time_index]

        for iVar in range(len(variable_list)):
            var_name = variable_list[iVar]
            dim_vals = all_dim_vals[var_name]
            field_var = nc_file.variables[var_name]
            var_has_time_dim[iVar] = 'Time' in field_var.dimensions
            
            if not var_has_time_dim[iVar] and iTime > 0:
                # already written out
                continue
            
            field_ndims = len(dim_vals)
            
            field = np.zeros(blockDim,dtype=outType)

            try:
                missing_val = field_var.missing_value
            except:
                missing_val = -9999999790214767953607394487959552.000000

            # Build data
            for iBlock in np.arange(0, nBlocks):
                blockStart = iBlock * blocking
                blockEnd = min( (iBlock + 1) * blocking, blockDim )

                if var_has_time_dim[iVar]:
                    assert(field_ndims != 1)
                    
                    if field_ndims == 2:
                        field_block = field_var[local_indices[time_index], blockStart:blockEnd]
                    elif field_ndims == 3:
                        field_block = field_var[local_indices[time_index], blockStart:blockEnd, dim_vals[2]]
                    elif field_ndims == 4:
                        field_block = field_var[local_indices[time_index], blockStart:blockEnd, dim_vals[2], dim_vals[3]]
                    elif field_ndims == 5:
                        field_block = field_var[local_indices[time_index], blockStart:blockEnd, dim_vals[2], dim_vals[3], dim_vals[4]]
                else:
                    if field_ndims == 1:
                        field_block = field_var[blockStart:blockEnd]
                    elif field_ndims == 2:
                        field_block = field_var[blockStart:blockEnd, dim_vals[1]]
                    elif field_ndims == 3:
                        field_block = field_var[blockStart:blockEnd, dim_vals[1], dim_vals[2]]
                    elif field_ndims == 4:
                        field_block = field_var[blockStart:blockEnd, dim_vals[1], dim_vals[2], dim_vals[3]]
                                        
                field_block[field_block == missing_val] = np.nan
                field[blockStart:blockEnd] = field_block

                if use_progress_bar:
                    field_bar.update((iTime*len(variable_list) + iVar)*nBlocks + iBlock)

            field = field[valid_mask]

            if var_has_time_dim[iVar]:
                file_index = time_index
            else:
                file_index = 0
            out_file_prefix = "vtk_files/time_series/%s.%d"%(var_name, file_index)

            npolys = len(field)
            vtkFile = VtkFile(out_file_prefix, VtkPolyData)
            vtkFile.openElement(vtkFile.ftype.name)
            vtkFile.openPiece(npoints=len(vertices[0]),npolys=npolys)
         
            vtkFile.openElement("Points")
            vtkFile.addData("points", vertices)
            vtkFile.closeElement("Points")
        
            vtkFile.openElement("Polys")
            vtkFile.addData("connectivity", connectivity)
            vtkFile.addData("offsets", offsets)
            vtkFile.closeElement("Polys")
           
            vtkFile.openData("Cell", scalars = var_name)
            vtkFile.addData(var_name, field)
            vtkFile.closeData("Cell")
            
            vtkFile.closePiece()
            vtkFile.closeElement(vtkFile.ftype.name)
            
            vtkFile.appendData(vertices)
            vtkFile.appendData(connectivity)
            vtkFile.appendData(offsets)
            vtkFile.appendData(field)
                  
            vtkFile.save()

            del vtkFile
            del field
            del field_ndims
            del field_var
            del dim_vals


    nc_file.close()
    if use_progress_bar:
        field_bar.finish()

    # Write pvd file, based on: http://www.paraview.org/Wiki/ParaView/Data_formats
    for iVar in range(len(variable_list)):
        var_name = variable_list[iVar]
        pvd_file = open('vtk_files/%s.pvd'%(var_name), 'w')
        pvd_file.write('<?xml version="1.0"?>\n')
        pvd_file.write('<VTKFile type="Collection" version="0.1"\n')
        pvd_file.write('\tbyte_order="LittleEndian"\n')
        pvd_file.write('\tcompressor="vtkZLibDataCompressor">\n')
        pvd_file.write('<Collection>\n')
        if var_has_time_dim[iVar]:
            for time_index in time_global_indices:
                pvd_file.write('<DataSet timestep="%d" group="" part="0"\n'%(time_index))
                pvd_file.write('\tfile="time_series/%s.%d.vtp"/>\n'%(var_name, time_index))
        else:
            pvd_file.write('<DataSet timestep="%d" group="" part="0"\n'%(0))
            pvd_file.write('\tfile="time_series/%s.%d.vtp"/>\n'%(var_name, 0))
        pvd_file.write('</Collection>\n')
        pvd_file.write('</VTKFile>\n')
#}}}


if __name__ == "__main__":
    if use_progress_bar:
        print " -- Using progress bars --"
    else:
        print " -- Progress bars are not available--"
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-f", "--file_pattern", dest="filename_pattern", help="MPAS Filename pattern.", metavar="FILE", required=True)
    parser.add_argument("-m", "--mesh_file", dest="mesh_filename", help="MPAS Mesh filename. If not set, it will use the first file in the -f flag as the mesh file.")
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

    (time_indices, time_global_indices, time_file_names) = setup_time_indices(args.filename_pattern)

    if not args.mesh_filename:
        args.mesh_filename = time_file_names[0]

    if not os.path.exists('vtk_files/time_series'):
        os.makedirs('vtk_files/time_series')

    # Setting dimension values:
    nc_file = NetCDFFile(time_file_names[0], 'r')
    (all_dim_vals, cellVars, vertexVars, edgeVars) = setup_dimension_values_and_sort_vars( nc_file, args.variable_list)
    nc_file.close()

    summarize_extraction(args.mesh_filename, time_global_indices, cellVars, vertexVars, edgeVars)

    # Handle cell variables
    if len(cellVars) > 0:
        print " -- Extracting cell fields --"

        nc_file = NetCDFFile(args.mesh_filename, 'r')

        # Build vertex list
        vertices = build_location_list_xyz( nc_file, 'xVertex', 'yVertex', 'zVertex', use_32bit )

        # Build cell list
        (connectivity, offsets, valid_mask) = build_cell_lists( nc_file, args.blocking )

        nc_file.close()

        build_field_time_series( time_indices, time_global_indices, time_file_names, args.blocking,
                                 all_dim_vals, 'nCells', cellVars, vertices, connectivity, offsets,
                                 valid_mask, use_32bit )

        print ""
        del vertices
        del connectivity
        del offsets
        del valid_mask

    if len(vertexVars) > 0:
        print " -- Extracting vertex fields --"

        nc_file = NetCDFFile(args.mesh_filename, 'r')

        # Build vertex list
        vertices = build_location_list_xyz( nc_file, 'xCell', 'yCell', 'zCell', use_32bit )

        # Build cell list
        (connectivity, offsets, valid_mask) = build_dual_cell_lists( nc_file, args.blocking )

        nc_file.close()

        build_field_time_series( time_indices, time_global_indices, time_file_names, args.blocking,
                                 all_dim_vals, 'nVertices', vertexVars, vertices, connectivity, offsets,
                                 valid_mask, use_32bit )

        print ""
        del vertices
        del connectivity
        del offsets
        del valid_mask

    if len(edgeVars) > 0:
        print " -- Extracting edge fields --"

        nc_file = NetCDFFile(args.mesh_filename, 'r')

        # Build vertex list
        verticesC = build_location_list_xyz( nc_file, 'xCell', 'yCell', 'zCell', use_32bit )
        verticesV = build_location_list_xyz( nc_file, 'xVertex', 'yVertex', 'zVertex', use_32bit )
        # compine the two into a single tuple
        vertices = (np.append(verticesC[0],verticesV[0]),
                    np.append(verticesC[1],verticesV[1]),
                    np.append(verticesC[2],verticesV[2]))
        del verticesC, verticesV

        # Build cell list
        (connectivity, offsets, valid_mask) = build_edge_cell_lists( nc_file, args.blocking )

        nc_file.close()

        build_field_time_series( time_indices, time_global_indices, time_file_names, args.blocking,
                                 all_dim_vals, 'nEdges', edgeVars, vertices, connectivity, offsets,
                                 valid_mask, use_32bit )


# vim: set expandtab:
