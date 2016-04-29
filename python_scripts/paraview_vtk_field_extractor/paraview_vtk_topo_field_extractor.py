#!/usr/bin/env python
"""
Name: paraview_vtk_topo_field_extractor.py
Authors: Xylar Asay-Davis
Date: 04/24/2016
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

        i_file = i_file + 1
        nc_file.close()
        if use_progress_bar:
            time_bar.update(i_file)

    if use_progress_bar:
        time_bar.finish()

    return (local_indices, file_names)
#}}}

def parse_extra_dim(dimName, indexString, time_series_file, mesh_file):#{{{
    if indexString == '':
        return np.zeros(0,int)
        
    if (mesh_file is not None) and (dimName in mesh_file.dimensions):
        nc_file = mesh_file
    else:
        nc_file = time_series_file
    dimSize = len(nc_file.dimensions[dimName])
    if ',' in indexString:
        indices = []
        for index in indexString.split(','):
            indices.append(int(index))
        indices = np.array(indices,int)
    elif ':' in indexString:
        index_list = indexString.split(':')
        if len(index_list) in [2,3]:
            if index_list[0] == '':
                first = 0
            else:
                first = int(index_list[0])
            if index_list[1] == '':
                last = dimSize
            else:
                last = int(index_list[1])
            if (len(index_list) == 2) or (index_list[2] == ''):
                step = 1
            else:
                step = int(index_list[2])
            indices = np.arange(first,last,step)
        else:
            print "Improperly formatted extra dimension:", dimName, indexString
            return None
    else:
        indices = np.array([int(indexString)])

    valid = np.all(np.logical_and(indices >= 0,indices <= dimSize-1))
    if not valid:
        print "Index (or indices) out of bounds for extra dimension:", dimName, indexString
        return None

    return indices

#}}}

def parse_extra_dims(topo_dim, topo_cell_indices_name, dimension_list, time_series_file, mesh_file):#{{{
    extra_dims = {}

    # topo_dim is a special extra dimension where we store an array of cell
    # indices to the bottom of the topography
    if topo_cell_indices_name:
        if (mesh_file is not None) and (topo_cell_indices_name in mesh_file.variables):
            topo_cell_indices = mesh_file.variables[topo_cell_indices_name][:]-1
        else:
            topo_cell_indices = time_series_file.variables[topo_cell_indices_name][:]-1
    else:
        index = len(mesh_file.dimensions[topo_dim])-1
        nCells = len(mesh_file.dimensions['nCells'])
        topo_cell_indices = index*np.ones(nCells, int)

    if dimension_list is not None:
        for dim_item in dimension_list:
            (dimName,indexString) = dim_item.split('=')
            indices = parse_extra_dim(dimName, indexString, time_series_file, mesh_file)
            if indices is not None:
                extra_dims[dimName] = indices
                
    return (topo_cell_indices, extra_dims)

#}}}

def setup_dimension_values_and_sort_vars(time_series_file, mesh_file, variable_list, topo_dim, extra_dims):#{{{
    all_dim_vals = {}
    cellVars = []

    if variable_list == 'all':
        variables = []
        files = [time_series_file]
        if mesh_file is not None:
            files.append(mesh_file)
        for nc_file in files:
            for var in nc_file.variables:
                dims = nc_file.variables[var].dimensions
                if 'nCells' in dims:
                    variables.append(str(var))
    else:
        variables = variable_list.split(',')

    # make sure the variables are unique
    variables = list(set(variables))

    promptDimNames = []
    for var in variables:
        if (mesh_file is not None) and (var in mesh_file.variables):
            nc_file = mesh_file
        else:
            nc_file = time_series_file
        field_dims = nc_file.variables[var].dimensions
        for dim in field_dims:
            if dim not in ['Time', 'nCells', 'nEdges', 'nVertices', topo_dim]:
                if dim not in extra_dims:
                    promptDimNames.append(str(dim))

    # make sure the extra dimension names are unique
    promptDimNames = list(set(promptDimNames))

    display_prompt = True
    for dim in promptDimNames:
        captured_input = True
        if display_prompt:
            print ""
            print "Need to define additional dimension values"
            display_prompt = False

        if (mesh_file is not None) and (dim in mesh_file.dimensions):
            nc_file = mesh_file
        else:
            nc_file = time_series_file
        dimSize = len(nc_file.dimensions[dim])
        valid = False
        while not valid:
            print "Valid range for dimension %s between 0 and %d"%(dim, dimSize-1)
            indexString = raw_input("Enter a value for dimension %s: "%(dim))
            indices = parse_extra_dim(dim, indexString, time_series_file, mesh_file)
            valid = indices is not None
            if valid:
                extra_dims[dim] = indices
            else:
                print " -- Invalid value, please re-enter --"


    for var in variables:
        captured_input = False

        dim_vals = []
        # first try the mesh file, then the time series
        if (mesh_file is not None) and (var in mesh_file.variables):
            nc_file = mesh_file
        else:
            nc_file = time_series_file
        field_var = nc_file.variables[var]

        # Setting dimension values:
        field_dims = field_var.dimensions
        indices = []
        for dim in field_dims:
            if dim not in ['Time', 'nCells', 'nEdges', 'nVertices', topo_dim]:
                indices.append(extra_dims[dim])
        if len(indices) == 0:
            dim_vals = None
        else:
            dim_vals = np.array(np.meshgrid(indices,indexing='ij'))

        cellVars.append(var)

        all_dim_vals[var] = dim_vals
        del dim_vals

        if captured_input:
            print ""
    return (all_dim_vals, cellVars)
#}}}

def summarize_extraction(mesh_file, time_indices, cellVars):#{{{
    print ""
    print "Extracting a total of %d time levels."%(len(time_indices))
    print "Using file '%s' as the mesh file for this extraction."%(mesh_file)
    print ""
    print "The following variables will be extracted from the input file(s)."
    print ""

    if len(cellVars) > 0:
        print "   Variables with 'nCells' as a dimension:"
        for var in cellVars:
            print "      name: %s"%(var)

    print ""
#}}}

def build_point_and_polygon_lists( nc_file, output_32bit ):#{{{
    nCells = len(nc_file.dimensions['nCells'])
    nEdges = len(nc_file.dimensions['nEdges'])

    nEdgesOnCell = nc_file.variables['nEdgesOnCell'][:]
    verticesOnCell = nc_file.variables['verticesOnCell'][:,:]-1
    edgesOnCell = nc_file.variables['edgesOnCell'][:,:]-1
    X_var = nc_file.variables['xVertex']
    Y_var = nc_file.variables['yVertex']
    Z_var = nc_file.variables['zVertex']

    nPoints = np.sum(nEdgesOnCell)
    nPolygons = nCells+nEdges

    if output_32bit:
        dtype = 'f4'
    else:
        dtype = 'f8'

    X = np.zeros(nPoints,dtype)
    Y = np.zeros(nPoints,dtype)
    Z = np.zeros(nPoints,dtype)
    pointToCellMap = np.zeros(nPoints,int)

    edgeVertToPointMap = -1*np.ones((nEdges,4),int)
    pointCountOnEdge = np.zeros((nEdges),int)

    offsets = np.zeros(nPolygons, dtype=int)
    # a polygon with nEdgesOnCell vertices per cell plus a polygon with 4 vertices per edge
    connectivity = np.zeros(nPoints+nEdges*4, dtype=int)

    # Build cells
    if use_progress_bar:
        widgets = ['Build cell connectivity: ', Percentage(), ' ', Bar(), ' ', ETA()]
        cell_bar = ProgressBar(widgets=widgets, maxval=nCells).start()
    else:
        print "Build cell connectivity..."

    outIndex = 0

    connectivity[0:nPoints] = np.arange(nPoints)

    for iCell in range(nCells):
        nVerts = nEdgesOnCell[iCell]
        iVerts = verticesOnCell[iCell,0:nVerts]
        iEdges = edgesOnCell[iCell,0:nVerts]
        
        X[outIndex:outIndex+nVerts] = X_var[iVerts]
        Y[outIndex:outIndex+nVerts] = Y_var[iVerts]
        Z[outIndex:outIndex+nVerts] = Z_var[iVerts]

        pointToCellMap[outIndex:outIndex+nVerts] = iCell
        for index in range(len(iEdges)):
            iEdge = iEdges[index]
            offset = pointCountOnEdge[iEdge]
            edgeVertToPointMap[iEdges[index], offset] = outIndex + np.mod(index-1, nVerts)
            edgeVertToPointMap[iEdges[index], offset+1] = outIndex + index
            pointCountOnEdge[iEdge] += 2
            

        outIndex += nVerts
        offsets[iCell] = outIndex

        if use_progress_bar:
            cell_bar.update(iCell)

    if use_progress_bar:
        cell_bar.finish()

    # Build cells
    if use_progress_bar:
        widgets = ['Build edge connectivity: ', Percentage(), ' ', Bar(), ' ', ETA()]
        cell_bar = ProgressBar(widgets=widgets, maxval=nEdges).start()
    else:
        print "Build cell connectivity..."

    for iEdge in range(nEdges):
        connectivity[outIndex] = edgeVertToPointMap[iEdge,0]
        connectivity[outIndex+1] = edgeVertToPointMap[iEdge,1]
        if(pointCountOnEdge[iEdge] > 2):
            connectivity[outIndex+2] = edgeVertToPointMap[iEdge,2]
            connectivity[outIndex+3] = edgeVertToPointMap[iEdge,3]
        else:
            connectivity[outIndex+2] = edgeVertToPointMap[iEdge,1]
            connectivity[outIndex+3] = edgeVertToPointMap[iEdge,0]
        
        outIndex += 4
        offsets[nCells + iEdge] = outIndex

        if use_progress_bar:
            cell_bar.update(iEdge)

    if use_progress_bar:
        cell_bar.finish()


    return ((X,Y,Z), connectivity, offsets, pointToCellMap)

#}}}


def build_field_time_series( local_time_indices, file_names, mesh_file, all_dim_vals,
                             variable_list, vertices, connectivity, offsets,
                             topo_dim, topo_cell_indices, pointToCellMap, 
                             output_32bit, combine_output, append ):#{{{
    def write_pvd_header(prefix):
        pvd_file = open('vtk_files/%s.pvd'%(prefix), 'w')
        pvd_file.write('<?xml version="1.0"?>\n')
        pvd_file.write('<VTKFile type="Collection" version="0.1"\n')
        pvd_file.write('\tbyte_order="LittleEndian"\n')
        pvd_file.write('\tcompressor="vtkZLibDataCompressor">\n')
        return pvd_file

    def get_hyperslab_name_and_dims(var_name):
        extra_dim_vals = all_dim_vals[var_name]
        if(extra_dim_vals is None):
            return ([var_name],None)
        if(extra_dim_vals.size == 0):
            return ([],None)
        out_var_names = []
        maxval = np.amax(extra_dim_vals,axis=1)
        pad = np.ones(maxval.shape,int)
        mask = maxval > 0
        pad[mask] = np.array(np.floor(np.log10(maxval[mask])),int)+1
        for iHyperSlab in range(extra_dim_vals.shape[1]):
            out_var_name = var_name
            for iVal in range(extra_dim_vals.shape[0]):
                val = extra_dim_vals[iVal,iHyperSlab]
                template = '%%s_%%0%dd'%pad[iVal]
                out_var_name = template%(out_var_name,val)
                out_var_names.append(out_var_name)
        return (out_var_names, extra_dim_vals)

    def write_vtp_header(prefix, activeVarIndex, varIndices):
        vtkFile = VtkFile("vtk_files/%s"%prefix, VtkPolyData)
        vtkFile.openElement(vtkFile.ftype.name)
        vtkFile.openPiece(npoints=nPoints,npolys=nPolygons)

        vtkFile.openElement("Points")
        vtkFile.addData("points", vertices)
        vtkFile.closeElement("Points")

        vtkFile.openElement("Polys")
        vtkFile.addData("connectivity", connectivity)
        vtkFile.addData("offsets", offsets)
        vtkFile.closeElement("Polys")

        vtkFile.openData("Point", scalars = variable_list[activeVarIndex])
        for iVar in varIndices:
            var_name = variable_list[iVar]
            (out_var_names, dim_list) = get_hyperslab_name_and_dims(var_name)
            for out_var_name in out_var_names:
                vtkFile.addHeader(out_var_name, outType, nPoints, 1)
        vtkFile.closeData("Point")

        vtkFile.closePiece()
        vtkFile.closeElement(vtkFile.ftype.name)

        vtkFile.appendData(vertices)
        vtkFile.appendData(connectivity)
        vtkFile.appendData(offsets)

        return vtkFile


    if len(variable_list) == 0:
        return

    if output_32bit:
        outType = 'float32'
    else:
        outType = 'float64'


    # Get dimension info to allocate the size of Colors
    time_series_file = NetCDFFile(file_names[0], 'r')

    if mesh_file is not None:
        # nGeomDim may not exist in time series file
        nCells = len(mesh_file.dimensions['nCells'])
    else:
        nCells = len(time_series_file.dimensions['nCells'])
        
    if (mesh_file is not None) and (topo_dim in mesh_file.dimensions):
        nTopoLevels = len(mesh_file.dimensions[topo_dim])
    else:
        nTopoLevels = len(time_series_file.dimensions[topo_dim])

    nPolygons = len(offsets)
    nPoints = len(vertices[0])
    nTimes = len(local_time_indices)
    nVars = len(variable_list)
    
    var_has_time_dim = np.zeros(nVars,bool)
    nHyperSlabs = 0
    for iVar in range(nVars):
        var_name = variable_list[iVar]
        if (mesh_file is not None) and (var_name in mesh_file.variables):
            # we can't support time dependence in the mesh file
            assert('Time' not in mesh_file.variables[var_name].dimensions)
            var_has_time_dim[iVar] = False
        else:
            var_has_time_dim[iVar] = 'Time' in time_series_file.variables[var_name].dimensions

        extra_dim_vals = all_dim_vals[var_name]
        if (extra_dim_vals is None) or (extra_dim_vals.size == 0):
            nHyperSlabs += 1
        else:
            nHyperSlabs += extra_dim_vals.shape[1]

    time_series_file.close()

    any_var_has_time_dim = np.any(var_has_time_dim)

    try:
        os.makedirs('vtk_files')
    except OSError:
        pass

    if any_var_has_time_dim:
        try:
            os.makedirs('vtk_files/time_series')
        except OSError:
            pass
    else:
        # there is no point in combining output if no fields have Time dim
        combine_output = False
        nTimes = 1

    # Output time series
    if use_progress_bar:
        widgets = ['Writing time series: ', Percentage(), ' ', Bar(), ' ', ETA()]
        field_bar = ProgressBar(widgets=widgets, maxval=nTimes*nHyperSlabs ).start()
    else:
        print "Writing time series...."

    suffix = 'Cells'
    if any_var_has_time_dim:
        if combine_output or np.all(var_has_time_dim):
            out_prefix = "topoFieldsOn%s"%suffix
        else:
            out_prefix = "timeDependentTopoFieldsOn%s"%suffix
        # start the pvd file
        pvd_file = write_pvd_header(out_prefix)
        pvd_file.write('<Collection>\n')

    if not combine_output and not np.all(var_has_time_dim):
        out_prefix = "staticTopoFieldsOn%s"%suffix
        varIndices = np.arange(nVars)[var_has_time_dim == False]
        timeIndependentFile = write_vtp_header(out_prefix, varIndices[0], varIndices)


    prev_file = ""
    for time_index in range(nTimes):

        if prev_file != file_names[time_index]:
            if prev_file != "":
                time_series_file.close()
            time_series_file = NetCDFFile(file_names[time_index], 'r')
            prev_file = file_names[time_index]

        if any_var_has_time_dim:
            # write the header for the vtp file
            vtp_file_prefix = "time_series/%s.%d"%(out_prefix, time_index)
            file_name = 'vtk_files/%s.vtp'%(vtp_file_prefix)
            if append and os.path.exists(file_name):
                continue

            if combine_output:
                varIndices = np.arange(nVars)
            else:
                varIndices = np.arange(nVars)[var_has_time_dim]
            timeDependentFile = write_vtp_header(vtp_file_prefix, varIndices[0], varIndices)

            # add time step to pdv file
            pvd_file.write('<DataSet timestep="%d" group="" part="0"\n'%(time_index))
            pvd_file.write('\tfile="%s.vtp"/>\n'%(vtp_file_prefix))

        if time_index == 0:
            varIndices = np.arange(nVars)
        else:
            # only the time-dependent variables
            varIndices = np.arange(nVars)[var_has_time_dim]


        iHyperSlabProgress = 0
        for iVar in varIndices:
            has_time = var_has_time_dim[iVar]
            if not has_time and time_index > 0:
                continue

            var_name = variable_list[iVar]
            if (mesh_file is not None) and (var_name in mesh_file.variables):
                nc_file = mesh_file
            else:
                nc_file = time_series_file
            field_var = nc_file.variables[var_name]

            field_ndims = len(field_var.dimensions)
            

            (out_var_names, dim_list) = get_hyperslab_name_and_dims(var_name)                    
                        
            if has_time or combine_output:
                vtkFile = timeDependentFile
            else:
                vtkFile = timeIndependentFile
            for iHyperSlab in range(len(out_var_names)):
                
                try:
                    missing_val = field_var.missing_value
                except:
                    missing_val = -9999999790214767953607394487959552.000000

                # Build data
                dim_vals = []
                extraDimIndex = 0
                for iDim in range(field_ndims):
                    dim = field_var.dimensions[iDim]
                    if dim == 'Time':
                        dim_vals.append(local_time_indices[time_index])
                    elif dim == 'nCells':
                        dim_vals.append(np.arange(nCells))
                    elif dim == topo_dim:
                        dim_vals.append(np.arange(nTopoLevels))
                    else:
                        dim_vals.append(dim_list[extraDimIndex,iHyperSlab])
                        extraDimIndex += 1
                        
                if field_ndims == 1:
                    field = field_var[dim_vals[0]]
                elif field_ndims == 2:
                    field = field_var[dim_vals[0], dim_vals[1]]
                elif field_ndims == 3:
                    field = field_var[dim_vals[0], dim_vals[1], dim_vals[2]]
                elif field_ndims == 4:
                    field = field_var[dim_vals[0], dim_vals[1], dim_vals[2], dim_vals[3]]
                elif field_ndims == 5:
                    field = field_var[dim_vals[0], dim_vals[1], dim_vals[2], dim_vals[3], dim_vals[4]]

                if topo_dim in field_var.dimensions:
                    field = field[np.arange(nCells),topo_cell_indices]
                
                # convert to the same type as field before masking with NaNs
                field = np.array(field, dtype=outType)
                field[field == missing_val] = np.nan
                
                # map field from cells to points
                field = field[pointToCellMap]

                vtkFile.appendData(field)

                if use_progress_bar:
                    field_bar.update(time_index*nHyperSlabs + iHyperSlabProgress)
                    iHyperSlabProgress += 1

                del field
                del field_ndims
                del field_var

        if any_var_has_time_dim:
            timeDependentFile.save()
            del timeDependentFile

        if time_index == 0 and not combine_output and not np.all(var_has_time_dim):
            timeIndependentFile.save()
            del timeIndependentFile

    time_series_file.close()
    if use_progress_bar:
        field_bar.finish()

    if any_var_has_time_dim:
        # finish the pdv file
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
    parser.add_argument("-t", "--topo_dim", dest="topo_dim", help="Dimension and range for topography dimension", required=True)
    parser.add_argument("-i", "--topo_cell_index", dest="topo_cell_index", help="Index array indicating the bottom of the domain (default is the topo_dim-1 for all cells)")
    parser.add_argument("-v", "--variable_list", dest="variable_list", help="List of variables on cells to extract ('all' for all variables on cells)", metavar="VAR", required=True)
    parser.add_argument("-3", "--32bit", dest="output_32bit", help="If set, the vtk files will be written using 32bit floats.", action="store_true")
    parser.add_argument("-c", "--combine", dest="combine_output", help="If set, time-independent fields are written to each file along with time-dependent fields.", action="store_true")
    parser.add_argument("-a", "--append", dest="append", help="If set, only vtp files that do not already exist are written out.", action="store_true")
    parser.add_argument("-d", "--dim_list", dest="dimension_list", nargs="+", help="A list of dimensions and associated values.")
    args = parser.parse_args()

    if not args.output_32bit:
        use_32bit = False
    else:
        use_32bit = True

    (time_indices, time_file_names) = setup_time_indices(args.filename_pattern)

    separate_mesh_file = True
    if not args.mesh_filename:
        args.mesh_filename = time_file_names[0]
        separate_mesh_file = False

    # Setting dimension values:
    time_series_file = NetCDFFile(time_file_names[0], 'r')
    if separate_mesh_file:
        mesh_file = NetCDFFile(args.mesh_filename, 'r')
    else:
        mesh_file = None
    (topo_cell_indices, extra_dims) = parse_extra_dims(args.topo_dim, 
        args.topo_cell_index, args.dimension_list, time_series_file, mesh_file)
    (all_dim_vals, cellVars) = setup_dimension_values_and_sort_vars(
            time_series_file, mesh_file, args.variable_list, args.topo_dim, extra_dims)
    time_series_file.close()
    if(mesh_file is not None):
        mesh_file.close()

    summarize_extraction(args.mesh_filename, time_indices, cellVars)

    # Handle cell variables
    if len(cellVars) > 0:
        print " -- Extracting cell fields --"

        mesh_file = NetCDFFile(args.mesh_filename, 'r')

        # Build vertex list
        # Build cell list
        (vertices, connectivity, offsets, pointToCellMap) = build_point_and_polygon_lists( 
           mesh_file, use_32bit )

        if not separate_mesh_file:
            mesh_file.close()
            mesh_file = None

        build_field_time_series( time_indices, time_file_names, mesh_file,
                                 all_dim_vals, cellVars, vertices, connectivity, offsets,
                                 args.topo_dim, topo_cell_indices, pointToCellMap, 
                                 use_32bit, args.combine_output, args.append )
        if separate_mesh_file:
            mesh_file.close()

        print ""
        del vertices
        del connectivity
        del offsets


# vim: set expandtab:
