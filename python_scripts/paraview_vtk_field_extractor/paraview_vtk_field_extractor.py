#!/usr/bin/env python
"""
Name: paraview_vtk_field_extractor.py
Authors: Doug Jacobsen, Xylar Asay-Davis
Date: 03/01/2016

This script is used to extract a field from a time series of NetCDF files as
VTK files for plotting in paraview.

It can extract a field across multiple files by passing in a regular expression
for the filename patter. As an example, one can run the script using:

`./paraview_vtk_field_extractor.py -v areaCell,latVertex -f "hist.comp.*.nc"`

To extract a time series of areaCell,latVertex that spans multiple files.
By default, time-independent fields on cells are written to a file
vtk_files/staticFieldsOnCells.vtp
and time-dependent fields on cells are written to
vtk_files/timeDependentFieldsOnCells.pvd
vtk_files/time_series/timeDependentFieldsOnCells.0.vtp
vtk_files/time_series/timeDependentFieldsOnCells.1.vtp
...
and similarly for edges and vertices.  Time-independent fields can be
included in each time step of the time-dependent fields for with
the --combine flag.  This allows time-dependent and -independent fields
to be combined in filters within Paraview at the expense of considerable
additional storage space.


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

def parse_extra_dims(dimension_list, time_series_file, mesh_file):#{{{
    if not dimension_list:
        return {}

    extra_dims = {}
    for dim_item in dimension_list:
        (dimName,indexString) = dim_item.split('=')
        indices = parse_extra_dim(dimName, indexString, time_series_file, mesh_file)
        if indices is not None:
            extra_dims[dimName] = indices
    return extra_dims

#}}}

def setup_dimension_values_and_sort_vars(time_series_file, mesh_file, variable_list, extra_dims):#{{{
    all_dim_vals = {}
    cellVars = []
    vertexVars = []
    edgeVars = []

    if variable_list == 'all':
        variables = []
        files = [time_series_file]
        if mesh_file is not None:
            files.append(mesh_file)
        for nc_file in files:
            for var in nc_file.variables:
                dims = nc_file.variables[var].dimensions
                supported = False
                for dim in ['nCells', 'nEdges', 'nVertices']:
                    if dim in dims:
                        supported = True
                if supported:
                    variables.append(str(var))
        # make sure the variables are unique
        variables = list(set(variables))
    else:
        variables = variable_list.split(',')

    extraDimNames = []
    promptDimNames = []
    for var in variables:
        if (mesh_file is not None) and (var in mesh_file.variables):
            nc_file = mesh_file
        else:
            nc_file = time_series_file
        field_dims = nc_file.variables[var].dimensions
        for dim in field_dims:
            if dim not in ['Time', 'nCells', 'nEdges', 'nVertices']:
                extraDimNames.append(str(dim))
                if dim not in extra_dims:
                    promptDimNames.append(str(dim))

    # make sure the extra dimension names are unique
    extraDimNames = list(set(extraDimNames))
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
            if dim not in ['Time', 'nCells', 'nEdges', 'nVertices']:
                indices.append(extra_dims[dim])
        if len(indices) == 0:
            dim_vals = None
        else:
            dim_vals = np.array(np.meshgrid(indices,indexing='ij'))

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

def summarize_extraction(mesh_file, time_indices, cellVars, vertexVars, edgeVars):#{{{
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

def build_field_time_series( local_time_indices, file_names, mesh_file, blocking, all_dim_vals,
                             blockDimName, variable_list, vertices, connectivity, offsets,
                             valid_mask, output_32bit, combine_output, append ):#{{{
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

        vtkFile.openData("Cell", scalars = variable_list[activeVarIndex])
        for iVar in varIndices:
            var_name = variable_list[iVar]
            (out_var_names, dim_list) = get_hyperslab_name_and_dims(var_name)
            for out_var_name in out_var_names:
                vtkFile.addHeader(out_var_name, outType, nPolygons, 1)
        vtkFile.closeData("Cell")

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
        # blockDim may not exist in time series file
        blockDim = len(mesh_file.dimensions[blockDimName])
    else:
        blockDim = len(time_series_file.dimensions[blockDimName])

    # Pre-compute the number of blocks
    nBlocks = 1 + blockDim / blocking


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

    suffix = blockDimName[1:]
    if any_var_has_time_dim:
        if combine_output or np.all(var_has_time_dim):
            out_prefix = "fieldsOn%s"%suffix
        else:
            out_prefix = "timeDependentFieldsOn%s"%suffix
        # start the pvd file
        pvd_file = write_pvd_header(out_prefix)
        pvd_file.write('<Collection>\n')

    if not combine_output and not np.all(var_has_time_dim):
        out_prefix = "staticFieldsOn%s"%suffix
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
            (out_var_names, dim_list) = get_hyperslab_name_and_dims(var_name)
            if has_time or combine_output:
                vtkFile = timeDependentFile
            else:
                vtkFile = timeIndependentFile
            for iHyperSlab in range(len(out_var_names)):
                if dim_list is not None:
                    dim_vals = dim_list[:,iHyperSlab]

                if (mesh_file is not None) and (var_name in mesh_file.variables):
                    nc_file = mesh_file
                else:
                    nc_file = time_series_file

                field_var = nc_file.variables[var_name]

                field_ndims = len(field_var.dimensions)

                field = np.zeros(blockDim,dtype=outType)

                try:
                    missing_val = field_var.missing_value
                except:
                    missing_val = -9999999790214767953607394487959552.000000

                # Build data
                for iBlock in np.arange(0, nBlocks):
                    blockStart = iBlock * blocking
                    blockEnd = min( (iBlock + 1) * blocking, blockDim )

                    if has_time:
                        assert(field_ndims != 1)

                        if field_ndims == 2:
                            field_block = field_var[local_time_indices[time_index], blockStart:blockEnd]
                        elif field_ndims == 3:
                            field_block = field_var[local_time_indices[time_index], blockStart:blockEnd, dim_vals[0]]
                        elif field_ndims == 4:
                            field_block = field_var[local_time_indices[time_index], blockStart:blockEnd, dim_vals[0], dim_vals[1]]
                        elif field_ndims == 5:
                            field_block = field_var[local_time_indices[time_index], blockStart:blockEnd, dim_vals[0], dim_vals[1], dim_vals[2]]
                    else:
                        if field_ndims == 1:
                            field_block = field_var[blockStart:blockEnd]
                        elif field_ndims == 2:
                            field_block = field_var[blockStart:blockEnd, dim_vals[0]]
                        elif field_ndims == 3:
                            field_block = field_var[blockStart:blockEnd, dim_vals[0], dim_vals[1]]
                        elif field_ndims == 4:
                            field_block = field_var[blockStart:blockEnd, dim_vals[0], dim_vals[1], dim_vals[2]]

                    # convert to the same type as field before masking with NaNs
                    field_block = np.array(field_block, dtype=field.dtype)
                    field_block[field_block == missing_val] = np.nan
                    field[blockStart:blockEnd] = field_block


                field = field[valid_mask]

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
    parser.add_argument("-b", "--blocking", dest="blocking", help="Size of blocks when reading MPAS file", metavar="BLK")
    parser.add_argument("-v", "--variable_list", dest="variable_list", help="List of variables to extract", metavar="VAR", required=True)
    parser.add_argument("-3", "--32bit", dest="output_32bit", help="If set, the vtk files will be written using 32bit floats.", action="store_true")
    parser.add_argument("-c", "--combine", dest="combine_output", help="If set, time-independent fields are written to each file along with time-dependent fields.", action="store_true")
    parser.add_argument("-a", "--append", dest="append", help="If set, only vtp files that do not already exist are written out.", action="store_true")
    parser.add_argument("-d", "--dim_list", dest="dimension_list", nargs="+", help="A list of dimensions and associated values.")
    args = parser.parse_args()

    if not args.output_32bit:
        use_32bit = False
    else:
        use_32bit = True

    if not args.blocking:
        args.blocking = int(10000)
    else:
        args.blocking = int(args.blocking)

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
    extra_dims = parse_extra_dims(args.dimension_list, time_series_file, mesh_file)
    (all_dim_vals, cellVars, vertexVars, edgeVars) = setup_dimension_values_and_sort_vars(
            time_series_file, mesh_file, args.variable_list, extra_dims)
    time_series_file.close()
    if(mesh_file is not None):
        mesh_file.close()

    summarize_extraction(args.mesh_filename, time_indices, cellVars, vertexVars, edgeVars)

    # Handle cell variables
    if len(cellVars) > 0:
        print " -- Extracting cell fields --"

        mesh_file = NetCDFFile(args.mesh_filename, 'r')

        # Build vertex list
        vertices = build_location_list_xyz( mesh_file, 'xVertex', 'yVertex', 'zVertex', use_32bit )

        # Build cell list
        (connectivity, offsets, valid_mask) = build_cell_lists( mesh_file, args.blocking )

        if not separate_mesh_file:
            mesh_file.close()
            mesh_file = None

        build_field_time_series( time_indices, time_file_names, mesh_file, args.blocking,
                                 all_dim_vals, 'nCells', cellVars, vertices, connectivity, offsets,
                                 valid_mask, use_32bit, args.combine_output, args.append )
        if separate_mesh_file:
            mesh_file.close()

        print ""
        del vertices
        del connectivity
        del offsets
        del valid_mask

    if len(vertexVars) > 0:
        print " -- Extracting vertex fields --"

        mesh_file = NetCDFFile(args.mesh_filename, 'r')

        # Build vertex list
        vertices = build_location_list_xyz( mesh_file, 'xCell', 'yCell', 'zCell', use_32bit )

        # Build cell list
        (connectivity, offsets, valid_mask) = build_dual_cell_lists( mesh_file, args.blocking )

        if not separate_mesh_file:
            mesh_file.close()
            mesh_file = None

        build_field_time_series( time_indices, time_file_names, mesh_file, args.blocking,
                                 all_dim_vals, 'nVertices', vertexVars, vertices, connectivity, offsets,
                                 valid_mask, use_32bit, args.combine_output, args.append )

        if separate_mesh_file:
            mesh_file.close()

        print ""
        del vertices
        del connectivity
        del offsets
        del valid_mask

    if len(edgeVars) > 0:
        print " -- Extracting edge fields --"

        mesh_file = NetCDFFile(args.mesh_filename, 'r')

        # Build vertex list
        verticesC = build_location_list_xyz( mesh_file, 'xCell', 'yCell', 'zCell', use_32bit )
        verticesV = build_location_list_xyz( mesh_file, 'xVertex', 'yVertex', 'zVertex', use_32bit )
        # compine the two into a single tuple
        vertices = (np.append(verticesC[0],verticesV[0]),
                    np.append(verticesC[1],verticesV[1]),
                    np.append(verticesC[2],verticesV[2]))
        del verticesC, verticesV

        # Build cell list
        (connectivity, offsets, valid_mask) = build_edge_cell_lists( mesh_file, args.blocking )

        if not separate_mesh_file:
            mesh_file.close()
            mesh_file = None

        build_field_time_series( time_indices, time_file_names, mesh_file, args.blocking,
                                 all_dim_vals, 'nEdges', edgeVars, vertices, connectivity, offsets,
                                 valid_mask, use_32bit, args.combine_output, args.append )

        if separate_mesh_file:
            mesh_file.close()


# vim: set expandtab:
