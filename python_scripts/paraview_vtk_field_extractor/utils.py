#!/usr/bin/env python
"""
Name: utils.py
Authors: Xylar Asay-Davis
Date: 06/20/2016

Utility library for various scripts used to extract vtk geometry from NetCDF
files on MPAS grids.
"""

from pyevtk.vtk import VtkFile, VtkPolyData

import sys, glob
import numpy

from netCDF4 import Dataset as NetCDFFile

try:
    from progressbar import ProgressBar, Percentage, Bar, ETA
    use_progress_bar = True
except:
    use_progress_bar = False

def is_valid_mesh_var(mesh_file, variable_name):
    if mesh_file is None:
        return False

    if variable_name not in mesh_file.variables:
        return False

    return 'Time' not in mesh_file.variables[variable_name].dimensions

def get_var(variable_name, mesh_file, time_series_file):
    if is_valid_mesh_var(mesh_file, variable_name):
        return mesh_file.variables[variable_name]
    else:
        return time_series_file.variables[variable_name]

# This function finds a list of NetCDF files containing time-dependent
# MPAS data and extracts the time indices in each file.  The routine
# insures that each time is unique.
def setup_time_indices(fn_pattern):#{{{
    # Build file list and time indices
    if ';' in fn_pattern:
        file_list = []
        for pattern in fn_pattern.split(';'):
            file_list.extend(glob.glob(pattern))
    else:
        file_list = glob.glob(fn_pattern)
    file_list.sort()

    local_indices = []
    file_names = []
    all_times = []

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
    for file_name in file_list:
        nc_file = NetCDFFile(file_name, 'r')


        try:
            xtime = nc_file.variables['xtime'][:,:]
            local_times = []
            for index in range(xtime.shape[0]):
              local_times.append(''.join(xtime[index,:]))

        except:
            local_times = ['0']

        if(len(local_times) == 0):
            local_times = ['0']
        nTime = len(local_times)

        for time_idx in range(nTime):
            if local_times[time_idx] not in all_times:
                local_indices.append(time_idx)
                file_names.append(file_name)
                all_times.append(local_times[time_idx])

        i_file = i_file + 1
        nc_file.close()
        if use_progress_bar:
            time_bar.update(i_file)

    if use_progress_bar:
        time_bar.finish()

    return (local_indices, file_names)
#}}}

# Parses the indices to be extracted along a given dimension.
# The index_string can be fomatted as follows:
#   <blank> -- no indices are to be extracted
#   n       -- the index n is to be extracted
#   m,n,p   -- the list of indices is to be extracted
#   m:n     -- all indices from m to n are to be extracted (including m but
#              excluding n, in the typical python indexing convention)
def parse_extra_dim(dim_name, index_string, time_series_file, mesh_file):#{{{
    if index_string == '':
        return numpy.zeros(0,int)

    if (mesh_file is not None) and (dim_name in mesh_file.dimensions):
        nc_file = mesh_file
    else:
        nc_file = time_series_file
    dim_size = len(nc_file.dimensions[dim_name])
    if ',' in index_string:
        indices = []
        for index in index_string.split(','):
            indices.append(int(index))
        indices = numpy.array(indices,int)
    elif ':' in index_string:
        index_list = index_string.split(':')
        if len(index_list) in [2,3]:
            if index_list[0] == '':
                first = 0
            else:
                first = int(index_list[0])
            if index_list[1] == '':
                last = dim_size
            else:
                last = int(index_list[1])
            if (len(index_list) == 2) or (index_list[2] == ''):
                step = 1
            else:
                step = int(index_list[2])
            indices = numpy.arange(first,last,step)
        else:
            print "Improperly formatted extra dimension:", dim_name, index_string
            return None
    else:
        indices = numpy.array([int(index_string)])

    valid = numpy.all(numpy.logical_and(indices >= 0,indices <= dim_size-1))
    if not valid:
        print "Index (or indices) out of bounds for extra dimension:", dim_name, index_string
        return None

    return indices

#}}}

# Parses a list of dimensions and corresponding indices separated by equals signs.
# Optionally, a max_index_count (typically 1) can be provided, indicating that
# indices beyond max_index_count-1 will be ignored in each dimension.
# Optionally, topo_dim contains the name of a dimension associated with the
# surface or bottom topography (e.g. nVertLevels for MPAS-Ocean)
# If too_dim is provided, topo_cell_indices_name can optionally be either
# a constant value for the index vertical index to the topography or
# the name of a field with dimension nCells that contains the vertical index of
# the topography.
def parse_extra_dims(dimension_list, time_series_file, mesh_file,
                     max_index_count=None):#{{{
    if not dimension_list:
        return {}

    extra_dims = {}
    for dim_item in dimension_list:
        (dimName,index_string) = dim_item.split('=')
        indices = parse_extra_dim(dimName, index_string, time_series_file, mesh_file)
        if indices is not None:
            if max_index_count is None or len(indices) <= max_index_count:
                extra_dims[dimName] = indices
            else:
                extra_dims[dimName] = indices[0:max_index_count]

    return extra_dims
#}}}


# Creates a list of variables names to be extracted.  Prompts for indices
# of any extra dimensions that were not specified on the command line.
# extra_dims should be a dictionary of indices along extra dimensions (as
# opposed to "basic" dimensions).  basic_dims is a list of dimension names
# that should be excluded from extra_dims.  include_dims is a list of
# possible dimensions, one of which must be in each vairable to be extracted
# (used in expanding command line placeholders "all", "allOnCells", etc.)

def setup_dimension_values_and_sort_vars(time_series_file, mesh_file, variable_list, extra_dims,
                                         basic_dims=['nCells', 'nEdges', 'nVertices', 'Time'],
                                         include_dims=['nCells', 'nEdges', 'nVertices']):#{{{

    def add_var(variables, variable_name, include_dims, exclude_dims=None):
        if variable_name in variable_names:
            return

        dims = variables[variable_name].dimensions
        supported = False
        for dim in include_dims:
            if dim in dims:
                supported = True
        if (exclude_dims is not None):
            for dim in exclude_dims:
                if dim in dims:
                    supported = False
        if supported:
            variable_names.append(variable_name)

    all_dim_vals = {}
    cellVars = []
    vertexVars = []
    edgeVars = []

    if variable_list == 'all':
        variable_names = []
        exclude_dims = ['Time']
        for variable_name in time_series_file.variables:
            add_var(time_series_file.variables, str(variable_name), include_dims, exclude_dims=None)
        if mesh_file is not None:
            for variable_name in mesh_file.variables:
                add_var(mesh_file.variables, str(variable_name), include_dims, exclude_dims)
    else:
        variable_names = variable_list.split(',')

    for suffix in ['Cells','Edges','Vertices']:
        include_dim = 'n%s'%suffix
        if ('allOn%s'%suffix in variable_names) and (include_dim in include_dims):
            variable_names.remove('allOn%s'%suffix)
            exclude_dims = ['Time']
            for variable_name in time_series_file.variables:
                add_var(time_series_file.variables, str(variable_name),
                        include_dims=[include_dim], exclude_dims=None)
            if mesh_file is not None:
                for variable_name in mesh_file.variables:
                    add_var(mesh_file.variables, str(variable_name),
                            include_dims=[include_dim], exclude_dims=exclude_dims)

    variable_names.sort()

    promptDimNames = []
    display_prompt = True
    for variable_name in variable_names:
        if is_valid_mesh_var(mesh_file, variable_name):
            nc_file = mesh_file
        else:
            nc_file = time_series_file
        field_dims = nc_file.variables[variable_name].dimensions
        for dim in field_dims:
            if ((dim in basic_dims) or (dim in extra_dims) or (dim in promptDimNames)):
                # this dimension has already been accounted for
                continue
            promptDimNames.append(str(dim))

            if display_prompt:
                print ""
                print "Need to define additional dimension values"
                display_prompt = False

            dim_size = len(nc_file.dimensions[dim])
            valid = False
            while not valid:
                print "Valid range for dimension %s between 0 and %d"%(dim, dim_size-1)
                index_string = raw_input("Enter a value for dimension %s: "%(dim))
                indices = parse_extra_dim(str(dim), index_string, time_series_file, mesh_file)
                valid = indices is not None
                if valid:
                    extra_dims[str(dim)] = indices
                else:
                    print " -- Invalid value, please re-enter --"

    empty_dims = []
    for dim in extra_dims:
        if len(extra_dims[dim]) == 0:
            empty_dims.append(dim)

    for variable_name in variable_names:

        field_dims = get_var(variable_name, mesh_file, time_series_file).dimensions
        skip = False
        for dim in field_dims:
            if dim in empty_dims:
                skip = True
                break
        if skip:
            continue

        # Setting dimension values:
        dim_vals = []
        indices = []
        for dim in field_dims:
            if dim not in basic_dims:
                indices.append(extra_dims[str(dim)])
        if len(indices) == 0:
            dim_vals = None
        else:
            # expand full list of indices all extra dimensions
            dim_vals = numpy.array(numpy.meshgrid(indices,indexing='ij'))

        if "nCells" in field_dims:
            cellVars.append(variable_name)
        elif "nVertices" in field_dims:
            vertexVars.append(variable_name)
        elif "nEdges" in field_dims:
            edgeVars.append(variable_name)

        all_dim_vals[variable_name] = dim_vals
        del dim_vals

    return (all_dim_vals, cellVars, vertexVars, edgeVars)
#}}}

# Print a summary of the time levels, mesh file, transects file (optional)
# and variables to be extracted.
def summarize_extraction(mesh_file, time_indices, cellVars, vertexVars, edgeVars,
                         transects_file=None):#{{{
    print ""
    print "Extracting a total of %d time levels."%(len(time_indices))
    print "Using file '%s' as the mesh file for this extraction."%(mesh_file)
    if transects_file is not None:
        print "Using file '%s' as the transects file."%(transects_file)
    print ""
    print ""
    print "The following variables will be extracted from the input file(s)."
    print ""

    if len(cellVars) > 0:
        print "   Variables with 'nCells' as a dimension:"
        for variable_name in cellVars:
            print "      name: %s"%(variable_name)

    if len(vertexVars) > 0:
        print "   Variables with 'nVertices' as a dimension:"
        for variable_name in vertexVars:
            print "      name: %s"%(variable_name)

    if len(edgeVars) > 0:
        print "   Variables with 'nEdges' as adimension:"
        for variable_name in edgeVars:
            print "      name: %s"%(variable_name)

    print ""
#}}}

def write_pvd_header(path, prefix):#{{{
    pvd_file = open('%s/%s.pvd'%(path, prefix), 'w')
    pvd_file.write('<?xml version="1.0"?>\n')
    pvd_file.write('<VTKFile type="Collection" version="0.1"\n')
    pvd_file.write('\tbyte_order="LittleEndian"\n')
    pvd_file.write('\tcompressor="vtkZLibDataCompressor">\n')
    return pvd_file
#}}}

def get_hyperslab_name_and_dims(var_name, all_dim_vals):#{{{
    extra_dim_vals = all_dim_vals[var_name]
    if(extra_dim_vals is None):
        return ([var_name],None)
    if(extra_dim_vals.size == 0):
        return ([],None)
    out_var_names = []
    maxval = numpy.amax(extra_dim_vals,axis=1)
    pad = numpy.ones(maxval.shape,int)
    mask = maxval > 0
    pad[mask] = numpy.array(numpy.floor(numpy.log10(maxval[mask])),int)+1
    for iHyperSlab in range(extra_dim_vals.shape[1]):
        out_var_name = var_name
        for iVal in range(extra_dim_vals.shape[0]):
            val = extra_dim_vals[iVal,iHyperSlab]
            template = '%%s_%%0%dd'%pad[iVal]
            out_var_name = template%(out_var_name,val)
            out_var_names.append(out_var_name)
    return (out_var_names, extra_dim_vals)
#}}}

def write_vtp_header(path, prefix, active_var_index, var_indices, variable_list, all_dim_vals,
                     vertices, connectivity, offsets, nPoints, nPolygons, outType,
                     cellData=True, pointData=False):#{{{
    vtkFile = VtkFile("%s/%s"%(path,prefix), VtkPolyData)
    vtkFile.openElement(vtkFile.ftype.name)
    vtkFile.openPiece(npoints=nPoints,npolys=nPolygons)

    vtkFile.openElement("Points")
    vtkFile.addData("points", vertices)
    vtkFile.closeElement("Points")

    vtkFile.openElement("Polys")
    vtkFile.addData("connectivity", connectivity)
    vtkFile.addData("offsets", offsets)
    vtkFile.closeElement("Polys")

    if(cellData):
        vtkFile.openData("Cell", scalars = variable_list[active_var_index])
        for iVar in var_indices:
            var_name = variable_list[iVar]
            (out_var_names, dim_list) = get_hyperslab_name_and_dims(var_name, all_dim_vals)
            for out_var_name in out_var_names:
                vtkFile.addHeader(out_var_name, outType, nPolygons, 1)
        vtkFile.closeData("Cell")
    if(pointData):
        vtkFile.openData("Point", scalars = variable_list[active_var_index])
        for iVar in var_indices:
            var_name = variable_list[iVar]
            (out_var_names, dim_list) = get_hyperslab_name_and_dims(var_name, all_dim_vals)
            for out_var_name in out_var_names:
                vtkFile.addHeader(out_var_name, outType, nPoints, 1)
        vtkFile.closeData("Point")

    vtkFile.closePiece()
    vtkFile.closeElement(vtkFile.ftype.name)

    vtkFile.appendData(vertices)
    vtkFile.appendData(connectivity)
    vtkFile.appendData(offsets)

    return vtkFile
#}}}


def build_cell_lists( nc_file, blocking ):#{{{
    nCells = len(nc_file.dimensions['nCells'])

    nEdgesOnCell_var = nc_file.variables['nEdgesOnCell']
    verticesOnCell_var = nc_file.variables['verticesOnCell']

    offsets = numpy.cumsum(nEdgesOnCell_var[:], dtype=int)
    connectivity = numpy.zeros(offsets[-1], dtype=int)
    valid_mask = numpy.ones(nCells,bool)

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

    valid_mask = numpy.all(cellsOnVertex >= 0, axis=1)
    connectivity = cellsOnVertex[valid_mask,:].ravel()
    validCount = numpy.count_nonzero(valid_mask)
    offsets = vertexDegree*numpy.arange(1,validCount+1)

    return (connectivity, offsets, valid_mask)

#}}}

def build_edge_cell_lists( nc_file, blocking ):#{{{
    nCells = len(nc_file.dimensions['nCells'])
    nEdges = len(nc_file.dimensions['nEdges'])

    cellsOnEdge_var = nc_file.variables['cellsOnEdge']
    verticesOnEdge_var = nc_file.variables['verticesOnEdge']

    valid_mask = numpy.zeros(nEdges,bool)

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
    for iBlock in numpy.arange(0, nBlocks):
        blockStart = iBlock * blocking
        blockEnd = min( (iBlock + 1) * blocking, nEdges )
        blockCount = blockEnd - blockStart

        verticesOnEdge = verticesOnEdge_var[blockStart:blockEnd,:] - 1
        cellsOnEdge = cellsOnEdge_var[blockStart:blockEnd,:] - 1

        vertices = numpy.zeros((blockCount,4))
        vertices[:,0] = cellsOnEdge[:,0]
        vertices[:,1] = verticesOnEdge[:,0]
        vertices[:,2] = cellsOnEdge[:,1]
        vertices[:,3] = verticesOnEdge[:,1]
        valid = vertices >= 0
        vertices[:,1] += nCells
        vertices[:,3] += nCells
        validCount = numpy.sum(numpy.array(valid,int),axis=1)

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

    connectivity = numpy.array(connectivity, dtype=int)
    offsets = numpy.array(offsets, dtype=int)

    return (connectivity, offsets, valid_mask)

#}}}

def build_location_list_xyz( nc_file, xName, yName, zName, output_32bit ):#{{{

    X = nc_file.variables[xName][:]
    Y = nc_file.variables[yName][:]
    Z = nc_file.variables[zName][:]
    if output_32bit:
        X = numpy.array(X,'f4')
        Y = numpy.array(Y,'f4')
        Z = numpy.array(Z,'f4')
    return (X,Y,Z)

#}}}


def get_field_sign(field_name):
    if field_name[0] == '-':
        field_name = field_name[1:]
        sign = -1
    else:
        sign = 1

    return (field_name, sign)

def read_field(field_var, extra_dim_vals, time_index, cell_indices,
               outType, sign=1, vert_dim_name=None, nLevels=None):#{{{

    def read_field_with_dims(field_var, dim_vals, field):#{{{
        outDims = len(field.shape)
        inDims = len(dim_vals)
        if outDims <= 0 or outDims > 2:
            print 'reading field into variable with dimension', outDims, 'not supported.'
            exit(1)
        if inDims <= 0 or inDims > 5:
            print 'reading field with dimension', inDims, 'not supported.'
            exit(1)

        if outDims == 1:
            if inDims == 1:
                field[:] = field_var[dim_vals[0]]
            elif inDims == 2:
                field[:] = field_var[dim_vals[0], dim_vals[1]]
            elif inDims == 3:
                field[:] = field_var[dim_vals[0], dim_vals[1], dim_vals[2]]
            elif inDims == 4:
                field[:] = field_var[dim_vals[0], dim_vals[1], dim_vals[2], dim_vals[3]]
            elif inDims == 5:
                field[:] = field_var[dim_vals[0], dim_vals[1], dim_vals[2], dim_vals[3], dim_vals[4]]
        elif outDims == 2:
            if inDims == 1:
                field[:,:] = field_var[dim_vals[0]]
            elif inDims == 2:
                field[:,:] = field_var[dim_vals[0], dim_vals[1]]
            elif inDims == 3:
                field[:,:] = field_var[dim_vals[0], dim_vals[1], dim_vals[2]]
            elif inDims == 4:
                field[:,:] = field_var[dim_vals[0], dim_vals[1], dim_vals[2], dim_vals[3]]
            elif inDims == 5:
                field[:,:] = field_var[dim_vals[0], dim_vals[1], dim_vals[2], dim_vals[3], dim_vals[4]]
#}}}


    dim_vals = []
    extraDimIndex = 0
    for dim in field_var.dimensions:
        if dim == 'Time':
            dim_vals.append(time_index)
        elif dim == 'nCells':
            dim_vals.append(cell_indices)
        elif dim == vert_dim_name:
            dim_vals.append(numpy.arange(nLevels))
        else:
            dim_vals.append(extra_dim_vals[extraDimIndex])
            extraDimIndex += 1

    if vert_dim_name in field_var.dimensions:
        field = numpy.zeros((len(cell_indices),nLevels), dtype=outType)

    else:
        field = numpy.zeros((len(cell_indices)), dtype=outType)

    read_field_with_dims(field_var, dim_vals, field)

    return sign*field
#}}}


def compute_zInterface(minLevelCell, maxLevelCell, layerThicknessCell,
                           zMinCell, zMaxCell, dtype, cellsOnEdge=None):#{{{

    (nCells,nLevels) = layerThicknessCell.shape

    cellMask = numpy.ones((nCells,nLevels), bool)
    for iLevel in range(nLevels):
        if minLevelCell is not None:
            cellMask[:,iLevel] = numpy.logical_and(cellMask[:,iLevel], iLevel >= minLevelCell)
        if maxLevelCell is not None:
            cellMask[:,iLevel] = numpy.logical_and(cellMask[:,iLevel], iLevel <= maxLevelCell)

    zInterfaceCell = numpy.zeros((nCells,nLevels+1),dtype=dtype)
    for iLevel in range(nLevels):
        zInterfaceCell[:,iLevel+1] = (zInterfaceCell[:,iLevel]
            + cellMask[:,iLevel]*layerThicknessCell[:,iLevel])

    if zMinCell is not None:
        minLevel = minLevelCell.copy()
        minLevel[minLevel < 0] = nLevels-1
        zOffsetCell = zMinCell - zInterfaceCell[numpy.arange(0,nCells),minLevel]
    else:
        zOffsetCell = zMaxCell - zInterfaceCell[numpy.arange(0,nCells),maxLevelCell+1]

    for iLevel in range(nLevels+1):
        zInterfaceCell[:,iLevel] += zOffsetCell

    if cellsOnEdge is None:
         return zInterfaceCell
    else:
        nEdges = cellsOnEdge.shape[0]
        zInterfaceEdge = numpy.zeros((nEdges,nLevels+1),dtype=dtype)

        # Get a list of valid cells on edges and a mask of which are valid
        cellsOnEdgeMask = numpy.logical_and(cellsOnEdge >= 0, cellsOnEdge < nCells)
        cellIndicesOnEdge = []
        cellIndicesOnEdge.append(cellsOnEdge[cellsOnEdgeMask[:,0],0])
        cellIndicesOnEdge.append(cellsOnEdge[cellsOnEdgeMask[:,1],1])

        for iLevel in range(nLevels):
            edgeMask = numpy.zeros(nEdges, bool)
            layerThicknessEdge = numpy.zeros(nEdges, float)
            denom = numpy.zeros(nEdges, float)
            for index in range(2):
                mask = cellsOnEdgeMask[:,index]
                cellIndices = cellIndicesOnEdge[index]
                cellMaskLocal = cellMask[cellIndices,iLevel]

                edgeMask[mask] = numpy.logical_or(edgeMask[mask], cellMaskLocal)

                layerThicknessEdge[mask] += cellMaskLocal*layerThicknessCell[cellIndices,iLevel]
                denom[mask] += 1.0*cellMaskLocal

            layerThicknessEdge[edgeMask] /= denom[edgeMask]

            zInterfaceEdge[:,iLevel+1] = (zInterfaceEdge[:,iLevel]
                                       + edgeMask*layerThicknessEdge)

        if zMinCell is not None:
            refLevelEdge = numpy.zeros(nEdges, int)
            for index in range(2):
                mask = cellsOnEdgeMask[:,index]
                cellIndices = cellIndicesOnEdge[index]
                refLevelEdge[mask] = numpy.maximum(refLevelEdge[mask], minLevel[cellIndices])
        else:
            refLevelEdge = (nLevels-1)*numpy.ones(nEdges, int)
            for index in range(2):
                mask = cellsOnEdgeMask[:,index]
                cellIndices = cellIndicesOnEdge[index]
                refLevelEdge[mask] = numpy.minimum(refLevelEdge[mask], maxLevelCell[cellIndices]+1)


        zOffsetEdge = numpy.zeros(nEdges, float)
        # add the average of zInterfaceCell at each adjacent cell
        denom = numpy.zeros(nEdges, float)
        for index in range(2):
            mask = cellsOnEdgeMask[:,index]
            cellIndices = cellIndicesOnEdge[index]
            zOffsetEdge[mask] += zInterfaceCell[cellIndices,refLevelEdge[mask]]
            denom[mask] += 1.0

        mask = denom > 0.
        zOffsetEdge[mask] /= denom[mask]

        # subtract the depth of zInterfaceEdge at the level of the bottom
        zOffsetEdge -= zInterfaceEdge[numpy.arange(nEdges),refLevelEdge]

        for iLevel in range(nLevels+1):
            zInterfaceEdge[:,iLevel] += zOffsetEdge

        return (zInterfaceCell, zInterfaceEdge)

#}}}


# vim: set expandtab:
