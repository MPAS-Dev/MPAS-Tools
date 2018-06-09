#!/usr/bin/env python
"""
Name: utils.py
Authors: Xylar Asay-Davis

Utility library for various scripts used to extract vtk geometry from NetCDF
files on MPAS grids.
"""

try:
    from evtk.vtk import VtkFile, VtkPolyData
except ImportError:
    from pyevtk.vtk import VtkFile, VtkPolyData

import sys, glob
import numpy

from netCDF4 import Dataset as NetCDFFile

try:
    from progressbar import ProgressBar, Percentage, Bar, ETA
    use_progress_bar = True
except:
    use_progress_bar = False


def open_netcdf(file_name):
    nc_file = NetCDFFile(file_name, 'r')
    # turn off auto mask (if applicable)
    try:
        nc_file.set_auto_mask(False)
    except AttributeError:
        pass
    return nc_file


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


def setup_time_indices(fn_pattern, xtimeName):  # {{{
    """
    This function finds a list of NetCDF files containing time-dependent
    MPAS data and extracts the time indices in each file.  The routine
    insures that each time is unique.
    """
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
        widgets = ['Build time indices: ', Percentage(), ' ', Bar(), ' ',
                   ETA()]
        time_bar = ProgressBar(widgets=widgets, maxval=len(file_list)).start()
    else:
        print "Build time indices..."

    i_file = 0
    for file_name in file_list:
        try:
            nc_file = open_netcdf(file_name)
        except IOError:
            print "Warning: could not open {}".format(file_name)
            continue

        if 'Time' not in nc_file.dimensions or xtimeName is None:
            local_times = ['0']
        else:
            if xtimeName not in nc_file.variables:
                raise ValueError("xtime variable name {} not found in "
                                 "{}".format(xtimeName, file_name))
            local_times = []
            xtime = nc_file.variables[xtimeName][:, :]
            for index in range(xtime.shape[0]):
                local_times.append(''.join(xtime[index, :]))

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

    return (local_indices, file_names)  # }}}


# Parses the indices to be extracted along a given dimension.
# The index_string can be fomatted as follows:
#   <blank> -- no indices are to be extracted
#   n       -- the index n is to be extracted
#   m,n,p   -- the list of indices is to be extracted
#   m:n     -- all indices from m to n are to be extracted (including m but
#              excluding n, in the typical python indexing convention)
def parse_extra_dim(dim_name, index_string, time_series_file, mesh_file):#{{{
    """
    Parses the indices to be extracted along a given dimension.
    The index_string can be fomatted as follows:
      <blank> -- no indices are to be extracted
      n       -- the index n is to be extracted
      m,n,p   -- the list of indices is to be extracted
      m:n     -- all indices from m to n are to be extracted (including m but
                 excluding n, in the typical python indexing convention)
      m:n:s   -- all indices from m to n are to be extracted (including m but
                 excluding n, in the typical python indexing convention) with
                 stride s between indices

    Parameters
    ----------
    dim_name : str
        The name of the dimension to be indexed

    index_string : str
        An index string indicating with indices are to be extracted

    time_series_file : str
        The name of a time series file that can be used to determine the size
        of the dimension if ``mesh_file=None``.

    mesh_file : str
        The name of a mesh file that can be used to determine the size
        of the dimension, or ``None`` if the time series file should be used

    Returns
    -------
    indices : list of str
        Indices into the given dimension formatted as zero-padded strings
        (if indices are numerical, as opposed to the name of an index variable)
    """
    if index_string == '':
        return []

    if (mesh_file is not None) and (dim_name in mesh_file.dimensions):
        nc_file = mesh_file
    else:
        nc_file = time_series_file

    dim_size = len(nc_file.dimensions[dim_name])

    indices, numerical_indices = parse_index_string(index_string, dim_size)

    # zero-pad integer indices
    if len(numerical_indices) > 0:
        max_index = numpy.amax(numerical_indices)
        pad = int(numpy.log10(max(max_index,1)))+1
        template = '%%0%dd'%pad
        for i in range(len(indices)):
            try:
                val = int(indices[i])
            except ValueError:
                continue
            indices[i] = template%(val)

    return indices

#}}}


def parse_time_indices(index_string, time_indices, time_file_names):  # {{{
    """
    Parses the indices to be extracted along the Time dimension.
    The index_string can be fomatted as follows:
      <blank> -- no indices are to be extracted
      n       -- the index n is to be extracted
      m,n,p   -- the list of indices is to be extracted
      m:n     -- all indices from m to n are to be extracted (including m but
                 excluding n, in the typical python indexing convention)
      m:n:s   -- all indices from m to n are to be extracted (including m but
                 excluding n, in the typical python indexing convention) with
                 stride s between indices

    Parameters
    ----------
    index_string : str
        An index string indicating with indices are to be extracted

    time_indices : list of int
        The local time indices in each input NetCDF file

    time_file_names : list of str
        The name of the file associated with each time index

    Returns
    -------
    time_indices : list of int
        The local time indices in each input NetCDF file after reindexing

    time_file_names : list of str
        The name of the file associated with each time index after reindexing

    """
    dim_size = len(time_indices)
    indices, numerical_indices = parse_index_string(index_string, dim_size)

    time_indices = [time_indices[index] for index in numerical_indices]
    time_file_names = [time_file_names[index] for index in numerical_indices]
    return time_indices, time_file_names  # }}}


def parse_index_string(index_string, dim_size):  # {{{
    """
    Parses an index string into a list of indices.
    The index_string can be fomatted as follows:
      <blank> -- no indices are to be extracted
      n       -- the index n is to be extracted
      m,n,p   -- the list of indices is to be extracted
      m:n     -- all indices from m to n are to be extracted (including m but
                 excluding n, in the typical python indexing convention)
      m:n:s   -- all indices from m to n are to be extracted (including m but
                 excluding n, in the typical python indexing convention) with
                 stride s between indices

    Parameters
    ----------
    index_string : str
        An index string indicating with indices are to be extracted

    dim_size : int
        The size of the dimension

    Returns
    -------
    indices : list of int
        The indices corresponding to the given index string.
    """
    if index_string == '':
        return [], []

    if ',' in index_string:
        indices = [index for index in index_string.split(',')]
    elif ':' in index_string:
        index_list = index_string.split(':')
        if len(index_list) not in [2, 3]:
            raise ValueError("Improperly formatted index string: "
                             "{}".format(index_string))
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
        indices = [str(index) for index in numpy.arange(first, last, step)]
    else:
        indices = [index_string]

    numerical_indices = []
    for index in indices:
        try:
            val = int(index)
        except ValueError:
            continue

        if val < 0 or val >= dim_size:
            raise ValueError("Index (or indices) out of bounds 0 <= "
                             "index < {}: {}".format(dim_size, index_string))

        numerical_indices.append(val)

    return indices, numerical_indices  # }}}


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
        indices = []
        for dim in field_dims:
            if dim not in ['Time', 'nCells', 'nEdges', 'nVertices']:
                indices.append(extra_dims[dim])
        if len(indices) == 0:
            dim_vals = None
        elif len(indices) == 1:
            dim_vals = []
            for index0 in indices[0]:
               dim_vals.append([index0])
        elif len(indices) == 2:
            dim_vals = []
            for index0 in indices[0]:
                for index1 in indices[1]:
                    dim_vals.append([index0,index1])
        elif len(indices) == 3:
            dim_vals = []
            for index0 in indices[0]:
                for index1 in indices[1]:
                    for index2 in indices[2]:
                        dim_vals.append([index0,index1,index2])
        else:
           print "variable %s has too many extra dimensions and will be skipped."%variable_name
           continue

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


def write_pvd_header(path, prefix):  # {{{
    pvd_file = open('{}/{}.pvd'.format(path, prefix), 'w')
    pvd_file.write('<?xml version="1.0"?>\n')
    pvd_file.write('<VTKFile type="Collection" version="0.1"\n')
    pvd_file.write('\tbyte_order="LittleEndian"\n')
    pvd_file.write('\tcompressor="vtkZLibDataCompressor">\n')
    return pvd_file  # }}}


def get_hyperslab_name_and_dims(var_name, extra_dim_vals):#{{{
    if(extra_dim_vals is None):
        return ([var_name],None)
    if(len(extra_dim_vals) == 0):
        return ([],None)
    out_var_names = []
    for hyper_slab in extra_dim_vals:
        pieces = [var_name]
        pieces.extend(hyper_slab)
        out_var_names.append('_'.join(pieces))
    return (out_var_names, extra_dim_vals)
#}}}


def write_vtp_header(path, prefix, active_var_index, var_indices,
                     variable_list, all_dim_vals, vertices, connectivity,
                     offsets, nPoints, nPolygons, outType, cellData=True,
                     pointData=False, xtime=None):  # {{{
    vtkFile = VtkFile("{}/{}".format(path, prefix), VtkPolyData)

    if xtime is not None:
        vtkFile.openElement("metadata")
        vtkFile.openElement("xtime")
        vtkFile.xml.addText(xtime)
        vtkFile.closeElement("xtime")
        vtkFile.closeElement("metadata")

    vtkFile.openElement(vtkFile.ftype.name)
    vtkFile.openPiece(npoints=nPoints, npolys=nPolygons)

    vtkFile.openElement("Points")
    vtkFile.addData("points", vertices)
    vtkFile.closeElement("Points")

    vtkFile.openElement("Polys")
    vtkFile.addData("connectivity", connectivity)
    vtkFile.addData("offsets", offsets)
    vtkFile.closeElement("Polys")

    if(cellData):
        vtkFile.openData("Cell", scalars=variable_list[active_var_index])
        for iVar in var_indices:
            var_name = variable_list[iVar]
            (out_var_names, dim_list) = \
                get_hyperslab_name_and_dims(var_name, all_dim_vals[var_name])
            for out_var_name in out_var_names:
                vtkFile.addHeader(out_var_name, outType, nPolygons, 1)
        vtkFile.closeData("Cell")
    if(pointData):
        vtkFile.openData("Point", scalars=variable_list[active_var_index])
        for iVar in var_indices:
            var_name = variable_list[iVar]
            (out_var_names, dim_list) = \
                get_hyperslab_name_and_dims(var_name, all_dim_vals[var_name])
            for out_var_name in out_var_names:
                vtkFile.addHeader(out_var_name, outType, nPoints, 1)
        vtkFile.closeData("Point")

    vtkFile.closePiece()
    vtkFile.closeElement(vtkFile.ftype.name)

    vtkFile.appendData(vertices)
    vtkFile.appendData(connectivity)
    vtkFile.appendData(offsets)

    return vtkFile  # }}}


def build_cell_geom_lists(nc_file, output_32bit, lonlat):  # {{{

    print "Build geometry for fields on cells..."

    vertices = _build_location_list_xyz(nc_file, 'Vertex', output_32bit,
                                        lonlat)

    if lonlat:
        lonCell = numpy.rad2deg(nc_file.variables['lonCell'][:])

    nCells = len(nc_file.dimensions['nCells'])

    nEdgesOnCell = nc_file.variables['nEdgesOnCell'][:]
    verticesOnCell = nc_file.variables['verticesOnCell'][:, :] - 1

    validVertices = numpy.zeros(verticesOnCell.shape, bool)
    for vIndex in range(validVertices.shape[1]):
        validVertices[:, vIndex] = vIndex < nEdgesOnCell

    if lonlat:
        vertices, verticesOnCell = _fix_lon_lat_vertices(vertices,
                                                         verticesOnCell,
                                                         validVertices,
                                                         lonCell)

    connectivity = verticesOnCell[validVertices]
    offsets = numpy.cumsum(nEdgesOnCell, dtype=int)
    valid_mask = numpy.ones(nCells, bool)

    return (vertices, connectivity, offsets, valid_mask)  # }}}


def build_vertex_geom_lists(nc_file, output_32bit, lonlat):  # {{{
    print "Build geometry for fields on vertices...."

    vertices = _build_location_list_xyz(nc_file, 'Cell', output_32bit, lonlat)

    if lonlat:
        lonVertex = numpy.rad2deg(nc_file.variables['lonVertex'][:])

    vertexDegree = len(nc_file.dimensions['vertexDegree'])

    cellsOnVertex = nc_file.variables['cellsOnVertex'][:, :] - 1

    valid_mask = numpy.all(cellsOnVertex >= 0, axis=1)

    cellsOnVertex = cellsOnVertex[valid_mask, :]

    if lonlat:
        # all remaining entries in cellsOnVertex are valid
        validVertices = numpy.ones(cellsOnVertex.shape, bool)
        vertices, cellsOnVertex = _fix_lon_lat_vertices(vertices,
                                                        cellsOnVertex,
                                                        validVertices,
                                                        lonVertex[valid_mask])

    connectivity = cellsOnVertex.ravel()
    validCount = cellsOnVertex.shape[0]
    offsets = vertexDegree*numpy.arange(1, validCount+1)

    return (vertices, connectivity, offsets, valid_mask)  # }}}


def build_edge_geom_lists(nc_file, output_32bit, lonlat):  # {{{
    (xCell, yCell, zCell) = \
        _build_location_list_xyz(nc_file, 'Cell', output_32bit, lonlat)
    (xVertex, yVertex, zVertex) = \
        _build_location_list_xyz(nc_file, 'Vertex', output_32bit, lonlat)

    vertices = (numpy.append(xCell, xVertex),
                numpy.append(yCell, yVertex),
                numpy.append(zCell, zVertex))

    if lonlat:
        lonEdge = numpy.rad2deg(nc_file.variables['lonEdge'][:])

    nEdges = len(nc_file.dimensions['nEdges'])
    nCells = len(nc_file.dimensions['nCells'])

    cellsOnEdge = nc_file.variables['cellsOnEdge'][:] - 1
    verticesOnEdge = nc_file.variables['verticesOnEdge'][:] - 1

    vertsOnCell = -1*numpy.ones((nEdges, 4), int)
    vertsOnCell[:, 0] = cellsOnEdge[:, 0]
    vertsOnCell[:, 1] = verticesOnEdge[:, 0]
    vertsOnCell[:, 2] = cellsOnEdge[:, 1]
    vertsOnCell[:, 3] = verticesOnEdge[:, 1]

    validVerts = vertsOnCell >= 0

    vertsOnCell[:, 1] += nCells
    vertsOnCell[:, 3] += nCells

    validCount = numpy.sum(numpy.array(validVerts, int), axis=1)

    # exclude any isolated points or lines -- we need only triangles or quads
    valid_mask = validCount >= 3

    vertsOnCell = vertsOnCell[valid_mask, :]
    validVerts = validVerts[valid_mask, :]

    if lonlat:
        vertices, cellsOnVertex = _fix_lon_lat_vertices(vertices,
                                                        vertsOnCell,
                                                        validVerts,
                                                        lonEdge[valid_mask])

    connectivity = vertsOnCell[validVerts]
    validCount = numpy.sum(numpy.array(validVerts, int), axis=1)
    offsets = numpy.cumsum(validCount, dtype=int)

    return (vertices, connectivity, offsets, valid_mask)  # }}}


def get_field_sign(field_name):
    if field_name[0] == '-':
        field_name = field_name[1:]
        sign = -1
    else:
        sign = 1

    return (field_name, sign)

def read_field(var_name, mesh_file, time_series_file, extra_dim_vals, time_index,
               block_indices, outType, sign=1):#{{{

    def read_field_with_dims(field_var, dim_vals, temp_shape, outType, index_arrays):#{{{
        temp_field = numpy.zeros(temp_shape, dtype=outType)
        inDims = len(dim_vals)
        if inDims <= 0 or inDims > 5:
            print 'reading field %s with %s dimensions not supported.'%(var_name, inDims)
            sys.exit(1)

        if inDims == 1:
            temp_field = numpy.array(field_var[dim_vals[0]], dtype=outType)
        elif inDims == 2:
            temp_field = numpy.array(field_var[dim_vals[0], dim_vals[1]],
                                     dtype=outType)
        elif inDims == 3:
            temp_field = numpy.array(field_var[dim_vals[0], dim_vals[1],
                                               dim_vals[2]],
                                     dtype=outType)
        elif inDims == 4:
            temp_field = numpy.array(field_var[dim_vals[0], dim_vals[1],
                                               dim_vals[2], dim_vals[3]],
                                     dtype=outType)
        elif inDims == 5:
            temp_field = numpy.array(field_var[dim_vals[0], dim_vals[1],
                                               dim_vals[2], dim_vals[3],
                                               dim_vals[4]],
                                     dtype=outType)

        outDims = len(temp_field.shape)

        if outDims <= 0 or outDims > 4:
            print 'something went wrong reading field %s, resulting in a temp array with %s dimensions.'%(var_name, outDims)
            sys.exit(1)
        block_indices = numpy.arange(temp_field.shape[0])
        if outDims == 1:
            field = temp_field
        elif outDims == 2:
            field = temp_field[block_indices,index_arrays[0]]
        elif outDims == 3:
            field = temp_field[block_indices,index_arrays[0],index_arrays[1]]
        elif outDims == 4:
            field = temp_field[block_indices,index_arrays[0],index_arrays[1],index_arrays[2]]

        return field


#}}}

    field_var = get_var(var_name, mesh_file, time_series_file)
    try:
        missing_val = field_var.missing_value
    except:
        missing_val = -9999999790214767953607394487959552.000000

    dim_vals = []
    extra_dim_index = 0
    shape = field_var.shape
    temp_shape = ()

    index_arrays = []

    for i in range(field_var.ndim):
        dim =  field_var.dimensions[i]
        if dim == 'Time':
            dim_vals.append(time_index)
        elif dim in ['nCells', 'nEdges', 'nVertices']:
            dim_vals.append(block_indices)
            temp_shape = temp_shape + (len(block_indices),)
        else:
            extra_dim_val = extra_dim_vals[extra_dim_index]
            try:
                index = int(extra_dim_val)
                dim_vals.append(index)
            except ValueError:
                # we have an index array so we need to read the whole range
                # first and then sample the result
                dim_vals.append(numpy.arange(shape[i]))
                temp_shape = temp_shape + (shape[i],)

                index_array_var = get_var(extra_dim_val, mesh_file, time_series_file)

                # read the appropriate indices from the index_array_var
                index_array = numpy.maximum(0,numpy.minimum(shape[i]-1, index_array_var[block_indices]-1))

                index_arrays.append(index_array)

            extra_dim_index += 1


    field = read_field_with_dims(field_var, dim_vals, temp_shape, outType,
                                 index_arrays)

    field[field == missing_val] = numpy.nan

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


def _build_location_list_xyz(nc_file, suffix, output_32bit, lonlat):  # {{{

    if lonlat:
        X = numpy.rad2deg(nc_file.variables['lon{}'.format(suffix)][:])
        Y = numpy.rad2deg(nc_file.variables['lat{}'.format(suffix)][:])
        Z = numpy.zeros(X.shape)
    else:
        X = nc_file.variables['x{}'.format(suffix)][:]
        Y = nc_file.variables['y{}'.format(suffix)][:]
        Z = nc_file.variables['z{}'.format(suffix)][:]
    if output_32bit:
        X = numpy.array(X, 'f4')
        Y = numpy.array(Y, 'f4')
        Z = numpy.array(Z, 'f4')
    return (X, Y, Z)

# }}}


def _fix_lon_lat_vertices(vertices, verticesOnCell, validVertices,
                          lonCell):  # {{{

    nCells = verticesOnCell.shape[0]
    nVertices = len(vertices[0])

    xVertex = vertices[0]
    xVertex = vertices[0]
    xDiff = xVertex[verticesOnCell] - lonCell.reshape(nCells, 1)

    # which cells have vertices that are out of range?
    outOfRange = numpy.logical_and(validVertices,
                                   numpy.logical_or(xDiff >= 180.,
                                                    xDiff < -180.))

    cellsOutOfRange = numpy.any(outOfRange, axis=1)

    valid = validVertices[cellsOutOfRange, :]

    verticesToChange = numpy.zeros(verticesOnCell.shape, bool)
    verticesToChange[cellsOutOfRange, :] = valid

    xDiff = xDiff[cellsOutOfRange, :][valid]
    voc = verticesOnCell[cellsOutOfRange, :][valid]
    nVerticesToAdd = numpy.count_nonzero(valid)
    verticesToAdd = numpy.arange(nVerticesToAdd) + nVertices
    xv = xVertex[voc]
    verticesOnCell[verticesToChange] = verticesToAdd

    # where the lon. difference between the cell center and the vertex
    # is out of range (presumably because of the periodic boundary),
    # move it to be within 180 degrees.
    mask = xDiff >= 180.
    xv[mask] -= 360.
    mask = xDiff < -180.
    xv[mask] += 360.

    vertices = (numpy.append(xVertex, xv),
                numpy.append(vertices[1], vertices[1][voc]),
                numpy.append(vertices[2], vertices[2][voc]))

    return vertices, verticesOnCell  # }}}

# vim: set expandtab:
