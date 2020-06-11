"""
Extract a field from a time series of NetCDF files as VTK files for plotting in
paraview.

It can extract a field across multiple files by passing in a regular expression
for the filename patter. As an example, one can run the script using:

    ./paraview_vtk_field_extractor.py -v areaCell,latVertex -f "hist.comp.*.nc"

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

Indices for extra dimensions can either be supplied at runtime at a prompt
or through the -d flag.  For each extra dimension, the use can specifiy
nothing (an empty string, meaning skip any fields with this dimension),
a single index, or a comma-separated list of indices or a range of indices
indices (separated by 1 or 2 colons).  For example,

    -d maxEdges= nVertLeves=0:10:2 nParticles=0,2,4,6,8

will ignore any fields with dimension maxEdges, extract every other layer from
the first 10 vertical levels (each into its own field) and extract the five
specified particles.

An index array can also be specified in this way (and these can be mixed with
integer indices in a comma-separated list but not in a colon-separated range):

    -d nVertLeves=0,maxLevelCell

will extract fields from the first vertical level and the vertical level with
index given by maxLevelCell.

The extractor includes optional support for extracting geometry appropriate
for displaying variables at the depth of a topographic feature (typically the
top or bottom of the domain) for MPAS components with a spatially variable
top or bottom index (e.g. `maxLevelCell` in MPAS-Ocean).  This is accomplished
with flags such as:

    --topo_dim=nVertLevels --topo_cell_index=maxLevelCell

Fields on cells are sampled at the topographic index and the geometry includes
polygons corresponding to edges so that vertical faces between adjacent cells
can be displayed.  Fields are extracted as normal except that they are sampled
as point data rather than cell data, allowing computations in ParaView to
display the topography.  A mask field is also included indicating which parts
of edge polygons correspond to the boundary of the domain (boundaryMask == 1)
and which parts of cell and edge polygons are interior (boundaryMask == 0).
Together, this can be used to plot topography by using a calculator filter like
the following:

    coords*(1.0 + 100.0/mag(coords)*((1 - boundaryMask)*(-bottomDepth)
                                     + 10.0*boundaryMask))

If this is entered into a Calculator Filter in ParaView with the "coordinate
result" box checked, the result will to display the MPAS-Ocean topography,
exaggerated by a factor of 100, with a value equivalent to 10 m along boundary
points of edge polygons (a "water-tight" surface).
"""

from __future__ import absolute_import, division, print_function, \
    unicode_literals

from pyevtk.vtk import VtkFile, VtkPolyData

import sys
import os
import glob
import numpy
import argparse
from datetime import datetime
from netCDF4 import date2num

from builtins import input

from netCDF4 import Dataset as NetCDFFile

import xarray
import json
from geometric_features import FeatureCollection
import logging
from io import StringIO

try:
    from progressbar import ProgressBar, Percentage, Bar, ETA
    use_progress_bar = True
except ImportError:
    use_progress_bar = False

from mpas_tools.conversion import mask, cull
from mpas_tools.io import write_netcdf


def extract_vtk(filename_pattern, variable_list='all', dimension_list=None,
                mesh_filename=None, blocking=10000, output_32bit=False,
                combine=False, append=False, out_dir=None, xtime='xtime',
                lonlat=False, time=None, ignore_time=False, topo_dim=None,
                topo_cell_index=None, include_mesh_vars=False,
                fc_region_mask=None, temp_dir='./culled_region'):
    """
    Extract fields from a time series of NetCDF files as VTK files for plotting
    in paraview.

    To extract fields across multiple files, passing in a regular expression
    for the filename pattern, for exmaple ``filename_pattern="hist.comp.*.nc"``.

    By default, time-independent fields on cells are written to a file

    .. code-block::

      vtk_files/staticFieldsOnCells.vtp

    and time-dependent fields on cells are written to

    .. code-block::

      vtk_files/timeDependentFieldsOnCells.pvd
      vtk_files/time_series/timeDependentFieldsOnCells.0.vtp
      vtk_files/time_series/timeDependentFieldsOnCells.1.vtp
      ...

    and similarly for edges and vertices.  Time-independent fields can be
    included in each time step of the time-dependent fields for with
    ``combine=True``.  This allows time-dependent and -independent fields
    to be combined in filters within Paraview at the expense of considerable
    additional storage space.

    Indices for extra dimensions can either be supplied at runtime at a prompt
    (if ``dimension_list=None``) or via a list of strings with the dimensions
    and associated indices.  For each extra dimension, you can specifiy
    nothing for the indices (an empty string, meaning skip any fields with this
    dimension), a single index, a comma-separated list of indices, or a range
    of indices (separated by 1 or 2 colons).  For example,

    .. code-block::

      dimension_list=['maxEdges=', 'nVertLeves=0:10:2', 'nParticles=0,2,4,6,8']

    will ignore any fields with dimension ``maxEdges``, extract every other
    layer from the first 10 vertical levels (each into its own field) and
    extract the 5 specified particles.

    An index array can also be specified in this way (and these can be mixed
    with integer indices in a comma-separated list but not in a colon-separated
    range):

    .. code-block::

      dimension_list=['nVertLeves=0,maxLevelCell']

    will extract fields from the first vertical level and the vertical level
    with index given by ``maxLevelCell``.

    The extractor includes optional support for extracting geometry appropriate
    for displaying variables at the depth of a topographic feature (typically
    the top or bottom of the domain) for MPAS components with a spatially
    variable top or bottom index (e.g. ``maxLevelCell`` in MPAS-Ocean).  This is
    accomplished with arguments such as:

    .. code-block::

      topo_dim='nVertLevels', topo_cell_index='maxLevelCell'

    Fields on cells are sampled at the topographic index and the geometry
    includes polygons corresponding to edges so that vertical faces between
    adjacent cells can be displayed.  Fields are extracted as normal except
    that they are sampled as point data rather than cell data, allowing
    computations in ParaView to display the topography.  A mask field is also
    included indicating which parts of edge polygons correspond to the boundary
    of the domain (``boundaryMask == 1``) and which parts of cell and edge
    polygons are interior (``boundaryMask == 0``).  Together, this can be used
    to plot topography by using a calculator filter like the following:

    .. code-block::

      coords*(1.0 + 100.0/mag(coords)*((1 - boundaryMask)*(-bottomDepth)
                                       + 10.0*boundaryMask))

    If this is entered into a Calculator Filter in ParaView with the "coordinate
    result" box checked, the result will to display the MPAS-Ocean topography,
    exaggerated by a factor of 100, with a value equivalent to 10 m along
    boundary points of edge polygons (a "water-tight" surface).

    Parameters
    ----------
    filename_pattern : str
        MPAS Filename pattern

    variable_list: list of str, optional
        List of variables to extract ('all' for all variables, 'allOnCells'
        for all variables on cells, etc.)"

    dimension_list: list of str, optional
        A list of dimensions and associated indices

    mesh_filename: str, optional
        MPAS Mesh filename. By default, the first file matching
        ``filename_pattern`` will be used

    blocking: int, optional
        Size of blocks when reading MPAS file

    output_32bit: bool, optional
        Whether the vtk files will be written using 32bit floats

    combine: bool, optional
        Whether time-independent fields are written to each file along with
        time-dependent fields

    append: bool, optional
        Whether only vtp files that do not already exist are written out

    out_dir: str, optional
        The output directory

    xtime: str, optional"
        The name of the xtime variable or 'none' to extract Time dim without
        xtime

    lonlat: bool, optional
        Whether the resulting points are in lon-lat space, not Cartesian

    time: str, optional
        Indices for the time dimension

    ignore_time: bool, optional
        Whether to ignore the Time dimension if it exists for files with a Time
        dimension but no xtime variable (e.g. mesh file)

    topo_dim: str, optional
        Dimension and range for topography dimension

    topo_cell_index: str, optional
        Index array indicating the bottom of the domain (default is the
        topo_dim-1 for all cells)

    include_mesh_vars: bool, optional
        Whether to include mesh variables as well as time-series variables
        in the extraction

    fc_region_mask: geometric_features.FeatureCollection, optional
        A feature collection used to define a mask.  The MPAS data is culled to
        lie within the mask before conversion to VTK proceeds

    temp_dir: str, optional
        If fc_region_mask is supplied, a temporary directory where the culled
        mesh and time series files are stored
    """

    if ignore_time:
        xtime = None

    (time_indices, time_file_names) = setup_time_indices(
        filename_pattern, xtime)

    if time is not None:
        time_indices, time_file_names = \
            parse_time_indices(time, time_indices, time_file_names)

    separate_mesh_file = True
    if mesh_filename is None:
        mesh_filename = time_file_names[0]
        separate_mesh_file = False

    if fc_region_mask is not None:
        mesh_filename, time_file_names = _cull_files(
            fc_region_mask, temp_dir, mesh_filename, time_file_names,
            separate_mesh_file, variable_list, include_mesh_vars, xtime)
        separate_mesh_file = True

    # Setting dimension values:
    time_series_file = open_netcdf(time_file_names[0])
    if separate_mesh_file:
        mesh_file = open_netcdf(mesh_filename)
    else:
        mesh_file = None
    extra_dims, topo_cell_indices = \
        parse_extra_dims(dimension_list, time_series_file,
                         mesh_file, topo_dim=topo_dim,
                         topo_cell_index_name=topo_cell_index)
    basic_dims = ['nCells', 'nEdges', 'nVertices', 'Time']
    if topo_dim is not None:
        basic_dims.append(topo_dim)

    (all_dim_vals, cellVars, vertexVars, edgeVars) = \
        setup_dimension_values_and_sort_vars(
            time_series_file, mesh_file,  variable_list, extra_dims,
            include_mesh_vars, basic_dims=basic_dims)
    time_series_file.close()
    if mesh_file is not None:
        mesh_file.close()

    summarize_extraction(mesh_filename, time_indices, cellVars,
                         vertexVars, edgeVars)

    # Handle cell variables
    if len(cellVars) > 0:
        print(" -- Extracting cell fields --")

        mesh_file = open_netcdf(mesh_filename)

        # Build cell geometry
        if topo_dim is None:
            (vertices, connectivity, offsets, valid_mask) = \
                build_cell_geom_lists(mesh_file, output_32bit, lonlat)
            cell_to_point_map = None
            boundary_mask = None
        else:
            (vertices, connectivity, offsets, valid_mask, cell_to_point_map,
             boundary_mask) = build_topo_point_and_polygon_lists(
                     mesh_file, output_32bit, lonlat)

        if not separate_mesh_file:
            mesh_file.close()
            mesh_file = None

        build_field_time_series(time_indices, time_file_names, mesh_file,
                                out_dir, blocking, all_dim_vals,
                                'nCells', cellVars, vertices, connectivity,
                                offsets, valid_mask, output_32bit,
                                combine, append, xtime,
                                topo_dim=topo_dim,
                                topo_cell_indices=topo_cell_indices,
                                cell_to_point_map=cell_to_point_map,
                                boundary_mask=boundary_mask)
        if separate_mesh_file:
            mesh_file.close()

        print("")

    if len(vertexVars) > 0:
        print(" -- Extracting vertex fields --")

        mesh_file = open_netcdf(mesh_filename)

        # Build vertex geometry
        (vertices, connectivity, offsets, valid_mask) = \
            build_vertex_geom_lists(mesh_file, output_32bit, lonlat)

        if not separate_mesh_file:
            mesh_file.close()
            mesh_file = None

        build_field_time_series(time_indices, time_file_names, mesh_file,
                                out_dir, blocking, all_dim_vals,
                                'nVertices', vertexVars, vertices,
                                connectivity, offsets, valid_mask, output_32bit,
                                combine, append, xtime)

        if separate_mesh_file:
            mesh_file.close()

        print("")

    if len(edgeVars) > 0:
        print(" -- Extracting edge fields --")

        mesh_file = open_netcdf(mesh_filename)

        # Build cell list
        (vertices, connectivity, offsets, valid_mask) = \
            build_edge_geom_lists(mesh_file, output_32bit, lonlat)

        if not separate_mesh_file:
            mesh_file.close()
            mesh_file = None

        build_field_time_series(time_indices, time_file_names, mesh_file,
                                out_dir, blocking, all_dim_vals,
                                'nEdges', edgeVars, vertices, connectivity,
                                offsets, valid_mask, output_32bit,
                                combine, append, xtime)

        if separate_mesh_file:
            mesh_file.close()


def main():
    if use_progress_bar:
        print(" -- Using progress bars --")
    else:
        print(" -- Progress bars are not available--")
    parser = \
        argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-f", "--file_pattern", dest="filename_pattern",
                        help="MPAS Filename pattern.", metavar="FILE",
                        required=True)
    parser.add_argument("-m", "--mesh_file", dest="mesh_filename",
                        help="MPAS Mesh filename. If not set, it will use the "
                             "first file in the -f flag as the mesh file.")
    parser.add_argument("-b", "--blocking", dest="blocking", type=int,
                        help="Size of blocks when reading MPAS file",
                        metavar="BLK", default=32000)
    parser.add_argument("-v", "--variable_list", dest="variable_list",
                        help="List of variables to extract ('all' for all "
                             "variables, 'allOnCells' for all variables on "
                             "cells, etc.)", metavar="VAR", required=True)
    parser.add_argument("-3", "--32bit", dest="output_32bit",
                        help="If set, the vtk files will be written using "
                             "32bit floats.", action="store_true")
    parser.add_argument("-c", "--combine", dest="combine_output",
                        help="If set, time-independent fields are written to "
                             "each file along with time-dependent fields.",
                        action="store_true")
    parser.add_argument("-a", "--append", dest="append",
                        help="If set, only vtp files that do not already "
                             "exist are written out.", action="store_true")
    parser.add_argument("-d", "--dim_list", dest="dimension_list", nargs="+",
                        help="A list of dimensions and associated indices.",
                        metavar="DIM")
    parser.add_argument("-o", "--out_dir", dest="out_dir",
                        help="the output directory.", default='vtk_files',
                        metavar="DIR")
    parser.add_argument("-x", "--xtime", dest="xtime", default='xtime',
                        metavar="XTIME",
                        help="the name of the xtime variable or 'none' to "
                             "extract Time dim without xtime")
    parser.add_argument("-l", "--lonlat", dest="lonlat",
                        help="If set, the resulting points are in lon-lat "
                             "space, not Cartesian.", action="store_true")
    parser.add_argument("-t", "--time", dest="time",
                        help="Indices for the time dimension", metavar="TIME",
                        required=False)
    parser.add_argument("--ignore_time", dest="ignore_time", required=False,
                        action="store_true",
                        help="ignore the Time dimension if it exists "
                             "for files with a Time dimension but no xtime"
                             "variable (e.g. mesh file)")
    parser.add_argument("--topo_dim", dest="topo_dim", required=False,
                        help="Dimension and range for topography dimension")
    parser.add_argument("--topo_cell_index", dest="topo_cell_index",
                        required=False,
                        help="Index array indicating the bottom of the domain "
                             "(default is the topo_dim-1 for all cells)")
    parser.add_argument("--include_mesh_vars", dest="include_mesh_vars",
                        action="store_true",
                        help="Whether to extract mesh variables as well as"
                             "time-series variables")
    parser.add_argument("--region_mask", dest="region_mask", required=False,
                        help="A geojson file defining a region that the data"
                             "should be masked to before extraction.  Make one"
                             "easily at https://goejson.io")
    parser.add_argument("--temp_dir", dest="temp_dir", required=False,
                        default="./culled_region",
                        help="If --region_mask is provided, a temporary "
                             "directory for the culled files")
    args = parser.parse_args()

    if args.region_mask is not None:
        fc_region_mask = FeatureCollection()
        with open(args.region_mask) as f:
            featuresDict = json.load(f)
            defaults = {'component': 'ocean',
                        'name': 'mask',
                        'object': 'region'}
            for feature in featuresDict['features']:
                if feature['geometry']['type'] not in ['Polygon',
                                                       'MultiPolygon']:
                    raise ValueError('All masking features must be regions '
                                     '(Polygons or MultiPolygons)')
                # assign the default values if they're not already present
                for key, value in defaults.items():
                    if key not in feature['properties']:
                        feature['properties'][key] = value
                fc_region_mask.add_feature(feature)
    else:
        fc_region_mask = None

    extract_vtk(filename_pattern=args.filename_pattern,
                variable_list=args.variable_list,
                dimension_list=args.dimension_list,
                mesh_filename=args.mesh_filename,
                blocking=args.blocking,
                output_32bit=args.output_32bit,
                combine=args.combine_output,
                append=args.append,
                out_dir=args.out_dir,
                xtime=args.xtime,
                lonlat=args.lonlat,
                time=args.time,
                ignore_time=args.ignore_time,
                topo_dim=args.topo_dim,
                topo_cell_index=args.topo_cell_index,
                include_mesh_vars=args.include_mesh_vars,
                fc_region_mask=fc_region_mask,
                temp_dir=args.temp_dir)


def build_field_time_series(local_time_indices, file_names, mesh_file,
                            out_dir, blocking, all_dim_vals, blockDimName,
                            variable_list, vertices, connectivity, offsets,
                            valid_mask, output_32bit, combine_output, append,
                            xtimeName, topo_dim=None, topo_cell_indices=None,
                            cell_to_point_map=None, boundary_mask=None):  # {{{

    if len(variable_list) == 0:
        return

    if output_32bit:
        outType = 'float32'
    else:
        outType = 'float64'

    # Get dimension info to allocate the size of Colors
    time_series_file = open_netcdf(file_names[0])

    if mesh_file is not None:
        # blockDim may not exist in time series file
        blockDim = len(mesh_file.dimensions[blockDimName])
    else:
        blockDim = len(time_series_file.dimensions[blockDimName])

    if boundary_mask is not None:
        variable_list.append('boundaryMask')
        all_dim_vals['boundaryMask'] = None
        pointData = True
        cellData = False
    else:
        pointData = False
        cellData = True

    # Pre-compute the number of blocks
    nBlocks = int(numpy.ceil(blockDim / blocking))

    nPolygons = len(offsets)
    nPoints = len(vertices[0])
    nTimes = len(local_time_indices)
    nVars = len(variable_list)

    var_has_time_dim = numpy.zeros(nVars, bool)
    nHyperSlabs = 0
    for iVar in range(nVars):
        var_name = variable_list[iVar]
        if boundary_mask is not None and var_name == 'boundaryMask':
            var_has_time_dim[iVar] = False
        elif xtimeName is not None:
            if var_name in time_series_file.variables:
                var_has_time_dim[iVar] = \
                    'Time' in time_series_file.variables[var_name].dimensions
            else:
                # we can't support time dependence in the mesh file
                assert('Time' not in mesh_file.variables[var_name].dimensions)
                var_has_time_dim[iVar] = False

        extra_dim_vals = all_dim_vals[var_name]
        if (extra_dim_vals is None) or (len(extra_dim_vals) == 0):
            nHyperSlabs += 1
        else:
            nHyperSlabs += len(extra_dim_vals)

    any_var_has_time_dim = numpy.any(var_has_time_dim)

    if topo_dim is not None:
        if (mesh_file is not None) and (topo_dim in mesh_file.dimensions):
            nTopoLevels = len(mesh_file.dimensions[topo_dim])
        else:
            nTopoLevels = len(time_series_file.dimensions[topo_dim])
    else:
        nTopoLevels = None

    time_series_file.close()

    try:
        os.makedirs(out_dir)
    except OSError:
        pass

    if any_var_has_time_dim:
        try:
            os.makedirs('{}/time_series'.format(out_dir))
        except OSError:
            pass
    else:
        # there is no point in combining output if no fields have Time dim
        combine_output = False
        nTimes = 1

    # Output time series
    if use_progress_bar:
        widgets = ['Writing time series: ', Percentage(), ' ', Bar(), ' ',
                   ETA()]
        field_bar = ProgressBar(widgets=widgets,
                                maxval=nTimes*nHyperSlabs).start()
    else:
        print("Writing time series....")

    suffix = blockDimName[1:]
    if any_var_has_time_dim:
        if combine_output or numpy.all(var_has_time_dim):
            out_prefix = "fieldsOn{}".format(suffix)
        else:
            out_prefix = "timeDependentFieldsOn{}".format(suffix)
        # start the pvd file
        pvd_file = write_pvd_header(out_dir, out_prefix)
        pvd_file.write('<Collection>\n')

    if not combine_output and not numpy.all(var_has_time_dim):
        static_prefix = "staticFieldsOn{}".format(suffix)
        varIndices = numpy.arange(nVars)[numpy.logical_not(var_has_time_dim)]
        timeIndependentFile = write_vtp_header(out_dir,
                                               static_prefix,
                                               varIndices[0],
                                               varIndices,
                                               variable_list,
                                               all_dim_vals,
                                               vertices,
                                               connectivity,
                                               offsets,
                                               nPoints,
                                               nPolygons,
                                               outType,
                                               cellData=cellData,
                                               pointData=pointData,
                                               xtime=None)

    prev_file = ""
    for time_index in range(nTimes):

        if prev_file != file_names[time_index]:
            if prev_file != "":
                time_series_file.close()
            time_series_file = open_netcdf(file_names[time_index])
            prev_file = file_names[time_index]

        if any_var_has_time_dim:
            if xtimeName is None:
                xtime = None
                years = float(time_index)
            else:
                if xtimeName == 'none':
                    xtime = '{}'.format(time_index)
                    years = float(time_index)
                else:
                    if xtimeName not in time_series_file.variables:
                        raise ValueError("xtime variable name {} not found in "
                                         "{}".format(xtimeName,
                                                     time_series_file))
                    var = time_series_file.variables[xtimeName]
                    if len(var.shape) == 2:
                        xtime = var[local_time_indices[time_index], :]
                        xtime = xtime.tostring().decode('utf-8').strip().strip(
                            '\x00')
                        date = datetime(int(xtime[0:4]), int(xtime[5:7]),
                                        int(xtime[8:10]), int(xtime[11:13]),
                                        int(xtime[14:16]), int(xtime[17:19]))
                        years = date2num(date, units='days since 0000-01-01',
                                         calendar='noleap')/365.
                    else:
                        xtime = var[local_time_indices[time_index]]
                        years = xtime/365.
                        xtime = str(xtime)

            # write the header for the vtp file
            vtp_file_prefix = "time_series/{}.{:d}".format(out_prefix,
                                                           time_index)
            file_name = '{}/{}.vtp'.format(out_dir, vtp_file_prefix)
            if append and os.path.exists(file_name):
                pvd_file.write('<DataSet timestep="{:.16f}" group="" '
                               'part="0"\n'.format(years))
                pvd_file.write('\tfile="{}.vtp"/>\n'.format(vtp_file_prefix))
                continue

            if combine_output:
                varIndices = numpy.arange(nVars)
            else:
                varIndices = numpy.arange(nVars)[var_has_time_dim]
            timeDependentFile = write_vtp_header(out_dir,
                                                 vtp_file_prefix,
                                                 varIndices[0],
                                                 varIndices,
                                                 variable_list,
                                                 all_dim_vals,
                                                 vertices,
                                                 connectivity,
                                                 offsets,
                                                 nPoints,
                                                 nPolygons,
                                                 outType,
                                                 cellData=cellData,
                                                 pointData=pointData,
                                                 xtime=xtime)

            # add time step to pdv file
            pvd_file.write('<DataSet timestep="{:.16f}" group="" '
                           'part="0"\n'.format(years))
            pvd_file.write('\tfile="{}.vtp"/>\n'.format(vtp_file_prefix))

        if time_index == 0 or combine_output:
            varIndices = numpy.arange(nVars)
        else:
            # only the time-dependent variables
            varIndices = numpy.arange(nVars)[var_has_time_dim]

        iHyperSlabProgress = 0
        for iVar in varIndices:
            has_time = var_has_time_dim[iVar]

            var_name = variable_list[iVar]
            (out_var_names, dim_list) = \
                get_hyperslab_name_and_dims(var_name, all_dim_vals[var_name])
            if has_time or combine_output:
                vtkFile = timeDependentFile
            else:
                vtkFile = timeIndependentFile
            for iHyperSlab in range(len(out_var_names)):
                if dim_list is not None:
                    dim_vals = dim_list[iHyperSlab]
                else:
                    dim_vals = None

                if boundary_mask is not None and var_name == 'boundaryMask':
                    field = numpy.array(boundary_mask, dtype=outType)
                else:
                    field = numpy.zeros(blockDim, dtype=outType)

                    for iBlock in numpy.arange(0, nBlocks):
                        blockStart = iBlock * blocking
                        blockEnd = min((iBlock + 1) * blocking, blockDim)
                        block_indices = numpy.arange(blockStart, blockEnd)
                        if topo_cell_indices is None:
                            block_topo_cell_indices = None
                        else:
                            block_topo_cell_indices = \
                                topo_cell_indices[block_indices]
                        field_block = read_field(
                            var_name, mesh_file, time_series_file,
                            dim_vals, local_time_indices[time_index],
                            block_indices, outType,  topo_dim=topo_dim,
                            topo_cell_indices=block_topo_cell_indices,
                            nTopoLevels=nTopoLevels)

                        field[blockStart:blockEnd] = field_block

                    field = field[valid_mask]

                    if cell_to_point_map is not None:
                        # map field from cells to points
                        field = field[cell_to_point_map]

                vtkFile.appendData(field)

                if use_progress_bar:
                    field_bar.update(time_index*nHyperSlabs +
                                     iHyperSlabProgress)
                    iHyperSlabProgress += 1

                del field

        if any_var_has_time_dim:
            timeDependentFile.save()
            del timeDependentFile

        if time_index == 0 and not combine_output and not \
                numpy.all(var_has_time_dim):
            timeIndependentFile.save()
            del timeIndependentFile

    time_series_file.close()
    if use_progress_bar:
        field_bar.finish()

    if any_var_has_time_dim:
        # finish the pdv file
        pvd_file.write('</Collection>\n')
        pvd_file.write('</VTKFile>\n')
        pvd_file.close()  # }}}


def open_netcdf(file_name):
    nc_file = NetCDFFile(file_name, 'r')
    # turn off auto mask (if applicable)
    try:
        nc_file.set_auto_mask(False)
    except AttributeError:
        pass
    return nc_file


def is_valid_mesh_var(mesh_file, variable_name):  # {{{
    if mesh_file is None:
        return False

    if variable_name not in mesh_file.variables:
        return False

    return 'Time' not in mesh_file.variables[variable_name].dimensions  # }}}


def get_var(variable_name, mesh_file, time_series_file):  # {{{
    if is_valid_mesh_var(mesh_file, variable_name):
        return mesh_file.variables[variable_name]
    else:
        return time_series_file.variables[variable_name]  # }}}


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
        print("No files to process.")
        print("Exiting...")
        sys.exit(0)

    if use_progress_bar:
        widgets = ['Build time indices: ', Percentage(), ' ', Bar(), ' ',
                   ETA()]
        time_bar = ProgressBar(widgets=widgets, maxval=len(file_list)).start()
    else:
        print("Build time indices...")

    i_file = 0
    allTIndex = 0
    for file_name in file_list:
        try:
            nc_file = open_netcdf(file_name)
        except IOError:
            print("Warning: could not open {}".format(file_name))
            continue

        if 'Time' not in nc_file.dimensions or xtimeName is None:
            local_times = ['0']
        else:
            local_times = []
            if xtimeName == 'none':
                # no xtime variable so just use integers converted to strings
                for index in range(len(nc_file.dimensions['Time'])):
                    local_times.append(allTIndex)
                    allTIndex += 1
            else:
                if xtimeName not in nc_file.variables:
                    raise ValueError("xtime variable name {} not found in "
                                     "{}".format(xtimeName, file_name))
                xtime = nc_file.variables[xtimeName]
                if len(xtime.shape) == 2:
                    xtime = xtime[:, :]
                    for index in range(xtime.shape[0]):
                        local_times.append(xtime[index, :].tostring())
                else:
                    local_times = xtime[:]

                if len(local_times) == 0:
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

    return local_indices, file_names  # }}}


def parse_extra_dim(dim_name, index_string, time_series_file, mesh_file):
    # {{{
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

    time_series_file : NetCDF4.Dataset
        The name of a time series file that can be used to determine the size
        of the dimension if ``mesh_file=None``.

    mesh_file : NetCDF4.Dataset
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
        pad = int(numpy.log10(max(max_index, 1)))+1
        template = '%%0%dd' % pad
        for i in range(len(indices)):
            try:
                val = int(indices[i])
            except ValueError:
                continue
            indices[i] = template % (val)

    return indices

# }}}


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
    index_string : str or list of int
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
    _, numerical_indices = parse_index_string(index_string, dim_size)

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
    index_string : str or list of int
        An index string indicating with indices are to be extracted

    dim_size : int
        The size of the dimension

    Returns
    -------
    indices : list of int
        The indices corresponding to the given index string.
    """
    if not isinstance(index_string, str):
        numerical_indices = index_string
        indices = []
        for index in numerical_indices:
            if index < 0 or index >= dim_size:
                raise ValueError("Index (or indices) out of bounds 0 <= "
                                 "index < {}: {}".format(dim_size,
                                                         index_string))
            indices.append('{}'.format(index))
    else:
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
                                 "index < {}: {}".format(dim_size,
                                                         index_string))

            numerical_indices.append(val)

    return indices, numerical_indices  # }}}


def parse_extra_dims(dimension_list, time_series_file, mesh_file,
                     topo_dim=None, topo_cell_index_name=None,
                     max_index_count=None):
    # {{{
    """
    Parses a list of dimensions and corresponding indices separated by equals
    signs. Optionally, a max_index_count (typically 1) can be provided,
    indicating that indices beyond max_index_count-1 will be ignored in each
    dimension. Optionally, topo_dim contains the name of a dimension associated
    with the surface or bottom topography (e.g. nVertLevels for MPAS-Ocean)
    If topo_dim is provided, topo_cell_index_name can optionally be either
    a constant value for the vertical index to the topography or the name of a
    field with dimension nCells that contains the vertical index of the
    topography.
    """

    extra_dims = {}
    topo_cell_indices = None

    if dimension_list is not None:
        for dim_item in dimension_list:
            print(dim_item)
            (dimName, index_string) = dim_item.split('=')
            indices = parse_extra_dim(dimName, index_string, time_series_file,
                                      mesh_file)
            if indices is not None:
                if max_index_count is None or len(indices) <= max_index_count:
                    extra_dims[dimName] = indices
                else:
                    extra_dims[dimName] = indices[0:max_index_count]

    if topo_dim is not None:
        if topo_cell_index_name is not None:
            if (mesh_file is not None) and \
                    (topo_cell_index_name in mesh_file.variables):
                topo_cell_indices = \
                    mesh_file.variables[topo_cell_index_name][:]-1
            else:
                topo_cell_indices = \
                    time_series_file.variables[topo_cell_index_name][:]-1
        else:
            index = len(mesh_file.dimensions[topo_dim])-1
            nCells = len(mesh_file.dimensions['nCells'])
            topo_cell_indices = index*numpy.ones(nCells, int)

    return extra_dims, topo_cell_indices
# }}}


def setup_dimension_values_and_sort_vars(
        time_series_file, mesh_file,  variable_list, extra_dims,
        include_mesh_vars, basic_dims=('nCells', 'nEdges', 'nVertices', 'Time'),
        include_dims=('nCells', 'nEdges', 'nVertices')):  # {{{
    """
    Creates a list of variables names to be extracted.  Prompts for indices
    of any extra dimensions that were not specified on the command line.
    extra_dims should be a dictionary of indices along extra dimensions (as
    opposed to "basic" dimensions).  basic_dims is a list of dimension names
    that should be excluded from extra_dims.  include_dims is a list of
    possible dimensions, one of which must be in each vairable to be extracted
    (used in expanding command line placeholders "all", "allOnCells", etc.)
    """

    time_series_variables = time_series_file.variables
    if mesh_file is None or not include_mesh_vars:
        mesh_variables = None
    else:
        mesh_variables = mesh_file.variables
    variable_names = _expand_variable_list(variable_list, time_series_variables,
                                           mesh_variables, include_dims)

    all_dim_vals = {}
    cellVars = []
    vertexVars = []
    edgeVars = []

    promptDimNames = []
    display_prompt = True
    for variable_name in variable_names:
        if is_valid_mesh_var(mesh_file, variable_name):
            nc_file = mesh_file
        else:
            nc_file = time_series_file
        field_dims = nc_file.variables[variable_name].dimensions
        for dim in field_dims:
            if ((dim in basic_dims) or (dim in extra_dims)
                    or (dim in promptDimNames)):
                # this dimension has already been accounted for
                continue
            promptDimNames.append(str(dim))

            if display_prompt:
                print("")
                print("Need to define additional dimension values")
                display_prompt = False

            dim_size = len(nc_file.dimensions[dim])
            valid = False
            while not valid:
                print("Valid range for dimension {} between 0 and {}"
                      "".format(dim, dim_size-1))
                index_string = input("Enter a value for dimension {}: "
                                     "".format(dim))
                indices = parse_extra_dim(str(dim), index_string,
                                          time_series_file, mesh_file)
                valid = indices is not None
                if valid:
                    extra_dims[str(dim)] = indices
                else:
                    print(" -- Invalid value, please re-enter --")

    empty_dims = []
    for dim in extra_dims:
        if len(extra_dims[dim]) == 0:
            empty_dims.append(dim)

    for variable_name in variable_names:

        field_dims = get_var(variable_name, mesh_file,
                             time_series_file).dimensions
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
            if dim not in basic_dims:
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
                    dim_vals.append([index0, index1])
        elif len(indices) == 3:
            dim_vals = []
            for index0 in indices[0]:
                for index1 in indices[1]:
                    for index2 in indices[2]:
                        dim_vals.append([index0, index1, index2])
        else:
            print("variable {} has too many extra dimensions and will be "
                  "skipped.".format(variable_name))
            continue

        if "nCells" in field_dims:
            cellVars.append(variable_name)
        elif "nVertices" in field_dims:
            vertexVars.append(variable_name)
        elif "nEdges" in field_dims:
            edgeVars.append(variable_name)

        all_dim_vals[variable_name] = dim_vals
        del dim_vals

    return all_dim_vals, cellVars, vertexVars, edgeVars
# }}}


def summarize_extraction(mesh_file, time_indices, cellVars, vertexVars,
                         edgeVars, transects_file=None):  # {{{
    """
    print a summary of the time levels, mesh file, transects file (optional)
    and variables to be extracted.
    """

    print("")
    print("Extracting a total of {} time levels.".format(len(time_indices)))
    print("Using file '{}' as the mesh file for this extraction."
          "".format(mesh_file))
    if transects_file is not None:
        print("Using file '{}' as the transects file.".format(transects_file))
    print("")
    print("")
    print("The following variables will be extracted from the input file(s).")
    print("")

    if len(cellVars) > 0:
        print("   Variables with 'nCells' as a dimension:")
        for variable_name in cellVars:
            print("      name: {}".format(variable_name))

    if len(vertexVars) > 0:
        print("   Variables with 'nVertices' as a dimension:")
        for variable_name in vertexVars:
            print("      name: {}".format(variable_name))

    if len(edgeVars) > 0:
        print("   Variables with 'nEdges' as adimension:")
        for variable_name in edgeVars:
            print("      name: {}".format(variable_name))

    print("")
# }}}


def write_pvd_header(path, prefix):  # {{{
    pvd_file = open('{}/{}.pvd'.format(path, prefix), 'w')
    pvd_file.write('<?xml version="1.0"?>\n')
    pvd_file.write('<VTKFile type="Collection" version="0.1"\n')
    pvd_file.write('\tbyte_order="LittleEndian"\n')
    pvd_file.write('\tcompressor="vtkZLibDataCompressor">\n')
    return pvd_file  # }}}


def get_hyperslab_name_and_dims(var_name, extra_dim_vals):  # {{{
    if extra_dim_vals is None:
        return [var_name], None
    if len(extra_dim_vals) == 0:
        return [], None
    out_var_names = []
    for hyper_slab in extra_dim_vals:
        pieces = [var_name]
        pieces.extend(hyper_slab)
        out_var_names.append('_'.join(pieces))
    return out_var_names, extra_dim_vals
# }}}


def write_vtp_header(path, prefix, active_var_index, var_indices,
                     variable_list, all_dim_vals, vertices, connectivity,
                     offsets, nPoints, nPolygons, outType, cellData=True,
                     pointData=False, xtime=None):  # {{{
    vtkFile = VtkFile("{}/{}".format(path, prefix), VtkPolyData)

    if xtime is not None:
        vtkFile.openElement(str("metadata"))
        vtkFile.openElement(str("xtime"))
        vtkFile.xml.addText(str(xtime))
        vtkFile.closeElement(str("xtime"))
        vtkFile.closeElement(str("metadata"))

    vtkFile.openElement(vtkFile.ftype.name)
    vtkFile.openPiece(npoints=nPoints, npolys=nPolygons)

    vtkFile.openElement(str("Points"))
    vtkFile.addData(str("points"), vertices)
    vtkFile.closeElement(str("Points"))

    vtkFile.openElement(str("Polys"))
    vtkFile.addData(str("connectivity"), connectivity)
    vtkFile.addData(str("offsets"), offsets)
    vtkFile.closeElement(str("Polys"))

    if cellData:
        vtkFile.openData(str("Cell"),
                         scalars=[str(var) for var in
                                  variable_list[active_var_index]])
        for iVar in var_indices:
            var_name = variable_list[iVar]
            (out_var_names, dim_list) = \
                get_hyperslab_name_and_dims(var_name, all_dim_vals[var_name])
            for out_var_name in out_var_names:
                vtkFile.addHeader(str(out_var_name), outType, nPolygons, 1)
        vtkFile.closeData(str("Cell"))
    if pointData:
        vtkFile.openData(str("Point"),
                         scalars=[str(var) for var in
                                  variable_list[active_var_index]])
        for iVar in var_indices:
            var_name = variable_list[iVar]
            (out_var_names, dim_list) = \
                get_hyperslab_name_and_dims(var_name, all_dim_vals[var_name])
            for out_var_name in out_var_names:
                vtkFile.addHeader(str(out_var_name), outType, nPoints, 1)
        vtkFile.closeData(str("Point"))

    vtkFile.closePiece()
    vtkFile.closeElement(vtkFile.ftype.name)

    vtkFile.appendData(vertices)
    vtkFile.appendData(connectivity)
    vtkFile.appendData(offsets)

    return vtkFile  # }}}


def build_topo_point_and_polygon_lists(nc_file, output_32bit, lonlat):  # {{{

    if output_32bit:
        dtype = 'f4'
    else:
        dtype = 'f8'

    xVertex, yVertex, zVertex = \
        _build_location_list_xyz(nc_file, 'Vertex', output_32bit, lonlat)

    nCells = len(nc_file.dimensions['nCells'])
    nEdges = len(nc_file.dimensions['nEdges'])
    maxEdges = len(nc_file.dimensions['maxEdges'])

    nEdgesOnCell = nc_file.variables['nEdgesOnCell'][:]
    verticesOnCell = nc_file.variables['verticesOnCell'][:, :]-1
    edgesOnCell = nc_file.variables['edgesOnCell'][:, :]-1
    verticesOnEdge = nc_file.variables['verticesOnEdge'][:] - 1
    cellsOnEdge = nc_file.variables['cellsOnEdge'][:] - 1

    # 4 points for each edge face
    nPoints = 4*nEdges
    # 1 polygon for each edge and cell
    nPolygons = nEdges + nCells

    X = numpy.zeros(nPoints, dtype)
    Y = numpy.zeros(nPoints, dtype)
    Z = numpy.zeros(nPoints, dtype)

    # The points on an edge are vertex 0, 1, 1, 0 on that edge, making a
    # vertical rectangle if the points are offset
    iEdges, voe = numpy.meshgrid(numpy.arange(nEdges), [0, 1, 1, 0],
                                 indexing='ij')
    iVerts = verticesOnEdge[iEdges, voe].ravel()
    X[:] = xVertex[iVerts]
    Y[:] = yVertex[iVerts]
    Z[:] = zVertex[iVerts]
    vertices = (X, Y, Z)

    verticesOnPolygon = -1*numpy.ones((nPolygons, maxEdges), int)
    verticesOnPolygon[0:nEdges, 0:4] = \
        numpy.arange(4*nEdges).reshape(nEdges, 4)

    # Build cells
    if use_progress_bar:
        widgets = ['Build cell connectivity: ', Percentage(), ' ', Bar(), ' ',
                   ETA()]
        bar = ProgressBar(widgets=widgets, maxval=nCells).start()
    else:
        print("Build cell connectivity...")

    outIndex = nEdges

    for iCell in range(nCells):
        neoc = nEdgesOnCell[iCell]
        eocs = edgesOnCell[iCell, 0:neoc]
        vocs = verticesOnCell[iCell, 0:neoc]
        for index in range(neoc):
            iVert = vocs[index]
            iEdge = eocs[index]
            # which vertex on the edge corresponds to iVert?
            coes = cellsOnEdge[iEdge, :]
            voes = verticesOnEdge[iEdge, :]

            if coes[0] == iCell:
                if voes[0] == iVert:
                    voe = 0
                else:
                    voe = 1
            else:
                if voes[0] == iVert:
                    voe = 3
                else:
                    voe = 2

            verticesOnPolygon[nEdges+iCell, index] = 4*iEdge + voe

        outIndex += neoc

        if use_progress_bar:
            bar.update(iCell)

    if use_progress_bar:
        bar.finish()

    validVerts = verticesOnPolygon >= 0

    if lonlat:
        lonEdge = numpy.rad2deg(nc_file.variables['lonEdge'][:])
        latEdge = numpy.rad2deg(nc_file.variables['latEdge'][:])
        lonCell = numpy.rad2deg(nc_file.variables['lonCell'][:])
        latCell = numpy.rad2deg(nc_file.variables['latCell'][:])
        lonPolygon = numpy.append(lonEdge, lonCell)
        latPolygon = numpy.append(latEdge, latCell)

        vertices, verticesOnPolygon = _fix_lon_lat_vertices(vertices,
                                                            verticesOnPolygon,
                                                            validVerts,
                                                            lonPolygon)

    if nc_file.on_a_sphere.strip() == 'NO' and \
            nc_file.is_periodic.strip() == 'YES':
        if lonlat:
            xcoord = lonPolygon
            ycoord = latPolygon
        else:
            xEdge = numpy.rad2deg(nc_file.variables['xEdge'][:])
            yEdge = numpy.rad2deg(nc_file.variables['yEdge'][:])
            xCell = numpy.rad2deg(nc_file.variables['xCell'][:])
            yCell = numpy.rad2deg(nc_file.variables['yCell'][:])
            xcoord = numpy.append(xEdge, xCell)
            ycoord = numpy.append(yEdge, yCell)

        vertices, verticesOnPolygon = _fix_periodic_vertices(vertices,
                                                             verticesOnPolygon,
                                                             validVerts,
                                                             xcoord, ycoord,
                                                             nc_file.x_period,
                                                             nc_file.y_period)

    nPoints = len(vertices[0])

    # we want to know the cells corresponding to each point.  The first two
    # points correspond to the first cell, the second two to the second cell
    # (if any).
    cell_to_point_map = -1*numpy.ones((nPoints,), int)
    boundary_mask = numpy.zeros((nPoints,), bool)

    # first cell on edge always exists
    coe = cellsOnEdge[:, 0].copy()
    for index in range(2):
        voe = verticesOnPolygon[0:nEdges, index]
        cell_to_point_map[voe] = coe
        boundary_mask[voe] = False

    # second cell on edge may not exist
    coe = cellsOnEdge[:, 1].copy()
    mask = coe == -1
    # use the first cell if the second doesn't exist
    coe[mask] = cellsOnEdge[:, 0][mask]
    for index in range(2, 4):
        voe = verticesOnPolygon[0:nEdges, index]
        cell_to_point_map[voe] = coe
        boundary_mask[voe] = mask

    # for good measure, make sure vertices on cell are also accounted for
    for index in range(maxEdges):
        iCells = numpy.arange(nCells)
        voc = verticesOnPolygon[nEdges:nEdges+nCells, index]
        mask = index < nEdgesOnCell
        cell_to_point_map[voc[mask]] = iCells[mask]
        boundary_mask[voc[mask]] = False

    connectivity = verticesOnPolygon[validVerts]
    validCount = numpy.sum(numpy.array(validVerts, int), axis=1)
    offsets = numpy.cumsum(validCount, dtype=int)
    valid_mask = numpy.ones(nCells, bool)

    return vertices, connectivity, offsets, valid_mask, \
        cell_to_point_map, boundary_mask.ravel()  # }}}


def build_cell_geom_lists(nc_file, output_32bit, lonlat):  # {{{

    print("Build geometry for fields on cells...")

    vertices = _build_location_list_xyz(nc_file, 'Vertex', output_32bit,
                                        lonlat)

    if lonlat:
        lonCell = numpy.rad2deg(nc_file.variables['lonCell'][:])
        latCell = numpy.rad2deg(nc_file.variables['latCell'][:])

    nCells = len(nc_file.dimensions['nCells'])

    nEdgesOnCell = nc_file.variables['nEdgesOnCell'][:]
    verticesOnCell = nc_file.variables['verticesOnCell'][:, :] - 1
    # MPAS-O sets non-masked values to total number of vertices instead of 0
    # (as produced in mesh workflow)
    verticesOnCell[numpy.where(verticesOnCell == len(vertices[0]))] = 0

    validVertices = numpy.zeros(verticesOnCell.shape, bool)
    for vIndex in range(validVertices.shape[1]):
        validVertices[:, vIndex] = vIndex < nEdgesOnCell

    if lonlat:
        vertices, verticesOnCell = _fix_lon_lat_vertices(vertices,
                                                         verticesOnCell,
                                                         validVertices,
                                                         lonCell)

    if nc_file.on_a_sphere.strip() == 'NO' and \
            nc_file.is_periodic.strip() == 'YES':
        if lonlat:
            xcoord = lonCell
            ycoord = latCell
        else:
            xcoord = nc_file.variables['xCell'][:]
            ycoord = nc_file.variables['yCell'][:]
        vertices, verticesOnCell = _fix_periodic_vertices(vertices,
                                                          verticesOnCell,
                                                          validVertices,
                                                          xcoord, ycoord,
                                                          nc_file.x_period,
                                                          nc_file.y_period)

    connectivity = verticesOnCell[validVertices]
    offsets = numpy.cumsum(nEdgesOnCell, dtype=int)
    valid_mask = numpy.ones(nCells, bool)

    return vertices, connectivity, offsets, valid_mask  # }}}


def build_vertex_geom_lists(nc_file, output_32bit, lonlat):  # {{{
    print("Build geometry for fields on vertices....")

    vertices = _build_location_list_xyz(nc_file, 'Cell', output_32bit, lonlat)

    if lonlat:
        lonVertex = numpy.rad2deg(nc_file.variables['lonVertex'][:])
        latVertex = numpy.rad2deg(nc_file.variables['latVertex'][:])

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

    if nc_file.on_a_sphere.strip() == 'NO' and \
            nc_file.is_periodic.strip() == 'YES':
        # all remaining entries in cellsOnVertex are valid
        validVertices = numpy.ones(cellsOnVertex.shape, bool)
        if lonlat:
            xcoord = lonVertex[valid_mask]
            ycoord = latVertex[valid_mask]
        else:
            xcoord = nc_file.variables['xVertex'][valid_mask]
            ycoord = nc_file.variables['yVertex'][valid_mask]
        vertices, cellsOnVertex = _fix_periodic_vertices(vertices,
                                                         cellsOnVertex,
                                                         validVertices,
                                                         xcoord, ycoord,
                                                         nc_file.x_period,
                                                         nc_file.y_period)

    connectivity = cellsOnVertex.ravel()
    validCount = cellsOnVertex.shape[0]
    offsets = vertexDegree*numpy.arange(1, validCount+1)

    return vertices, connectivity, offsets, valid_mask  # }}}


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
        latEdge = numpy.rad2deg(nc_file.variables['latEdge'][:])

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
        vertices, vertsOnCell = _fix_lon_lat_vertices(vertices,
                                                      vertsOnCell,
                                                      validVerts,
                                                      lonEdge[valid_mask])
    if nc_file.on_a_sphere.strip() == 'NO' and \
            nc_file.is_periodic.strip() == 'YES':
        if lonlat:
            xcoord = lonEdge[valid_mask]
            ycoord = latEdge[valid_mask]
        else:
            xcoord = nc_file.variables['xEdge'][valid_mask]
            ycoord = nc_file.variables['yEdge'][valid_mask]

        vertices, vertsOnCell = _fix_periodic_vertices(vertices,
                                                       vertsOnCell,
                                                       validVerts,
                                                       xcoord, ycoord,
                                                       nc_file.x_period,
                                                       nc_file.y_period)

    connectivity = vertsOnCell[validVerts]
    validCount = numpy.sum(numpy.array(validVerts, int), axis=1)
    offsets = numpy.cumsum(validCount, dtype=int)

    return vertices, connectivity, offsets, valid_mask  # }}}


def get_field_sign(field_name):
    if field_name[0] == '-':
        field_name = field_name[1:]
        sign = -1
    else:
        sign = 1

    return field_name, sign


def read_field(var_name, mesh_file, time_series_file, extra_dim_vals,
               time_index, block_indices, outType, sign=1,
               topo_dim=None, topo_cell_indices=None, nTopoLevels=None):  # {{{

    def read_field_with_dims(field_var, dim_vals, temp_shape, outType,
                             index_arrays):  # {{{
        temp_field = numpy.zeros(temp_shape, dtype=outType)
        inDims = len(dim_vals)
        if inDims <= 0 or inDims > 5:
            print('reading field {} with {} dimensions not supported.'
                  ''.format(var_name, inDims))
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

        if topo_dim is not None and topo_dim in field_var.dimensions:
            if len(temp_field.shape) != 2:
                raise ValueError('Field with dimensions {} not supported in '
                                 'topogrpahy extraction mode.'.format(
                                         field_var.dimensions))
            # sample the depth-dependent field at the index of the topography
            temp_field = temp_field[numpy.arange(temp_field.shape[0]),
                                    topo_cell_indices]

        outDims = len(temp_field.shape)

        if outDims <= 0 or outDims > 4:
            print('something went wrong reading field {}, resulting in a temp '
                  'array with {} dimensions.'.format(var_name, outDims))
            sys.exit(1)
        block_indices = numpy.arange(temp_field.shape[0])
        if outDims == 1:
            field = temp_field
        elif outDims == 2:
            field = temp_field[block_indices, index_arrays[0]]
        elif outDims == 3:
            field = temp_field[block_indices, index_arrays[0], index_arrays[1]]
        elif outDims == 4:
            field = temp_field[block_indices, index_arrays[0], index_arrays[1],
                               index_arrays[2]]

        return field  # }}}

    field_var = get_var(var_name, mesh_file, time_series_file)
    try:
        missing_val = field_var.missing_value
    except AttributeError:
        missing_val = -9999999790214767953607394487959552.000000

    dim_vals = []
    extra_dim_index = 0
    shape = field_var.shape
    temp_shape = ()

    index_arrays = []

    for i in range(field_var.ndim):
        dim = field_var.dimensions[i]
        if dim == 'Time':
            dim_vals.append(time_index)
        elif dim in ['nCells', 'nEdges', 'nVertices']:
            dim_vals.append(block_indices)
            temp_shape = temp_shape + (len(block_indices),)
        elif topo_dim is not None and dim == topo_dim:
            dim_vals.append(numpy.arange(nTopoLevels))
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

                index_array_var = get_var(extra_dim_val, mesh_file,
                                          time_series_file)

                # read the appropriate indices from the index_array_var
                index_array = numpy.maximum(0, numpy.minimum(
                        shape[i]-1, index_array_var[block_indices]-1))

                index_arrays.append(index_array)

            extra_dim_index += 1

    field = read_field_with_dims(field_var, dim_vals, temp_shape, outType,
                                 index_arrays)

    field[field == missing_val] = numpy.nan

    return sign*field  # }}}


def compute_zInterface(minLevelCell, maxLevelCell, layerThicknessCell,
                       zMinCell, zMaxCell, dtype, cellsOnEdge=None):
    # {{{

    (nCells, nLevels) = layerThicknessCell.shape

    cellMask = numpy.ones((nCells, nLevels), bool)
    for iLevel in range(nLevels):
        if minLevelCell is not None:
            cellMask[:, iLevel] = numpy.logical_and(cellMask[:, iLevel],
                                                    iLevel >= minLevelCell)
        if maxLevelCell is not None:
            cellMask[:, iLevel] = numpy.logical_and(cellMask[:, iLevel],
                                                    iLevel <= maxLevelCell)

    zInterfaceCell = numpy.zeros((nCells, nLevels+1), dtype=dtype)
    for iLevel in range(nLevels):
        zInterfaceCell[:, iLevel+1] = \
            zInterfaceCell[:, iLevel] \
            + cellMask[:, iLevel]*layerThicknessCell[:, iLevel]

    if zMinCell is not None:
        minLevel = minLevelCell.copy()
        minLevel[minLevel < 0] = nLevels-1
        zOffsetCell = zMinCell - zInterfaceCell[numpy.arange(0, nCells),
                                                minLevel]
    else:
        zOffsetCell = zMaxCell - zInterfaceCell[numpy.arange(0, nCells),
                                                maxLevelCell+1]

    for iLevel in range(nLevels+1):
        zInterfaceCell[:, iLevel] += zOffsetCell

    if cellsOnEdge is None:
        return zInterfaceCell
    else:
        nEdges = cellsOnEdge.shape[0]
        zInterfaceEdge = numpy.zeros((nEdges, nLevels+1), dtype=dtype)

        # Get a list of valid cells on edges and a mask of which are valid
        cellsOnEdgeMask = numpy.logical_and(cellsOnEdge >= 0,
                                            cellsOnEdge < nCells)
        cellIndicesOnEdge = list()
        cellIndicesOnEdge.append(cellsOnEdge[cellsOnEdgeMask[:, 0], 0])
        cellIndicesOnEdge.append(cellsOnEdge[cellsOnEdgeMask[:, 1], 1])

        for iLevel in range(nLevels):
            edgeMask = numpy.zeros(nEdges, bool)
            layerThicknessEdge = numpy.zeros(nEdges, float)
            denom = numpy.zeros(nEdges, float)
            for index in range(2):
                mask = cellsOnEdgeMask[:, index]
                cellIndices = cellIndicesOnEdge[index]
                cellMaskLocal = cellMask[cellIndices, iLevel]

                edgeMask[mask] = numpy.logical_or(edgeMask[mask],
                                                  cellMaskLocal)

                layerThicknessEdge[mask] += \
                    cellMaskLocal*layerThicknessCell[cellIndices, iLevel]
                denom[mask] += 1.0*cellMaskLocal

            layerThicknessEdge[edgeMask] /= denom[edgeMask]

            zInterfaceEdge[:, iLevel+1] = (zInterfaceEdge[:, iLevel]
                                           + edgeMask*layerThicknessEdge)

        if zMinCell is not None:
            refLevelEdge = numpy.zeros(nEdges, int)
            for index in range(2):
                mask = cellsOnEdgeMask[:, index]
                cellIndices = cellIndicesOnEdge[index]
                refLevelEdge[mask] = numpy.maximum(refLevelEdge[mask],
                                                   minLevel[cellIndices])
        else:
            refLevelEdge = (nLevels-1)*numpy.ones(nEdges, int)
            for index in range(2):
                mask = cellsOnEdgeMask[:, index]
                cellIndices = cellIndicesOnEdge[index]
                refLevelEdge[mask] = numpy.minimum(refLevelEdge[mask],
                                                   maxLevelCell[cellIndices]+1)

        zOffsetEdge = numpy.zeros(nEdges, float)
        # add the average of zInterfaceCell at each adjacent cell
        denom = numpy.zeros(nEdges, float)
        for index in range(2):
            mask = cellsOnEdgeMask[:, index]
            cellIndices = cellIndicesOnEdge[index]
            zOffsetEdge[mask] += zInterfaceCell[cellIndices,
                                                refLevelEdge[mask]]
            denom[mask] += 1.0

        mask = denom > 0.
        zOffsetEdge[mask] /= denom[mask]

        # subtract the depth of zInterfaceEdge at the level of the bottom
        zOffsetEdge -= zInterfaceEdge[numpy.arange(nEdges), refLevelEdge]

        for iLevel in range(nLevels+1):
            zInterfaceEdge[:, iLevel] += zOffsetEdge

        return zInterfaceCell, zInterfaceEdge  # }}}


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
    return X, Y, Z

# }}}


def _fix_lon_lat_vertices(vertices, verticesOnCell, validVertices,
                          lonCell):  # {{{

    nCells = verticesOnCell.shape[0]
    nVertices = len(vertices[0])

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


def _fix_periodic_vertices(vertices, verticesOnCell, validVertices,
                           xCell, yCell, xperiod, yperiod):  # {{{

    vertices, verticesOnCell = _fix_periodic_vertices_1D(
            vertices, verticesOnCell, validVertices, xCell, xperiod, dim=0)
    vertices, verticesOnCell = _fix_periodic_vertices_1D(
            vertices, verticesOnCell, validVertices, yCell, yperiod, dim=1)

    return vertices, verticesOnCell  # }}}


def _fix_periodic_vertices_1D(vertices, verticesOnCell, validVertices,
                              coordCell, coordPeriod, dim):  # {{{

    nCells = verticesOnCell.shape[0]
    nVertices = len(vertices[0])

    coordVertex = vertices[dim]

    coordDiff = coordVertex[verticesOnCell] - coordCell.reshape(nCells, 1)

    # which cells have vertices that are out of range?
    coordOutOfRange = numpy.logical_and(
            validVertices,
            numpy.logical_or(coordDiff > coordPeriod / 2.0,
                             coordDiff < -coordPeriod / 2.0))

    coordCellsOutOfRange = numpy.any(coordOutOfRange, axis=1)

    coordValid = validVertices[coordCellsOutOfRange, :]

    coordVerticesToChange = numpy.zeros(verticesOnCell.shape, bool)
    coordVerticesToChange[coordCellsOutOfRange, :] = coordValid

    coordDiff = coordDiff[coordCellsOutOfRange, :][coordValid]
    coordVOC = verticesOnCell[coordCellsOutOfRange, :][coordValid]

    coordNVerticesToAdd = numpy.count_nonzero(coordValid)

    coordVerticesToAdd = numpy.arange(coordNVerticesToAdd) + nVertices
    coordV = coordVertex[coordVOC]
    verticesOnCell[coordVerticesToChange] = coordVerticesToAdd

    # need to shift points outside periodic domain (assumes that mesh is only
    # within one period) can use mod if this is not the case in general
    coordMask = coordDiff > coordPeriod / 2.0
    coordV[coordMask] -= coordPeriod
    coordMask = coordDiff < -coordPeriod / 2.0
    coordV[coordMask] += coordPeriod

    outVertices = []
    for outDim in range(3):
        if outDim == dim:
            outVertices.append(numpy.append(vertices[outDim], coordV))
        else:
            outVertices.append(numpy.append(vertices[outDim],
                                            vertices[outDim][coordVOC]))

    return tuple(outVertices), verticesOnCell  # }}}


def _expand_variable_list(variable_list, time_series_variables,
                          mesh_variables, include_dims):


    if variable_list == 'all':
        variable_names = []
        exclude_dims = ['Time']
        for variable_name in time_series_variables:
            _add_var(time_series_variables, str(variable_name),
                     include_dims, variable_names, exc_dims=None)
        if mesh_variables is not None:
            for variable_name in mesh_variables:
                _add_var(mesh_variables, str(variable_name), include_dims,
                         variable_names, exclude_dims)
    elif isinstance(variable_list, str):
        variable_names = variable_list.split(',')
    else:
        variable_names = variable_list

    for suffix in ['Cells', 'Edges', 'Vertices']:
        include_dim = 'n{}'.format(suffix)
        all_on = 'allOn{}'.format(suffix)
        if (all_on in variable_names) and (include_dim in include_dims):
            variable_names.remove(all_on)
            exclude_dims = ['Time']
            for variable_name in time_series_variables:
                _add_var(time_series_variables, str(variable_name),
                         inc_dims=[include_dim], variable_names=variable_names,
                         exc_dims=None)
            if mesh_variables is not None:
                for variable_name in mesh_variables:
                    _add_var(mesh_variables, str(variable_name),
                             inc_dims=[include_dim],
                             variable_names=variable_names,
                             exc_dims=exclude_dims)

    variable_names.sort()
    return variable_names


def _add_var(variables, var_name, inc_dims, variable_names, exc_dims=None):
    if var_name in variable_names:
        return

    dims = variables[var_name].dimensions
    supported = False
    for d in inc_dims:
        if d in dims:
            supported = True
    if exc_dims is not None:
        for d in exc_dims:
            if d in dims:
                supported = False
    if supported:
        variable_names.append(var_name)


def _cull_files(fc_region_mask, temp_dir, mesh_filename, time_file_names,
    separate_mesh_file, variable_list, include_mesh_vars, xtime):

    mesh_vars = [
        'areaCell', 'cellsOnCell', 'edgesOnCell', 'indexToCellID',
        'latCell', 'lonCell', 'nEdgesOnCell', 'verticesOnCell',
        'xCell', 'yCell', 'zCell', 'angleEdge', 'cellsOnEdge', 'dcEdge',
        'dvEdge', 'edgesOnEdge', 'indexToEdgeID', 'latEdge',
        'lonEdge', 'nEdgesOnCell', 'nEdgesOnEdge', 'verticesOnEdge',
        'xEdge', 'yEdge', 'zEdge', 'areaTriangle',
        'cellsOnVertex', 'edgesOnVertex', 'indexToVertexID',
        'kiteAreasOnVertex', 'latVertex', 'lonVertex', 'xVertex', 'yVertex',
        'zVertex', 'weightsOnEdge']

    try:
        os.makedirs(temp_dir)
    except OSError:
        pass

    log_stream = StringIO()
    logger = logging.getLogger('_cull_files')
    for handler in logger.handlers:
        logger.removeHandler(handler)
    handler = logging.StreamHandler(log_stream)
    logger.addHandler(handler)
    handler.setLevel(logging.INFO)

    # Figure out the variable names we want to extract
    with open_netcdf(time_file_names[0]) as time_series_file:
        time_series_variables = time_series_file.variables
        if separate_mesh_file and include_mesh_vars:
            mesh_file = open_netcdf(mesh_filename)
            mesh_variables = mesh_file.variables
        else:
            mesh_file = None
            mesh_variables = None

        include_dims = ('nCells', 'nEdges', 'nVertices')
        variable_names = _expand_variable_list(variable_list,
                                               time_series_variables,
                                               mesh_variables, include_dims)

        if mesh_file is not None:
            mesh_file.close()

    print('Including variables: {}'.format(', '.join(variable_names)))

    with xarray.open_dataset(mesh_filename) as ds_mesh:
        ds_mesh = ds_mesh[mesh_vars]
        print('Making a region mask file')
        ds_mask = mask(dsMesh=ds_mesh, fcMask=fc_region_mask, logger=logger)
        write_netcdf(ds_mask, '{}/mask.nc'.format(temp_dir))
        print('Cropping mesh to region')
        out_mesh_filename = '{}/mesh.nc'.format(temp_dir)
        ds_culled = cull(dsIn=ds_mesh, dsInverse=ds_mask, logger=logger)
        write_netcdf(ds_culled, out_mesh_filename)

        region_masks = dict()
        cell_mask = ds_mask.regionCellMasks.sum(dim='nRegions') > 0
        region_masks['nCells'] = cell_mask
        region_masks['nVertices'] = \
                ds_mask.regionVertexMasks.sum(dim='nRegions') > 0
        coe = ds_mesh.cellsOnEdge - 1
        valid_cell_on_edge = numpy.logical_and(coe >= 0, cell_mask[coe])
        region_masks['nEdges'] = numpy.logical_or(
            valid_cell_on_edge.isel(TWO=0),
            valid_cell_on_edge.isel(TWO=1))

    if use_progress_bar:
        widgets = ['Cropping time series to region: ', Percentage(), ' ',
                   Bar(), ' ', ETA()]
        bar = ProgressBar(widgets=widgets, maxval=len(time_file_names)).start()
    else:
        print('Cropping time series to region')
        bar = None

    out_time_file_names = []
    for index, filename in enumerate(time_file_names):
        out_filename = '{}/time_series{:04d}.nc'.format(temp_dir, index)
        out_time_file_names.append(out_filename)
        ds_in = xarray.open_dataset(filename)
        ds_out = xarray.Dataset()
        if xtime is None:
            ds_in = ds_in[variable_names]
        else:
            ds_in = ds_in[variable_names + [xtime]]
            ds_out[xtime] = ds_in[xtime]
        for var in ds_in.data_vars:
            for dim in region_masks:
                if dim in ds_in[var].dims:
                    ds_out[var] = ds_in[var].where(region_masks[dim], drop=True)
        write_netcdf(ds_out, out_filename)
        if use_progress_bar:
            bar.update(index+1)
    bar.finish()

    logger.removeHandler(handler)
    handler.close()

    return out_mesh_filename, out_time_file_names

# vim: set expandtab:
