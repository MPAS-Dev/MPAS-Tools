#!/usr/bin/env python
"""
Name: paraview_vtk_field_extractor.py
Authors: Doug Jacobsen, Xylar Asay-Davis

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

Indices for extra dimensions can either be supplied at runtime at a prompt
or through the -d flag.  For each extra dimension, the use can specifiy
nothing (an empty string, meaning skip any fields with this dimension),
a single index, or a comma-separated list of indices or a range of indices
indices (separated by 1 or 2 colons).  For example,

`-d maxEdges= nVertLeves=0:10:2 nParticles=0,2,4,6,8`

will ignore any fields with dimension maxEdges, extract every other layer from
the first 10 vertical levels (each into its own field) and extract the five
specified particles.

An index array can also be specified in this way (and these can be mixed with
integer indices in a comma-separated list but not in a colon-separated range):

`-d nVertLeves=0,maxLevelCell`

will extract fields from the first vertical level and the vertical level with
index given by maxLevelCell.

Requirements:
This script requires access to the following non standard modules:
pyevtk (available from opengeostat channel)
netCDF4
numpy

Optional modules:
progressbar
"""
import os
import numpy as np

from netCDF4 import Dataset as NetCDFFile
import argparse

try:
    from progressbar import ProgressBar, Percentage, Bar, ETA
    use_progress_bar = True
except ImportError:
    use_progress_bar = False

import utils


def build_field_time_series(local_time_indices, file_names, mesh_file,
                            out_dir, blocking, all_dim_vals, blockDimName,
                            variable_list, vertices, connectivity, offsets,
                            valid_mask, output_32bit, combine_output, append,
                            xtimeName):  # {{{

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

    var_has_time_dim = np.zeros(nVars, bool)
    nHyperSlabs = 0
    for iVar in range(nVars):
        var_name = variable_list[iVar]
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

    time_series_file.close()

    any_var_has_time_dim = np.any(var_has_time_dim)

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
        print "Writing time series...."

    suffix = blockDimName[1:]
    if any_var_has_time_dim:
        if combine_output or np.all(var_has_time_dim):
            out_prefix = "fieldsOn{}".format(suffix)
        else:
            out_prefix = "timeDependentFieldsOn{}".format(suffix)
        # start the pvd file
        pvd_file = utils.write_pvd_header(out_dir, out_prefix)
        pvd_file.write('<Collection>\n')

    if not combine_output and not np.all(var_has_time_dim):
        out_prefix = "staticFieldsOn{}".format(suffix)
        varIndices = np.arange(nVars)[np.logical_not(var_has_time_dim)]
        timeIndependentFile = utils.write_vtp_header(out_dir,
                                                     out_prefix,
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
                                                     cellData=True,
                                                     pointData=False,
                                                     xtime=None)

    prev_file = ""
    for time_index in range(nTimes):

        if prev_file != file_names[time_index]:
            if prev_file != "":
                time_series_file.close()
            time_series_file = NetCDFFile(file_names[time_index], 'r')
            prev_file = file_names[time_index]

        if any_var_has_time_dim:
            if xtimeName not in time_series_file.variables:
                raise ValueError("xtime variable name {} not found in "
                                 "{}".format(xtimeName, time_series_file))
            var = time_series_file.variables[xtimeName]
            xtime = ''.join(var[local_time_indices[time_index], :])

            # write the header for the vtp file
            vtp_file_prefix = "time_series/{}.{:d}".format(out_prefix,
                                                           time_index)
            file_name = '{}/{}.vtp'.format(out_dir, vtp_file_prefix)
            if append and os.path.exists(file_name):
                pvd_file.write('<DataSet timestep="{:d}" group="" '
                               'part="0"\n'.format(time_index))
                pvd_file.write('\tfile="{}.vtp"/>\n'.format(vtp_file_prefix))
                continue

            if combine_output:
                varIndices = np.arange(nVars)
            else:
                varIndices = np.arange(nVars)[var_has_time_dim]
            timeDependentFile = utils.write_vtp_header(out_dir,
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
                                                       cellData=True,
                                                       pointData=False,
                                                       xtime=xtime)

            # add time step to pdv file
            pvd_file.write('<DataSet timestep="{:d}" group="" '
                           'part="0"\n'.format(time_index))
            pvd_file.write('\tfile="{}.vtp"/>\n'.format(vtp_file_prefix))

        if time_index == 0 or combine_output:
            varIndices = np.arange(nVars)
        else:
            # only the time-dependent variables
            varIndices = np.arange(nVars)[var_has_time_dim]

        iHyperSlabProgress = 0
        for iVar in varIndices:
            has_time = var_has_time_dim[iVar]

            var_name = variable_list[iVar]
            (out_var_names, dim_list) = \
                utils.get_hyperslab_name_and_dims(var_name,
                                                  all_dim_vals[var_name])
            if has_time or combine_output:
                vtkFile = timeDependentFile
            else:
                vtkFile = timeIndependentFile
            for iHyperSlab in range(len(out_var_names)):
                if dim_list is not None:
                    dim_vals = dim_list[iHyperSlab]
                else:
                    dim_vals = None

                field = np.zeros(blockDim, dtype=outType)

                for iBlock in np.arange(0, nBlocks):
                    blockStart = iBlock * blocking
                    blockEnd = min((iBlock + 1) * blocking, blockDim)
                    cellIndices = np.arange(blockStart, blockEnd)
                    field_block = \
                        utils.read_field(var_name, mesh_file,
                                         time_series_file, dim_vals,
                                         local_time_indices[time_index],
                                         cellIndices, outType)

                    field[blockStart:blockEnd] = field_block

                field = field[valid_mask]

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
                np.all(var_has_time_dim):
            timeIndependentFile.save()
            del timeIndependentFile

    time_series_file.close()
    if use_progress_bar:
        field_bar.finish()

    if any_var_has_time_dim:
        # finish the pdv file
        pvd_file.write('</Collection>\n')
        pvd_file.write('</VTKFile>\n')  # }}}


if __name__ == "__main__":
    if use_progress_bar:
        print " -- Using progress bars --"
    else:
        print " -- Progress bars are not available--"
    parser = \
        argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-f", "--file_pattern", dest="filename_pattern",
                        help="MPAS Filename pattern.", metavar="FILE",
                        required=True)
    parser.add_argument("-m", "--mesh_file", dest="mesh_filename",
                        help="MPAS Mesh filename. If not set, it will use the "
                             "first file in the -f flag as the mesh file.")
    parser.add_argument("-b", "--blocking", dest="blocking",
                        help="Size of blocks when reading MPAS file",
                        metavar="BLK")
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
                        help="A list of dimensions and associated indices.")
    parser.add_argument("-o", "--out_dir", dest="out_dir",
                        help="the output directory.", default='vtk_files',
                        metavar="DIR")
    parser.add_argument("-x", "--xtime", dest="xtime",
                        help="the name of the xtime variable", default='xtime',
                        metavar="XTIME")
    parser.add_argument("-l", "--lonlat", dest="lonlat",
                        help="If set, the resulting points are in lon-lat "
                             "space, not Cartesian.", action="store_true")
    parser.add_argument("-t", "--time", dest="time",
                        help="Indices for the time dimension", metavar="TIME",
                             required=False)
    args = parser.parse_args()

    if not args.output_32bit:
        use_32bit = False
    else:
        use_32bit = True

    if not args.blocking:
        args.blocking = int(10000)
    else:
        args.blocking = int(args.blocking)

    (time_indices, time_file_names) = utils.setup_time_indices(
        args.filename_pattern, args.xtime)

    if args.time:
        time_indices, time_file_names = \
            utils.parse_time_indices(args.time, time_indices, time_file_names)

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
    extra_dims = utils.parse_extra_dims(args.dimension_list, time_series_file,
                                        mesh_file)
    (all_dim_vals, cellVars, vertexVars, edgeVars) = \
        utils.setup_dimension_values_and_sort_vars(time_series_file, mesh_file,
                                                   args.variable_list,
                                                   extra_dims)
    time_series_file.close()
    if(mesh_file is not None):
        mesh_file.close()

    utils.summarize_extraction(args.mesh_filename, time_indices, cellVars,
                               vertexVars, edgeVars)

    # Handle cell variables
    if len(cellVars) > 0:
        print " -- Extracting cell fields --"

        mesh_file = NetCDFFile(args.mesh_filename, 'r')

        # Build cell geometry
        (vertices, connectivity, offsets, valid_mask) = \
            utils.build_cell_geom_lists(mesh_file, use_32bit, args.lonlat)

        if not separate_mesh_file:
            mesh_file.close()
            mesh_file = None

        build_field_time_series(time_indices, time_file_names, mesh_file,
                                args.out_dir, args.blocking, all_dim_vals,
                                'nCells', cellVars, vertices, connectivity,
                                offsets, valid_mask, use_32bit,
                                args.combine_output, args.append, args.xtime)
        if separate_mesh_file:
            mesh_file.close()

        print ""

    if len(vertexVars) > 0:
        print " -- Extracting vertex fields --"

        mesh_file = NetCDFFile(args.mesh_filename, 'r')

        # Build vertex geometry
        (vertices, connectivity, offsets, valid_mask) = \
            utils.build_vertex_geom_lists(mesh_file, use_32bit, args.lonlat)

        if not separate_mesh_file:
            mesh_file.close()
            mesh_file = None

        build_field_time_series(time_indices, time_file_names, mesh_file,
                                args.out_dir, args.blocking, all_dim_vals,
                                'nVertices', vertexVars, vertices,
                                connectivity, offsets, valid_mask, use_32bit,
                                args.combine_output, args.append, args.xtime)

        if separate_mesh_file:
            mesh_file.close()

        print ""

    if len(edgeVars) > 0:
        print " -- Extracting edge fields --"

        mesh_file = NetCDFFile(args.mesh_filename, 'r')

        # Build cell list
        (vertices, connectivity, offsets, valid_mask) = \
            utils.build_edge_geom_lists(mesh_file, use_32bit, args.lonlat)

        if not separate_mesh_file:
            mesh_file.close()
            mesh_file = None

        build_field_time_series(time_indices, time_file_names, mesh_file,
                                args.out_dir, args.blocking, all_dim_vals,
                                'nEdges', edgeVars, vertices, connectivity,
                                offsets, valid_mask, use_32bit,
                                args.combine_output, args.append, args.xtime)

        if separate_mesh_file:
            mesh_file.close()


# vim: set expandtab:
