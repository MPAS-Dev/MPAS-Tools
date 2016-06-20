#!/usr/bin/env python
"""
Name: paraview_vtk_field_extractor.py
Authors: Doug Jacobsen, Xylar Asay-Davis
Date: 06/20/2016

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
import os
import numpy as np

from netCDF4 import Dataset as NetCDFFile
import argparse

try:
    from progressbar import ProgressBar, Percentage, Bar, ETA
    use_progress_bar = True
except:
    use_progress_bar = False

import utils


def build_field_time_series( local_time_indices, file_names, mesh_file, out_dir, blocking, all_dim_vals,
                             blockDimName, variable_list, vertices, connectivity, offsets,
                             valid_mask, output_32bit, combine_output, append ):#{{{

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
        if var_name in time_series_file.variables:
            var_has_time_dim[iVar] = 'Time' in time_series_file.variables[var_name].dimensions
        else:
            # we can't support time dependence in the mesh file
            assert('Time' not in mesh_file.variables[var_name].dimensions)
            var_has_time_dim[iVar] = False

        extra_dim_vals = all_dim_vals[var_name]
        if (extra_dim_vals is None) or (extra_dim_vals.size == 0):
            nHyperSlabs += 1
        else:
            nHyperSlabs += extra_dim_vals.shape[1]

    time_series_file.close()

    any_var_has_time_dim = np.any(var_has_time_dim)

    try:
        os.makedirs(out_dir)
    except OSError:
        pass

    if any_var_has_time_dim:
        try:
            os.makedirs('%s/time_series'%out_dir)
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
        pvd_file = utils.write_pvd_header(out_dir, out_prefix)
        pvd_file.write('<Collection>\n')

    if not combine_output and not np.all(var_has_time_dim):
        out_prefix = "staticFieldsOn%s"%suffix
        varIndices = np.arange(nVars)[var_has_time_dim == False]
        timeIndependentFile = utils.write_vtp_header(out_dir, out_prefix, varIndices[0], varIndices,
                                                     variable_list, all_dim_vals,
                                                     vertices, connectivity, offsets,
                                                     nPoints, nPolygons, outType,
                                                     cellData=True, pointData=False)


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
            file_name = '%s/%s.vtp'%(out_dir, vtp_file_prefix)
            if append and os.path.exists(file_name):
                continue

            if combine_output:
                varIndices = np.arange(nVars)
            else:
                varIndices = np.arange(nVars)[var_has_time_dim]
            timeDependentFile = utils.write_vtp_header(out_dir, vtp_file_prefix, varIndices[0], varIndices,
                                                       variable_list, all_dim_vals,
                                                       vertices, connectivity, offsets,
                                                       nPoints, nPolygons, outType,
                                                       cellData=True, pointData=False)

            # add time step to pdv file
            pvd_file.write('<DataSet timestep="%d" group="" part="0"\n'%(time_index))
            pvd_file.write('\tfile="%s.vtp"/>\n'%(vtp_file_prefix))

        if time_index == 0 or combine_output:
            varIndices = np.arange(nVars)
        else:
            # only the time-dependent variables
            varIndices = np.arange(nVars)[var_has_time_dim]


        iHyperSlabProgress = 0
        for iVar in varIndices:
            has_time = var_has_time_dim[iVar]

            var_name = variable_list[iVar]
            (out_var_names, dim_list) = utils.get_hyperslab_name_and_dims(var_name, all_dim_vals)
            if has_time or combine_output:
                vtkFile = timeDependentFile
            else:
                vtkFile = timeIndependentFile
            for iHyperSlab in range(len(out_var_names)):
                if dim_list is not None:
                    dim_vals = dim_list[:,iHyperSlab]
                else:
                    dim_vals = None

                field_var = utils.get_var(var_name, mesh_file, time_series_file)

                field = np.zeros(blockDim,dtype=outType)

                try:
                    missing_val = field_var.missing_value
                except:
                    missing_val = -9999999790214767953607394487959552.000000

                for iBlock in np.arange(0, nBlocks):
                    blockStart = iBlock * blocking
                    blockEnd = min( (iBlock + 1) * blocking, blockDim )
                    cellIndices = np.arange(blockStart,blockEnd)
                    field_block = utils.read_field(field_var, dim_vals,
                                     local_time_indices[time_index], cellIndices,
                                     outType)

                    field_block[field_block == missing_val] = np.nan
                    field[blockStart:blockEnd] = field_block


                field = field[valid_mask]

                vtkFile.appendData(field)

                if use_progress_bar:
                    field_bar.update(time_index*nHyperSlabs + iHyperSlabProgress)
                    iHyperSlabProgress += 1

                del field
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
    parser.add_argument("-v", "--variable_list", dest="variable_list", help="List of variables to extract ('all' for all variables, 'allOnCells' for all variables on cells, etc.)", metavar="VAR", required=True)
    parser.add_argument("-3", "--32bit", dest="output_32bit", help="If set, the vtk files will be written using 32bit floats.", action="store_true")
    parser.add_argument("-c", "--combine", dest="combine_output", help="If set, time-independent fields are written to each file along with time-dependent fields.", action="store_true")
    parser.add_argument("-a", "--append", dest="append", help="If set, only vtp files that do not already exist are written out.", action="store_true")
    parser.add_argument("-d", "--dim_list", dest="dimension_list", nargs="+", help="A list of dimensions and associated indices.")
    parser.add_argument("-o", "--out_dir", dest="out_dir", help="the output directory.", default='vtk_files', metavar="DIR")
    args = parser.parse_args()

    if not args.output_32bit:
        use_32bit = False
    else:
        use_32bit = True

    if not args.blocking:
        args.blocking = int(10000)
    else:
        args.blocking = int(args.blocking)

    (time_indices, time_file_names) = utils.setup_time_indices(args.filename_pattern)

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
    extra_dims = utils.parse_extra_dims(args.dimension_list, time_series_file, mesh_file)
    (all_dim_vals, cellVars, vertexVars, edgeVars) = utils.setup_dimension_values_and_sort_vars(
            time_series_file, mesh_file, args.variable_list, extra_dims)
    time_series_file.close()
    if(mesh_file is not None):
        mesh_file.close()

    utils.summarize_extraction(args.mesh_filename, time_indices, cellVars, vertexVars, edgeVars)

    # Handle cell variables
    if len(cellVars) > 0:
        print " -- Extracting cell fields --"

        mesh_file = NetCDFFile(args.mesh_filename, 'r')

        # Build vertex list
        vertices = utils.build_location_list_xyz( mesh_file, 'xVertex', 'yVertex', 'zVertex', use_32bit )

        # Build cell list
        (connectivity, offsets, valid_mask) = utils.build_cell_lists( mesh_file, args.blocking )

        if not separate_mesh_file:
            mesh_file.close()
            mesh_file = None

        build_field_time_series( time_indices, time_file_names, mesh_file, args.out_dir, args.blocking,
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
        vertices = utils.build_location_list_xyz( mesh_file, 'xCell', 'yCell', 'zCell', use_32bit )

        # Build cell list
        (connectivity, offsets, valid_mask) = utils.build_dual_cell_lists( mesh_file, args.blocking )

        if not separate_mesh_file:
            mesh_file.close()
            mesh_file = None

        build_field_time_series( time_indices, time_file_names, mesh_file, args.out_dir, args.blocking,
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
        verticesC = utils.build_location_list_xyz( mesh_file, 'xCell', 'yCell', 'zCell', use_32bit )
        verticesV = utils.build_location_list_xyz( mesh_file, 'xVertex', 'yVertex', 'zVertex', use_32bit )
        # compine the two into a single tuple
        vertices = (np.append(verticesC[0],verticesV[0]),
                    np.append(verticesC[1],verticesV[1]),
                    np.append(verticesC[2],verticesV[2]))
        del verticesC, verticesV

        # Build cell list
        (connectivity, offsets, valid_mask) = utils.build_edge_cell_lists( mesh_file, args.blocking )

        if not separate_mesh_file:
            mesh_file.close()
            mesh_file = None

        build_field_time_series( time_indices, time_file_names, mesh_file, args.out_dir, args.blocking,
                                 all_dim_vals, 'nEdges', edgeVars, vertices, connectivity, offsets,
                                 valid_mask, use_32bit, args.combine_output, args.append )

        if separate_mesh_file:
            mesh_file.close()


# vim: set expandtab:
