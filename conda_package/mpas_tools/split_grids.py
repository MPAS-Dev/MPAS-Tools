#!/usr/bin/env python
"""
Tool to split 2 previously merged MPAS non-contiguous meshes into separate files.
Typical usage is:
    split_grids.py -1 outfile1.nc -2 outfile2.nc infile
The optional arguments for nCells, nEdges, nVertices, and maxEdges should
generally not be required as this information is saved in the combined mesh file
as global attributes by the merge_grids.py script.
"""

import os
import sys
import json
import argparse

from datetime import datetime

from netCDF4 import Dataset


def parse_args(args=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('infile', metavar='MESHFILE',
                        help='Mesh file to split')

    parser.add_argument('-1', '--outfile1', default='mesh1.nc', metavar='FILENAME',
                        help='File name for first mesh output \n(default: %(default)s)')

    parser.add_argument('-2', '--outfile2', default='mesh2.nc', metavar='FILENAME',
                        help='File name for second mesh output \n(default: %(default)s)')

    parser.add_argument('--nCells', type=int,
                        help='The number of cells in the first mesh \n'
                             '(default: the value specified in MESHFILE global '
                             'attribute merge_point)')

    parser.add_argument('--nEdges', type=int,
                        help='The number of edges in the first mesh \n'
                             '(default: the value specified in MESHFILE global '
                             'attribute merge_point)')

    parser.add_argument('--nVertices', type=int,
                        help='The number of vertices in the first mesh \n'
                             '(default: the value specified in MESHFILE global '
                             'attribute merge_point)')

    parser.add_argument('--maxEdges', type=int, nargs=2, metavar=('MAXEDGES1', 'MAXEDGES2'),
                        help='The number of maxEdges in each mesh \n'
                             '(default: the value specified in MESHFILE global '
                             'attribute merge_point\n      OR: will use MESHFILE '
                             'maxEdges dimension and assume same for both)')

    return parser.parse_args(args)


def split_grids(infile=None, outfile1=None, outfile2=None,
                nCells=None, nEdges=None, nVertices=None, maxEdges=None, runner=None):
    """
    Split two previously merged MPAS non-contiguous meshes together into
    separate files. Typical usage is:

    .. code:: python

        split_grids(infile='infile.nc', outfile1='outfile1.nc', outfile2='outfile2.nc')

    The optional arguments for ``nCells``, ``nEdges``, ``nVertices``, and ``maxEdges``
    should generally not be required as this information sould have been saved in
    ``infiles``'s global attribute ``merge_point`` when created by
    :func:`mpas_tools.merge_grids.merge_grids`.

    Parameters
    ----------
    infile : str
        The file name for the mesh to split

    outfile1 : str
        The file name for the first split mesh

    outfile2 : str
        The file name for the second split mesh

    nCells : int, optional
        The number of cells in the first mesh (default: the value specified in
        infile global attribute merge_point)

    nEdges : int, optional
        The number of edges in the first mesh (default: the value specified in
        infile global attribute merge_point

    nVertices : int, optional
        The number of vertices in the first mesh (default: the value specified in
        infile global attribute merge_point

    maxEdges : list[int, int], optional
        A list of the number of max edges (int) in each mesh (default: the value
        specified in infile global attribute merge_point OR will use infile
        maxEdges dimension and assume same for both)

    runner : str, optional
        The command to write into the global history attribute of the outfile
    """
    now = datetime.now().strftime("%a %b %d %H:%M:%S %Y")
    if not runner:
        runner = '{}.split_grids(infile={}, outfile1={}, outfile2={}, nCells={},' \
                 'nEdges={}, nVertices={})'.format(os.path.splitext(__file__)[0], 
                                                   infile, outfile1, outfile2,
                                                   nCells, nEdges, nVertices)

    merge_point_args_missing = (nCells is None,
                                nEdges is None,
                                nVertices is None)

    print('Opening {} to split'.format(infile))
    with Dataset(infile) as nc_in:
        # NOTE: Because nCells, nEdges, and nVertices are optional arguments and
        #       the previous merge point can be specified in the mesh file, we
        #       need to do some complicated error handling.
        merge_point_in_file = 'merge_point' in nc_in.ncattrs()
        if not merge_point_in_file and any(merge_point_args_missing):
            raise ValueError('ERROR: Previous merge point under specified!\n'
                             '    nCells, nEdges, and nVertices options must all '
                             'be given, or merge_point global attribute must exist '
                             'in {}'.format(infile))
        elif merge_point_in_file and not any(merge_point_args_missing):
            print('Warning: command line arguments are overriding previous merge '
                  'point as specified in {} merge_point global'
                  ' attribute'.format(infile))
        elif merge_point_in_file:
            if not all(merge_point_args_missing):
                print('Warning: nCells, nEdges, and nVertices options must all '
                      'be given to override speification in {} merge_point global '
                      'attribute'.format(infile))
            try:
                mp = json.loads(nc_in.merge_point)
            except ValueError:
                raise ValueError('ERROR: {} merge_point global attribute is not valid JSON.\n'
                                 '    merge_point: {}'.format(infile, nc_in.merge_point))

            mp_keyset = set(mp)
            if {'nCells', 'nEdges', 'nVertices'} <= mp_keyset:
                nCells = mp['nCells']
                nEdges = mp['nEdges']
                nVertices = mp['nVertices']
            else:
                raise ValueError('ERROR: merge_point global attribute of {} must '
                                 'contain nCells, nEdges, and nVertices.\n'
                                 '    merge_point: {}'.format(infile, mp))
            if {'maxEdges1', 'maxEdges2'} <= mp_keyset:
                maxEdges = [mp['maxEdges1'], mp['maxEdges2']]

        print('Creating the mesh files:\n    {}\n    {}'.format(
                outfile1, outfile2))
        with Dataset(outfile1, 'w',  format="NETCDF3_CLASSIC") as mesh1, \
                Dataset(outfile2, 'w',  format="NETCDF3_CLASSIC") as mesh2:
            mesh1.createDimension('nCells', nCells)
            mesh1.createDimension('nEdges', nEdges)
            mesh1.createDimension('nVertices', nVertices)
            mesh1.createDimension('TWO', 2)
            mesh1.createDimension('vertexDegree',
                                  nc_in.dimensions['vertexDegree'].size)

            mesh2.createDimension('nCells', nc_in.dimensions['nCells'].size - nCells)
            mesh2.createDimension('nEdges', nc_in.dimensions['nEdges'].size - nEdges)
            mesh2.createDimension('nVertices', nc_in.dimensions['nVertices'].size - nVertices)
            mesh2.createDimension('TWO', 2)
            mesh2.createDimension('vertexDegree',
                                  nc_in.dimensions['vertexDegree'].size)

            if 'StrLen' in nc_in.dimensions:
                mesh1.createDimension('StrLen', nc_in.dimensions['StrLen'].size)
                mesh2.createDimension('StrLen', nc_in.dimensions['StrLen'].size)

            if maxEdges is None:
                maxEdges = [nc_in.dimensions['maxEdges'].size,
                            nc_in.dimensions['maxEdges'].size]

            mesh1.createDimension('maxEdges', maxEdges[0])
            mesh1.createDimension('maxEdges2', maxEdges[0] * 2)

            mesh2.createDimension('maxEdges', maxEdges[1])
            mesh2.createDimension('maxEdges2', maxEdges[1] * 2)

            mesh1.createDimension('nVertLevels', nc_in.dimensions['nVertLevels'].size)
            mesh1.createDimension('nVertInterfaces', nc_in.dimensions['nVertInterfaces'].size)
            mesh1.createDimension('Time', size=None)  # make unlimited

            mesh2.createDimension('nVertLevels', nc_in.dimensions['nVertLevels'].size)
            mesh2.createDimension('nVertInterfaces', nc_in.dimensions['nVertInterfaces'].size)
            mesh2.createDimension('Time', size=None)  # make unlimited

            print('Splitting variable:')
            for var in nc_in.variables:
                print('    {}'.format(var))
                var_in = nc_in.variables[var]

                var1 = mesh1.createVariable(var, var_in.dtype, var_in.dimensions)
                var2 = mesh2.createVariable(var, var_in.dtype, var_in.dimensions)

                slice1, slice2 = var_slice(var_in.dimensions, nc_in,
                                           nCells, nEdges, nVertices, maxEdges)

                var1[:] = nc_in.variables[var][slice1]
                var2[:] = nc_in.variables[var][slice2]

                # Adjust the indexes
                if var == 'indexToCellID':
                    var2[:] -= nCells
                elif var == 'indexToEdgeID':
                    var2[:] -= nEdges
                elif var == 'indexToVertexID':
                    var2[:] -= nVertices
                elif var in ['cellsOnCell', 'cellsOnEdge', 'cellsOnVertex']:
                    tmp = var2[...]
                    tmp[tmp > 0] -= nCells
                    var2[:] = tmp
                elif var in ['edgesOnCell',  'edgesOnEdge', 'edgesOnVertex']:
                    tmp = var2[...]
                    tmp[tmp > 0] -= nEdges
                    var2[:] = tmp
                elif var in ['verticesOnCell', 'verticesOnEdge']:
                    tmp = var2[...]
                    tmp[tmp > 0] -= nVertices
                    var2[:] = tmp

            attr_to_copy = ("on_a_sphere", "sphere_radius", "is_periodic")
            for attr in attr_to_copy:
                if attr in nc_in.ncattrs():
                    mesh1.setncattr(attr, nc_in.getncattr(attr))
                    mesh2.setncattr(attr, nc_in.getncattr(attr))
                else:
                    print("Warning: '{0}' global attribute not present in input "
                          "file. '{0}' will not be added to the two output "
                          "files.".format(attr))

            run_command = '{}: {} \n'.format(now, runner)
            if 'history' in nc_in.ncattrs():
                mesh1.history = maybe_encode(run_command + nc_in.history)
                mesh2.history = maybe_encode(run_command + nc_in.history)
            else:
                mesh1.history = maybe_encode(run_command)
                mesh2.history = maybe_encode(run_command)

    print('Split complete! Mesh files:\n    {}\n    {}'.format(outfile1, outfile2))


def var_slice(dimensions, nc_in, nCells, nEdges, nVertices, maxEdges):
    slice1 = ()
    slice2 = ()
    for dim in dimensions:
        if dim == 'nCells':
            slice1 += (slice(0, nCells),)
            slice2 += (slice(nCells, nc_in.dimensions['nCells'].size),)
        elif dim == 'nEdges':
            slice1 += (slice(0, nEdges),)
            slice2 += (slice(nEdges, nc_in.dimensions['nEdges'].size),)
        elif dim == 'nVertices':
            slice1 += (slice(0, nVertices),)
            slice2 += (slice(nVertices, nc_in.dimensions['nVertices'].size),)
        elif dim == 'maxEdges':
            slice1 += (slice(0, maxEdges[0]),)
            slice2 += (slice(0, maxEdges[1]),)
        elif dim == 'maxEdges2':
            slice1 += (slice(0, maxEdges[0]*2),)
            slice2 += (slice(0, maxEdges[1]*2),)
        else:
            slice1 += (slice(None),)
            slice2 += (slice(None),)

    return slice1, slice2


# NOTE: Python 2 and 3 string fun conflicting with NC_CHAR vs NC_STRING, see:
#       https://github.com/Unidata/netcdf4-python/issues/529
def maybe_encode(string, encoding='ascii'):
    try:
        return string.encode(encoding)
    except UnicodeEncodeError:
        return string


def main():
    arguments = parse_args()
    arguments.runner = ' '.join(sys.argv[:])
    split_grids(**vars(arguments))


if __name__ == '__main__':
    main()
