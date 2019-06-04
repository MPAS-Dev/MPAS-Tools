#!/usr/bin/env python
"""
Tool to merge 2 MPAS non-contiguous meshes together into a single file
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

    parser.add_argument('infile1', metavar='FILENAME1',
                        help='File name for first mesh to merge')

    parser.add_argument('infile2', metavar='FILENAME2',
                        help='File name for second mesh to merge')

    parser.add_argument('-o', dest='outfile', default='merged_mesh.nc', metavar='FILENAME',
                        help='The merged mesh file')

    return parser.parse_args(args)


def merge_grids(infile1=None, infile2=None, outfile=None, runner=None):
    """
    Merges two MPAS non-contiguous meshes together into a single file

    Parameters
    ----------
    infile1 : str
        The file name for the first mesh to merge

    infile2 : str
        The file name for the second mesh to merge

    outfile : str
        The file name for the first mesh to merge

    runner : str, optional
        The command to write into the global history attribute of the outfile
    """
    now = datetime.now().strftime("%a %b %d %H:%M:%S %Y")
    if not runner:
        runner = '{}.merge_grids(infile1={}, infile2={}, outfile={})'.format(
                os.path.splitext(__file__)[0], infile1, infile2, outfile)

    print('Opening files to merge:\n    {}\n    {}'.format(infile1, infile2))
    print('Creating the merged mesh file: {}'.format(outfile))
    with Dataset(infile1) as nc_in1, Dataset(infile2) as nc_in2, \
            Dataset(outfile, 'w', format="NETCDF3_CLASSIC") as mesh:
        nCells1 = nc_in1.dimensions['nCells'].size
        nEdges1 = nc_in1.dimensions['nEdges'].size
        nVertices1 = nc_in1.dimensions['nVertices'].size

        nCells2 = nc_in2.dimensions['nCells'].size
        nEdges2 = nc_in2.dimensions['nEdges'].size
        nVertices2 = nc_in2.dimensions['nVertices'].size

        if nc_in1.dimensions['vertexDegree'].size != nc_in2.dimensions['vertexDegree'].size:
            raise ValueError("ERROR: The two files have different lengths of the "
                             "vertexDegree dimension.")

        mesh.createDimension('nCells', nCells1 + nCells2)
        mesh.createDimension('nEdges', nEdges1 + nEdges2)
        mesh.createDimension('nVertices', nVertices1 + nVertices2)
        mesh.createDimension('TWO', 2)
        mesh.createDimension('vertexDegree', nc_in1.dimensions['vertexDegree'].size)
        if 'StrLen' in nc_in1.dimensions:
            mesh.createDimension('StrLen', nc_in1.dimensions['StrLen'].size)
        maxEdges = max(nc_in1.dimensions['maxEdges'].size, nc_in2.dimensions['maxEdges'].size)
        mesh.createDimension('maxEdges', maxEdges)
        mesh.createDimension('maxEdges2', maxEdges * 2)

        optionalDims = ('Time', 'nVertLevels', 'nVertInterfaces')
        for dim in optionalDims:
            if dim in nc_in1.dimensions and dim in nc_in2.dimensions:
                if len(nc_in1.dimensions[dim]) != len(nc_in2.dimensions[dim]):
                    raise ValueError("ERROR: The two files have different lengths "
                                     "of the {} dimension.".format(dim))
                if dim == 'Time':
                    mesh.createDimension('Time', size=None)  # make unlimited dimension
                else:
                    mesh.createDimension(dim, nc_in1.dimensions[dim].size)

        print('Merging variable:')
        vars1 = set(nc_in1.variables)
        vars2 = set(nc_in2.variables)
        # only copy variables common to both files
        for varname in (vars1 & vars2):
            print('    {}'.format(varname))
            if nc_in1.variables[varname].dimensions \
                    != nc_in2.variables[varname].dimensions:
                raise ValueError("ERROR: Variable {} has different dimensions in "
                                 "the two files.".format(varname))

            theVar = nc_in1.variables[varname]
            newVar = mesh.createVariable(varname, theVar.dtype, theVar.dimensions)
            # (Assuming here that nCells, nEdges, and nVertices are never both in a variable)
            # now assign value
            if 'nCells' in theVar.dimensions:
                tup1 = ()
                tup2 = ()
                tupMerge = ()
                for ind in range(len(theVar.dimensions)):
                    if theVar.dimensions[ind] == 'nCells':
                        tup1 += (slice(0, nCells1),)
                        tup2 += (slice(0, nCells2),)
                        tupMerge += (slice(nCells1, nCells1 + nCells2),)
                    else:
                        tup1 += (slice(None),)
                        tup2 += (slice(None),)
                        tupMerge += (slice(None),)
                newVar[tup1] = nc_in1.variables[varname][tup1]
                newVar[tupMerge] = nc_in2.variables[varname][tup2]
            elif 'nEdges' in theVar.dimensions:
                tup1 = ()
                tup2 = ()
                tupMerge = ()
                for ind in range(len(theVar.dimensions)):
                    if theVar.dimensions[ind] == 'nEdges':
                        tup1 += (slice(0, nEdges1),)
                        tup2 += (slice(0, nEdges2),)
                        tupMerge += (slice(nEdges1, nEdges1 + nEdges2),)
                    else:
                        tup1 += (slice(None),)
                        tup2 += (slice(None),)
                        tupMerge += (slice(None),)
                newVar[tup1] = nc_in1.variables[varname][tup1]
                newVar[tupMerge] = nc_in2.variables[varname][tup2]
            elif 'nVertices' in theVar.dimensions:
                tup1 = ()
                tup2 = ()
                tupMerge = ()
                for ind in range(len(theVar.dimensions)):
                    if theVar.dimensions[ind] == 'nVertices':
                        tup1 += (slice(0, nVertices1),)
                        tup2 += (slice(0, nVertices2),)
                        tupMerge += (slice(nVertices1, nVertices1 + nVertices2),)
                    else:
                        tup1 += (slice(None),)
                        tup2 += (slice(None),)
                        tupMerge += (slice(None),)
                newVar[tup1] = nc_in1.variables[varname][tup1]
                newVar[tupMerge] = nc_in2.variables[varname][tup2]
            else:
                # just take file 1's version
                newVar[:] = theVar[:]

            # Indexes need adjusting:
            if varname == "indexToCellID":
                newVar[nCells1:] += nCells1
            elif varname == "indexToEdgeID":
                newVar[nEdges1:] += nEdges1
            elif varname == "indexToVertexID":
                newVar[nVertices1:] += nVertices1
            elif varname == "cellsOnEdge":
                part2 = newVar[nEdges1:, :]
                part2[part2 > 0] += nCells1
                newVar[nEdges1:, :] = part2
            elif varname == "edgesOnCell":
                part2 = newVar[nCells1:, :]
                part2[part2 > 0] += nEdges1
                newVar[nCells1:, :] = part2
            elif varname == "edgesOnEdge":
                part2 = newVar[nEdges1:, :]
                part2[part2 > 0] += nEdges1
                newVar[nEdges1:, :] = part2
            elif varname == "cellsOnCell":
                part2 = newVar[nCells1:, :]
                part2[part2 > 0] += nCells1
                newVar[nCells1:, :] = part2
            elif varname == "verticesOnCell":
                part2 = newVar[nCells1:, :]
                part2[part2 > 0] += nVertices1
                newVar[nCells1:, :] = part2
            elif varname == "verticesOnEdge":
                part2 = newVar[nEdges1:, :]
                part2[part2 > 0] += nVertices1
                newVar[nEdges1:, :] = part2
            elif varname == "edgesOnVertex":
                part2 = newVar[nVertices1:, :]
                part2[part2 > 0] += nEdges1
                newVar[nVertices1:, :] = part2
            elif varname == "cellsOnVertex":
                part2 = newVar[nVertices1:, :]
                part2[part2 > 0] += nCells1
                newVar[nVertices1:, :] = part2

        attrToCopy = ("on_a_sphere", "sphere_radius", "is_periodic")
        for attr in attrToCopy:
            if attr in nc_in1.ncattrs() and attr in nc_in2.ncattrs():
                if nc_in1.getncattr(attr) == nc_in2.getncattr(attr):
                    mesh.setncattr(attr, nc_in1.getncattr(attr))
                else:
                    print(
                        "Warning: Value for '{0}' global attribute differs between "
                        "input files. '{0}' being skipped.".format(attr))
            else:
                print("Warning: '{0}' global attribute not present in both input "
                      "files. '{0}' being skipped.".format(attr))

        # Add merge info to allow exact splitting later
        mesh.merge_point = json.dumps({'nCells': nCells1,
                                       'nEdges': nEdges1,
                                       'nVertices': nVertices1,
                                       'maxEdges1': nc_in1.dimensions['maxEdges'].size,
                                       'maxEdges2': nc_in2.dimensions['maxEdges'].size
                                       })

        run_command = "{}: {} \n".format(now, runner)
        mesh.history = maybe_encode(run_command)

        print('Merge complete! Output file: {}.'.format(outfile))


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
    merge_grids(**vars(arguments))


if __name__ == '__main__':
    main()
