#!/usr/bin/env python

from __future__ import absolute_import, division, print_function, \
    unicode_literals

from optparse import OptionParser

import xarray

from mpas_tools.io import write_netcdf


def translate(mesh, xOffset=0., yOffset=0.):
    '''
    Translates the coordinate system of the planar MPAS mesh by an arbirary
    shift in x and/or y

    Parameters
    ----------
    mesh : ``xarray.Dataset``
        A planar mesh to translate

    xOffset : float, optional
        user-specified shift in the x-direction

    yOffset : float, optional
        user-specified shift in the y-direction

    '''

    mesh.xCell[:] += xOffset
    mesh.yCell[:] += yOffset
    mesh.xVertex[:] += xOffset
    mesh.yVertex[:] += yOffset
    mesh.xEdge[:] += xOffset
    mesh.yEdge[:] += yOffset


def center_on_mesh(mesh, otherMesh):
    '''
    Translates the coordinate system of the planar MPAS mesh by shifting the
    origin to the center of the domain described in a separate mesh

    Parameters
    ----------
    mesh : ``xarray.Dataset``
        A planar mesh to translate

    otherMesh : ``xarray.Dataset``
        Another planar mesh whose center will become the center of this mesh.
        Uses xCell,yCell or, if those fields do not exist, will secondly try
        x1,y1 fields
    '''

    mpasXcenter, mpasYcenter = get_center(mesh)

    if 'xCell' in otherMesh and 'yCell' in otherMesh:
        dataXcenter, dataYcenter = get_center(otherMesh, xVar='xCell',
                                              yVar='yCell')
    elif 'x1' in otherMesh and 'y1' in otherMesh:
        dataXcenter, dataYcenter = get_center(otherMesh, xVar='x1', yVar='y1')
    else:
        raise ValueError('reference mesh has neither xCell/yCell nor  x1/y1 '
                         'fields.')

    translate(mesh, dataXcenter-mpasXcenter, dataYcenter-mpasYcenter)


def center(mesh):
    '''
    Translates the coordinate system of the planar MPAS mesh by shifting the
    origin to the center of the domain

    Parameters
    ----------
    mesh : ``xarray.Dataset``
        A planar mesh to translate
    '''
    mpasXcenter, mpasYcenter = get_center(mesh)

    translate(mesh, -mpasXcenter, -mpasYcenter)


def get_center(mesh, xVar='xCell', yVar='yCell'):
    '''
    Find the center of the mesh
    '''

    xCenter = (mesh[xVar].min() + mesh[xVar].max()) * 0.5
    yCenter = (mesh[yVar].min() + mesh[yVar].max()) * 0.5

    return xCenter, yCenter


def main():

    print("== Gathering information.  (Invoke with --help for more details. "
          "All arguments are optional)")
    parser = OptionParser()
    parser.description = \
        "This script translates the coordinate system of the planar MPAS " \
        "mesh specified with the -f flag. \n" \
        "There are 3 possible methods to choose from:\n" \
        "1) shift the origin to the center of the domain\n" \
        "2) arbirary shift in x and/or y\n" \
        "3) shift to the center of the domain described in a separate file\n"
    parser.add_option("-f", "--file", dest="fileInName",
                      help="MPAS planar grid file name.", default="grid.nc",
                      metavar="FILENAME")
    parser.add_option("-d", "--datafile", dest="dataFileName",
                      help="data file name to which to match the domain "
                           "center of.  Uses xCell,yCell or, if those fields "
                           "do not exist, will secondly try x1,y1 fields.",
                      metavar="FILENAME")
    parser.add_option("-x", dest="xshift",
                      help="user-specified shift in the x-direction.",
                      type="float", default=0.0, metavar="SHIFT_VALUE")
    parser.add_option("-y", dest="yshift",
                      help="user-specified shift in the y-direction.",
                      type="float", default=0.0, metavar="SHIFT_VALUE")
    parser.add_option("-c", dest="center",
                      help="shift so origin is at center of domain",
                      action="store_true", default=False)
    for option in parser.option_list:
        if option.default != ("NO", "DEFAULT"):
            option.help += (" " if option.help else "") + "[default: %default]"
    options, args = parser.parse_args()

    print("Attempting to translate coordinates in file: {}".format(
            options.fileInName))

    if options.dataFileName is not None and \
            (options.xshift != 0. or options.yshift != 0.):
        raise ValueError('Specifying a datafile AND one or both of x/y shift '
                         'is invalid.  Please select one of those methods '
                         'only.')

    if options.center and (options.xshift != 0. or options.yshift != 0.):
        raise ValueError('Specifying a shift to center AND one or both of x/y '
                         'shift is invalid.  Please select one of those '
                         'methods only.')

    if options.dataFileName is not None and options.center:
        raise ValueError('Specifying a datafile AND a shift to center is '
                         'invalid.  Please select one of those methods only.')

    if not options.center and (options.xshift == 0.) and \
            (options.yshift == 0.) and options.dataFileName is None:
        raise ValueError('No translation method was specified.  Please select '
                         'one.  Run with -h for more information.')

    mesh = xarray.open_dataset(options.fileInName)
    if options.dataFileName is not None:
        print("  Translating coordinates in {} so the domain center matches "
              "the domain center in {}.\n\n".format(options.fileInName,
                                                    options.dataFileName))
        otherMesh = xarray.open_dataset(options.dataFileName)
        center_on_mesh(mesh, otherMesh)

    if options.xshift != 0. or options.yshift != 0.:
        print("  Translating coordinates in {} by user-specified values.  "
              "X-shift={:f}; Y-shift={:f}\n\n".format(options.fileInName,
                                                      options.xshift,
                                                      options.yshift))

        translate(mesh, options.xshift, options.yshift)

    if options.center:
        print("  Translating coordinates in %s so the origin is the center of "
              "the domain.\n\n")

        center(mesh)

    # close the file so we can re-open it for writing
    mesh.close()
    write_netcdf(mesh, options.fileInName)

    print("Translation completed.")


if __name__ == '__main__':
    main()
