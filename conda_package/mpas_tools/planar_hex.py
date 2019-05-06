#!/usr/bin/env python

from __future__ import absolute_import, division, print_function, \
    unicode_literals

import numpy
import xarray
import argparse

from mpas_tools.io import write_netcdf


def make_planar_hex_mesh(nx, ny, dc, nonperiodic_x,
                         nonperiodic_y, outFileName=None,
                         compareWithFileName=None,
                         format='NETCDF3_64BIT'):
    '''
    Builds an MPAS periodic, planar hexagonal mesh with the requested
    dimensions, optionally saving it to a file, and returs it as an
    ``xarray.Dataset``.

    Parameters
    ----------
    nx : int
        The number of cells in the x direction

    ny : even int
        The number of cells in the y direction (must be an even number for
        periodicity to work out)

    dc : float
        The distance in meters between adjacent cell centers.

    nonperiodic_x, nonperiodic_y : bool
        is the mesh non-periodic in x and y directions?

    outFileName : str, optional
        The name of a file to save the mesh to.  The mesh is not saved to a
        file if no file name is supplied.

    compareWithFileName : str, optional
        The name of a grid file to compare with to see if they are identical,
        used for testing purposes

    format : {'NETCDF4', 'NETCDF4_CLASSIC', 'NETCDF3_64BIT', 'NETCDF3_CLASSIC'}, optional
        The NetCDF format to use for output

    Returns
    -------
    mesh : ``xarray.Dataset``
        The mesh data set, available for further maniuplation such as culling
        cells or removing periodicity.
    '''

    mesh = initial_setup(nx, ny, dc, nonperiodic_x, nonperiodic_y)
    compute_indices_on_cell(mesh)
    if nonperiodic_x:
        mark_cull_cell_nonperiodic_x(mesh)
    if nonperiodic_y:
        mark_cull_cell_nonperiodic_y(mesh)
    compute_indices_on_edge(mesh)
    compute_indices_on_vertex(mesh)
    compute_weights_on_edge(mesh)
    compute_coordinates(mesh)
    add_one_to_indices(mesh)

    # drop some arrays that aren't stantard for MPAS but were used to compute
    # the hex mesh
    mesh = mesh.drop(['cellIdx', 'cellRow', 'cellCol'])
    mesh.attrs.pop('dc')

    if outFileName is not None:
        write_netcdf(mesh, outFileName, format=format)

    if compareWithFileName is not None:
        # used to make sure results are exactly identical to periodic_hex
        make_diff(mesh, compareWithFileName, 'diff.nc')

    return mesh


def initial_setup(nx, ny, dc, nonperiodic_x, nonperiodic_y):
    '''Setup the dimensions and add placeholders for some index variables'''
    if ny % 2 != 0:
        raise ValueError('ny must be divisible by 2 for the grid\'s '
                         'periodicity to work properly.')

    mesh = xarray.Dataset()

    if nonperiodic_x and nonperiodic_y:
        mesh.attrs['is_periodic'] = 'NO'
    else:
        mesh.attrs['is_periodic'] = 'YES'

    if nonperiodic_x:
        mesh.attrs['x_period'] = 0.
    else:
        mesh.attrs['x_period'] = nx * dc
    if nonperiodic_y:
        mesh.attrs['y_period'] = 0.
    else:
        mesh.attrs['y_period'] = ny * dc * numpy.sqrt(3.) / 2.

    mesh.attrs['dc'] = dc

    mesh.attrs['on_a_sphere'] = 'NO'
    mesh.attrs['sphere_radius'] = 0.

    if nonperiodic_x:
        nx = nx + 2
    if nonperiodic_y:
        ny = ny + 2

    nCells = nx * ny
    nEdges = 3 * nCells
    nVertices = 2 * nCells
    vertexDegree = 3
    maxEdges = 6

    # add some basic arrays to get all the dimensions in place
    indexToCellID = numpy.arange(nCells, dtype='i4')
    indexToEdgeID = numpy.arange(nEdges, dtype='i4')
    indexToVertexID = numpy.arange(nVertices, dtype='i4')

    cellIdx = indexToCellID.reshape(ny, nx)
    cellCol, cellRow = numpy.meshgrid(numpy.arange(nx, dtype='i4'),
                                      numpy.arange(ny, dtype='i4'))

    mesh['cellIdx'] = (('ny', 'nx'), cellIdx)
    mesh['cellRow'] = (('nCells'), cellRow.ravel())
    mesh['cellCol'] = (('nCells'), cellCol.ravel())

    mesh['indexToCellID'] = (('nCells'), indexToCellID)
    mesh['indexToEdgeID'] = (('nEdges'), indexToEdgeID)
    mesh['indexToVertexID'] = (('nVertices'), indexToVertexID)

    mesh['cullCell'] = (('nCells'), numpy.zeros(nCells, 'i4'))

    mesh['nEdgesOnCell'] = (('nCells',), 6 * numpy.ones((nCells,), 'i4'))
    mesh['cellsOnCell'] = (('nCells', 'maxEdges'),
                           numpy.zeros((nCells, maxEdges), 'i4'))
    mesh['edgesOnCell'] = (('nCells', 'maxEdges'),
                           numpy.zeros((nCells, maxEdges), 'i4'))
    mesh['verticesOnCell'] = (('nCells', 'maxEdges'),
                              numpy.zeros((nCells, maxEdges), 'i4'))

    mesh['nEdgesOnEdge'] = (('nEdges',), 10 * numpy.ones((nEdges,), 'i4'))
    mesh['cellsOnEdge'] = (('nEdges', 'TWO'),
                           numpy.zeros((nEdges, 2), 'i4'))
    mesh['edgesOnEdge'] = (('nEdges', 'maxEdges2'),
                           -1 * numpy.ones((nEdges, 2 * maxEdges), 'i4'))
    mesh['verticesOnEdge'] = (('nEdges', 'TWO'),
                              numpy.zeros((nEdges, 2), 'i4'))

    mesh['cellsOnVertex'] = (('nVertices', 'vertexDegree'),
                             numpy.zeros((nVertices, vertexDegree), 'i4'))
    mesh['edgesOnVertex'] = (('nVertices', 'vertexDegree'),
                             numpy.zeros((nVertices, vertexDegree), 'i4'))

    return mesh


def mark_cull_cell_nonperiodic_y(mesh):

    cullCell = mesh.cullCell
    nCells = mesh.sizes['nCells']
    nx = mesh.sizes['nx']
    cullCell[0:nx] = 1
    cullCell[nCells - nx:nCells + 1] = 1


def mark_cull_cell_nonperiodic_x(mesh):

    cullCell = mesh.cullCell
    nCells = mesh.sizes['nCells']
    nx = mesh.sizes['nx']
    cullCell[::nx] = 1
    cullCell[nx - 1:nCells + 1:nx] = 1


def compute_indices_on_cell(mesh):

    cellIdx = mesh.cellIdx
    cellRow = mesh.cellRow
    cellCol = mesh.cellCol

    indexToCellID = mesh.indexToCellID

    nx = mesh.sizes['nx']
    ny = mesh.sizes['ny']

    mx = numpy.mod(cellCol - 1, nx)
    my = numpy.mod(cellRow - 1, ny)
    px = numpy.mod(cellCol + 1, nx)
    py = numpy.mod(cellRow + 1, ny)

    mask = numpy.mod(cellRow, 2) == 0

    cellsOnCell = mesh.cellsOnCell
    cellsOnCell[:, 0] = cellIdx[cellRow, mx]
    cellsOnCell[:, 1] = cellIdx[my, mx].where(mask, cellIdx[my, cellCol])
    cellsOnCell[:, 2] = cellIdx[my, cellCol].where(mask, cellIdx[my, px])
    cellsOnCell[:, 3] = cellIdx[cellRow, px]
    cellsOnCell[:, 4] = cellIdx[py, cellCol].where(mask, cellIdx[py, px])
    cellsOnCell[:, 5] = cellIdx[py, mx].where(mask, cellIdx[py, cellCol])

    edgesOnCell = mesh.edgesOnCell
    edgesOnCell[:, 0] = 3 * indexToCellID
    edgesOnCell[:, 1] = 3 * indexToCellID + 1
    edgesOnCell[:, 2] = 3 * indexToCellID + 2
    edgesOnCell[:, 3] = 3 * cellsOnCell[:, 3]
    edgesOnCell[:, 4] = 3 * cellsOnCell[:, 4] + 1
    edgesOnCell[:, 5] = 3 * cellsOnCell[:, 5] + 2

    verticesOnCell = mesh.verticesOnCell
    verticesOnCell[:, 0] = 2 * indexToCellID
    verticesOnCell[:, 1] = 2 * indexToCellID + 1
    verticesOnCell[:, 2] = 2 * cellsOnCell[:, 2]
    verticesOnCell[:, 3] = 2 * cellsOnCell[:, 3] + 1
    verticesOnCell[:, 4] = 2 * cellsOnCell[:, 3]
    verticesOnCell[:, 5] = 2 * cellsOnCell[:, 4] + 1


def compute_indices_on_edge(mesh):
    edgesOnCell = mesh.edgesOnCell
    verticesOnCell = mesh.verticesOnCell
    indexToCellID = mesh.indexToCellID

    cellsOnEdge = mesh.cellsOnEdge
    for j in range(3):
        cellsOnEdge[edgesOnCell[:, j], 1] = indexToCellID
    for j in range(3, 6):
        cellsOnEdge[edgesOnCell[:, j], 0] = indexToCellID

    verticesOnEdge = mesh.verticesOnEdge
    verticesOnEdge[edgesOnCell[:, 0], 0] = verticesOnCell[:, 1]
    verticesOnEdge[edgesOnCell[:, 0], 1] = verticesOnCell[:, 0]
    verticesOnEdge[edgesOnCell[:, 1], 0] = verticesOnCell[:, 2]
    verticesOnEdge[edgesOnCell[:, 1], 1] = verticesOnCell[:, 1]
    verticesOnEdge[edgesOnCell[:, 2], 0] = verticesOnCell[:, 3]
    verticesOnEdge[edgesOnCell[:, 2], 1] = verticesOnCell[:, 2]

    edgesOnEdge = mesh.edgesOnEdge
    edgesOnEdge[edgesOnCell[:, 3], 0] = edgesOnCell[:, 4]
    edgesOnEdge[edgesOnCell[:, 3], 1] = edgesOnCell[:, 5]
    edgesOnEdge[edgesOnCell[:, 3], 2] = edgesOnCell[:, 0]
    edgesOnEdge[edgesOnCell[:, 3], 3] = edgesOnCell[:, 1]
    edgesOnEdge[edgesOnCell[:, 3], 4] = edgesOnCell[:, 2]

    edgesOnEdge[edgesOnCell[:, 4], 0] = edgesOnCell[:, 5]
    edgesOnEdge[edgesOnCell[:, 4], 1] = edgesOnCell[:, 0]
    edgesOnEdge[edgesOnCell[:, 4], 2] = edgesOnCell[:, 1]
    edgesOnEdge[edgesOnCell[:, 4], 3] = edgesOnCell[:, 2]
    edgesOnEdge[edgesOnCell[:, 4], 4] = edgesOnCell[:, 3]

    edgesOnEdge[edgesOnCell[:, 5], 0] = edgesOnCell[:, 0]
    edgesOnEdge[edgesOnCell[:, 5], 1] = edgesOnCell[:, 1]
    edgesOnEdge[edgesOnCell[:, 5], 2] = edgesOnCell[:, 2]
    edgesOnEdge[edgesOnCell[:, 5], 3] = edgesOnCell[:, 3]
    edgesOnEdge[edgesOnCell[:, 5], 4] = edgesOnCell[:, 4]

    edgesOnEdge[edgesOnCell[:, 0], 5] = edgesOnCell[:, 1]
    edgesOnEdge[edgesOnCell[:, 0], 6] = edgesOnCell[:, 2]
    edgesOnEdge[edgesOnCell[:, 0], 7] = edgesOnCell[:, 3]
    edgesOnEdge[edgesOnCell[:, 0], 8] = edgesOnCell[:, 4]
    edgesOnEdge[edgesOnCell[:, 0], 9] = edgesOnCell[:, 5]

    edgesOnEdge[edgesOnCell[:, 1], 5] = edgesOnCell[:, 2]
    edgesOnEdge[edgesOnCell[:, 1], 6] = edgesOnCell[:, 3]
    edgesOnEdge[edgesOnCell[:, 1], 7] = edgesOnCell[:, 4]
    edgesOnEdge[edgesOnCell[:, 1], 8] = edgesOnCell[:, 5]
    edgesOnEdge[edgesOnCell[:, 1], 9] = edgesOnCell[:, 0]

    edgesOnEdge[edgesOnCell[:, 2], 5] = edgesOnCell[:, 3]
    edgesOnEdge[edgesOnCell[:, 2], 6] = edgesOnCell[:, 4]
    edgesOnEdge[edgesOnCell[:, 2], 7] = edgesOnCell[:, 5]
    edgesOnEdge[edgesOnCell[:, 2], 8] = edgesOnCell[:, 0]
    edgesOnEdge[edgesOnCell[:, 2], 9] = edgesOnCell[:, 1]


def compute_indices_on_vertex(mesh):
    edgesOnCell = mesh.edgesOnCell
    verticesOnCell = mesh.verticesOnCell
    indexToCellID = mesh.indexToCellID

    cellsOnVertex = mesh.cellsOnVertex
    cellsOnVertex[verticesOnCell[:, 1], 2] = indexToCellID
    cellsOnVertex[verticesOnCell[:, 3], 0] = indexToCellID
    cellsOnVertex[verticesOnCell[:, 5], 1] = indexToCellID
    cellsOnVertex[verticesOnCell[:, 0], 0] = indexToCellID
    cellsOnVertex[verticesOnCell[:, 2], 1] = indexToCellID
    cellsOnVertex[verticesOnCell[:, 4], 2] = indexToCellID

    edgesOnVertex = mesh.edgesOnVertex
    edgesOnVertex[verticesOnCell[:, 0], 0] = edgesOnCell[:, 0]
    edgesOnVertex[verticesOnCell[:, 1], 0] = edgesOnCell[:, 0]
    edgesOnVertex[verticesOnCell[:, 2], 2] = edgesOnCell[:, 1]
    edgesOnVertex[verticesOnCell[:, 1], 2] = edgesOnCell[:, 1]
    edgesOnVertex[verticesOnCell[:, 2], 1] = edgesOnCell[:, 2]
    edgesOnVertex[verticesOnCell[:, 3], 1] = edgesOnCell[:, 2]


def compute_weights_on_edge(mesh):
    edgesOnCell = mesh.edgesOnCell

    nEdges = mesh.sizes['nEdges']
    maxEdges2 = mesh.sizes['maxEdges2']
    mesh['weightsOnEdge'] = (('nEdges', 'maxEdges2'),
                             numpy.zeros((nEdges, maxEdges2), 'f8'))
    weightsOnEdge = mesh.weightsOnEdge

    weights = (1. / numpy.sqrt(3.)) * numpy.array(
        [[1. / 3., 1. / 6., 0., 1. / 6., 1. / 3.],
         [1. / 3., -1. / 6., 0., 1. / 6., -1. / 3.],
         [-1. / 3., -1. / 6., 0., -1. / 6., -1. / 3.]])
    for i in range(3):
        for j in range(5):
            weightsOnEdge[edgesOnCell[:, i + 3], j] = weights[i, j]
    for i in range(3):
        for j in range(5):
            weightsOnEdge[edgesOnCell[:, i], j + 5] = weights[i, j]


def compute_coordinates(mesh):

    dc = mesh.attrs['dc']
    edgesOnCell = mesh.edgesOnCell
    verticesOnCell = mesh.verticesOnCell

    nCells = mesh.sizes['nCells']
    nEdges = mesh.sizes['nEdges']
    nVertices = mesh.sizes['nVertices']
    vertexDegree = mesh.sizes['vertexDegree']

    mesh['latCell'] = (('nCells'), numpy.zeros((nCells,), 'f8'))
    mesh['lonCell'] = (('nCells'), numpy.zeros((nCells,), 'f8'))

    mesh['latEdge'] = (('nEdges'), numpy.zeros((nEdges,), 'f8'))
    mesh['lonEdge'] = (('nEdges'), numpy.zeros((nEdges,), 'f8'))

    mesh['latVertex'] = (('nVertices'), numpy.zeros((nVertices,), 'f8'))
    mesh['lonVertex'] = (('nVertices'), numpy.zeros((nVertices,), 'f8'))

    cellRow = mesh.cellRow
    cellCol = mesh.cellCol
    mask = numpy.mod(cellRow, 2) == 0

    mesh['xCell'] = (dc * (cellCol + 0.5)).where(mask, dc * (cellCol + 1))
    mesh['yCell'] = dc * (cellRow + 1) * numpy.sqrt(3.) / 2.
    mesh['zCell'] = (('nCells'), numpy.zeros((nCells,), 'f8'))

    mesh['xEdge'] = (('nEdges'), numpy.zeros((nEdges,), 'f8'))
    mesh['yEdge'] = (('nEdges'), numpy.zeros((nEdges,), 'f8'))
    mesh['zEdge'] = (('nEdges'), numpy.zeros((nEdges,), 'f8'))

    mesh.xEdge[edgesOnCell[:, 0]] = mesh.xCell - 0.5 * dc
    mesh.yEdge[edgesOnCell[:, 0]] = mesh.yCell

    mesh.xEdge[edgesOnCell[:, 1]] = mesh.xCell - \
        0.5 * dc * numpy.cos(numpy.pi / 3.)
    mesh.yEdge[edgesOnCell[:, 1]] = mesh.yCell - \
        0.5 * dc * numpy.sin(numpy.pi / 3.)

    mesh.xEdge[edgesOnCell[:, 2]] = mesh.xCell + \
        0.5 * dc * numpy.cos(numpy.pi / 3.)
    mesh.yEdge[edgesOnCell[:, 2]] = mesh.yCell - \
        0.5 * dc * numpy.sin(numpy.pi / 3.)

    mesh['xVertex'] = (('nVertices'), numpy.zeros((nVertices,), 'f8'))
    mesh['yVertex'] = (('nVertices'), numpy.zeros((nVertices,), 'f8'))
    mesh['zVertex'] = (('nVertices'), numpy.zeros((nVertices,), 'f8'))

    mesh.xVertex[verticesOnCell[:, 0]] = mesh.xCell - 0.5 * dc
    mesh.yVertex[verticesOnCell[:, 0]] = mesh.yCell + dc * numpy.sqrt(3.) / 6.

    mesh.xVertex[verticesOnCell[:, 1]] = mesh.xCell - 0.5 * dc
    mesh.yVertex[verticesOnCell[:, 1]] = mesh.yCell - dc * numpy.sqrt(3.) / 6.

    mesh['angleEdge'] = (('nEdges'), numpy.zeros((nEdges,), 'f8'))
    mesh.angleEdge[edgesOnCell[:, 1]] = numpy.pi / 3.
    mesh.angleEdge[edgesOnCell[:, 2]] = 2. * numpy.pi / 3.

    mesh['dcEdge'] = (('nEdges'), dc * numpy.ones((nEdges,), 'f8'))
    mesh['dvEdge'] = mesh.dcEdge * numpy.sqrt(3.) / 3.

    mesh['areaCell'] = \
        (('nCells'), dc**2 * numpy.sqrt(3.) / 2. * numpy.ones((nCells,), 'f8'))

    mesh['areaTriangle'] = \
        (('nVertices'), dc**2 * numpy.sqrt(3.) /
         4. * numpy.ones((nVertices,), 'f8'))

    mesh['kiteAreasOnVertex'] = \
        (('nVertices', 'vertexDegree'),
         dc**2 * numpy.sqrt(3.) / 12. * numpy.ones((nVertices, vertexDegree),
         'f8'))

    mesh['meshDensity'] = (('nCells',), numpy.ones((nCells,), 'f8'))


def add_one_to_indices(mesh):
    '''Neede to adhere to Fortran indexing'''
    indexVars = ['indexToCellID', 'indexToEdgeID', 'indexToVertexID',
                 'cellsOnCell', 'edgesOnCell', 'verticesOnCell',
                 'cellsOnEdge', 'edgesOnEdge', 'verticesOnEdge',
                 'cellsOnVertex', 'edgesOnVertex']
    for var in indexVars:
        mesh[var] = mesh[var] + 1


def make_diff(mesh, refMeshFileName, diffFileName):

    refMesh = xarray.open_dataset(refMeshFileName)
    diff = xarray.Dataset()
    for variable in mesh.data_vars:
        if variable in refMesh:
            diff[variable] = mesh[variable] - refMesh[variable]
            print(diff[variable].name, float(numpy.abs(diff[variable]).max()))
        else:
            print('mesh has extra variable {}'.format(mesh[variable].name))

    for variable in refMesh.data_vars:
        if variable not in mesh:
            print('mesh mising variable {}'.format(refMesh[variable].name))

    for attr in refMesh.attrs:
        if attr not in mesh.attrs:
            print('mesh mising attribute {}'.format(attr))

    for attr in mesh.attrs:
        if attr not in refMesh.attrs:
            print('mesh has extra attribute {}'.format(attr))

    write_netcdf(diff, diffFileName)


def main():

    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--nx', dest='nx', type=int, required=True,
                        help='Cells in x direction')
    parser.add_argument('--ny', dest='ny', type=int, required=True,
                        help='Cells in y direction')
    parser.add_argument('--dc', dest='dc', type=float, required=True,
                        help='Distance between cell centers in meters')
    parser.add_argument('--npx', '--nonperiodic_x', dest='nonperiodic_x',
                        action="store_true",
                        help='non-periodic in x direction')
    parser.add_argument('--npy', '--nonperiodic_y', dest='nonperiodic_y',
                        action="store_true",
                        help='non-periodic in y direction')
    parser.add_argument('-o', '--outFileName', dest='outFileName', type=str,
                        required=False, default='grid.nc',
                        help='The name of the output file')

    args = parser.parse_args()

    make_planar_hex_mesh(args.nx, args.ny, args.dc,
                         args.nonperiodic_x, args.nonperiodic_y,
                         args.outFileName)


if __name__ == '__main__':
    main()
