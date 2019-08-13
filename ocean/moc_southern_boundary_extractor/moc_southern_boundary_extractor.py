#!/usr/bin/env python

'''
This script takes a mesh file (-m flag) and a file with MOC regions masks
(-f flag) produce by the MPAS mask creator.  The script produces a copy of
the contents of the MOC mask file, adding transects that mark the southern
boundary of each region in a file indicated with the -o flag.  The transect
is applied only to vertices and edges, not cells, because the need for southern
boundary transect data on cells is not foreseen.

Author: Xylar Asay-Davis
last modified: 5/22/2018
'''

from __future__ import absolute_import, division, print_function, \
    unicode_literals

import xarray
import argparse
import numpy


def getEdgeSequenceOnBoundary(startEdge, edgeSign, edgesOnVertex,
                              verticesOnEdge):
    # Follows the boundary from a starting edge to produce a sequence of
    # edges that form a closed loop.
    #
    # startEdge is an edge on the boundary that will be both the start and
    # end of the loop.
    #
    # isBoundaryEdge is a mask that indicates which edges are on the
    # boundary
    #
    # returns lists of edges, edge signs and vertices

    iEdge = startEdge
    edgeSequence = []
    vertexSequence = []
    while(True):
        assert(edgeSign[iEdge] == 1. or edgeSign[iEdge] == -1.)
        if edgeSign[iEdge] == 1.:
            v = 0
        else:
            v = 1
        iVertex = verticesOnEdge[iEdge, v]

        eov = edgesOnVertex[iVertex, :]

        # find the edge that is not iEdge but is on the boundary
        nextEdge = -1
        for edge in eov:
            if edge != iEdge and edgeSign[edge] != 0:
                nextEdge = edge
                break
        assert(nextEdge != -1)

        edgeSequence.append(iEdge)
        vertexSequence.append(iVertex)

        iEdge = nextEdge

        if iEdge == startEdge:
            break

    edgeSequence = numpy.array(edgeSequence)
    edgeSequenceSigns = edgeSign[edgeSequence]
    vertexSequence = numpy.array(vertexSequence)

    return (edgeSequence, edgeSequenceSigns, vertexSequence)


def extractSouthernBounary(mesh, mocMask, latBuffer):
    # Extrcts the southern boundary of each region mask in mocMask.  Mesh info
    # is taken from mesh.  latBuffer is a number of radians above the southern-
    # most point that should be considered to definitely be in the southern
    # boundary.

    southernBoundaryEdges = []
    southernBounderyEdgeSigns = []
    southernBoundaryVertices = []
    nCells = mesh.dims['nCells']
    nEdges = mesh.dims['nEdges']

    nRegions = mocMask.dims['nRegions']
    assert(mocMask.dims['nCells'] == nCells)

    # convert to python zero-based indices
    cellsOnEdge = mesh.variables['cellsOnEdge'].values-1
    verticesOnEdge = mesh.variables['verticesOnEdge'].values-1
    edgesOnVertex = mesh.variables['edgesOnVertex'].values-1

    latEdge = mesh.variables['latEdge'].values

    cellsOnEdgeInRange = numpy.logical_and(cellsOnEdge >= 0,
                                           cellsOnEdge < nCells)

    southernBoundaryEdges = []
    southernBounderyEdgeSigns = []
    southernBoundaryVertices = []

    for iRegion in range(nRegions):
        name = mocMask.regionNames[iRegion].values.astype('U')
        print(name)
        cellMask = mocMask.variables['regionCellMasks'][:, iRegion].values

        # land cells are outside not in the MOC region
        cellsOnEdgeMask = numpy.zeros(cellsOnEdge.shape, bool)
        # set mask values for cells that are in range (not land)
        cellsOnEdgeMask[cellsOnEdgeInRange] = \
            cellMask[cellsOnEdge[cellsOnEdgeInRange]] == 1

        print('  computing edge sign...')
        edgeSign = numpy.zeros(nEdges)
        # positive sign if the first cell on edge is in the region
        mask = numpy.logical_and(cellsOnEdgeMask[:, 0],
                                 numpy.logical_not(cellsOnEdgeMask[:, 1]))
        edgeSign[mask] = -1.
        # negative sign if the second cell on edge is in the region
        mask = numpy.logical_and(cellsOnEdgeMask[:, 1],
                                 numpy.logical_not(cellsOnEdgeMask[:, 0]))
        edgeSign[mask] = 1.
        isMOCBoundaryEdge = edgeSign != 0.
        edgesMOCBoundary = numpy.arange(nEdges)[isMOCBoundaryEdge]
        print('  done.')

        startEdge = numpy.argmin(latEdge[isMOCBoundaryEdge])
        startEdge = edgesMOCBoundary[startEdge]
        minLat = latEdge[startEdge]

        print('  getting edge sequence...')
        # follow the boundary from this point to get a loop of edges
        # Note: it is possible but unlikely that the southern-most point is
        # not within bulk region of the MOC mask if the region is not a single
        # shape
        edgeSequence, edgeSequenceSigns, vertexSequence = \
            getEdgeSequenceOnBoundary(startEdge, edgeSign, edgesOnVertex,
                                      verticesOnEdge)

        print('  done: {} edges in transect.'.format(len(edgeSequence)))

        aboveSouthernBoundary = latEdge[edgeSequence] > minLat + latBuffer

        # find and eliminate the longest contiguous sequence (if any) from the
        # edge sequence that is above the possible region of the soutehrn
        # boundary

        belowToAbove = \
            numpy.logical_and(numpy.logical_not(aboveSouthernBoundary[0:-1]),
                              aboveSouthernBoundary[1:])

        aboveToBelow = \
            numpy.logical_and(aboveSouthernBoundary[0:-1],
                              numpy.logical_not(aboveSouthernBoundary[1:]))

        startIndices = numpy.arange(1, len(edgeSequence))[belowToAbove]
        endIndices = numpy.arange(1, len(edgeSequence))[aboveToBelow]

        assert(len(startIndices) == len(endIndices))

        if len(startIndices) == 0:
            # the whole sequence is the southern boundary
            southernBoundaryEdges.append(edgeSequence)
            southernBounderyEdgeSigns.append(edgeSequenceSigns)
            southernBoundaryVertices.append(vertexSequence)
            continue

        # there are some parts of the sequence above the boundary. Let's
        # eliminate the longest one.

        aboveLength = endIndices - startIndices
        longest = numpy.argmax(aboveLength)
        # we want all the indices in the sequence *not* part of the longest
        indices = numpy.arange(endIndices[longest],
                               startIndices[longest] + len(edgeSequence))
        indices = numpy.mod(indices, len(edgeSequence))

        southernBoundaryEdges.append(edgeSequence[indices])
        southernBounderyEdgeSigns.append(edgeSequenceSigns[indices])

        # we want one extra vertex in the vertex sequence
        indices = numpy.arange(endIndices[longest],
                               startIndices[longest] + len(edgeSequence) + 1)
        indices = numpy.mod(indices, len(edgeSequence))

        southernBoundaryVertices.append(vertexSequence[indices])

    return (southernBoundaryEdges, southernBounderyEdgeSigns,
            southernBoundaryVertices)


def addTransectsToMOC(mesh, mocMask, southernBoundaryEdges,
                      southernBounderyEdgeSigns, southernBoundaryVertices):
    # Creates transect fields in mocMask from the edges, edge signs and
    # vertices defining the southern boundaries.  Mesh info (nEdges and
    # nVertices) is taken from the mesh file.

    nTransects = len(southernBoundaryEdges)

    nEdges = mesh.dims['nEdges']
    nVertices = mesh.dims['nVertices']

    maxEdgesInTransect = numpy.amax([len(southernBoundaryEdges[iTransect])
                                     for iTransect in range(nTransects)])

    maxVerticesInTransect = \
        numpy.amax([len(southernBoundaryVertices[iTransect])
                    for iTransect in range(nTransects)])

    transectEdgeMasks = numpy.zeros((nEdges, nTransects),
                                    numpy.int32)
    transectEdgeMaskSigns = numpy.zeros((nEdges, nTransects),
                                        numpy.int32)
    transectEdgeGlobalIDs = numpy.zeros((nTransects, maxEdgesInTransect),
                                        numpy.int32)
    transectVertexMasks = numpy.zeros((nVertices, nTransects),
                                      numpy.int32)
    transectVertexGlobalIDs = numpy.zeros((nTransects, maxVerticesInTransect),
                                          numpy.int32)

    for iTransect in range(nTransects):
        transectEdgeMasks[southernBoundaryEdges[iTransect], iTransect] = 1

        transectEdgeMaskSigns[southernBoundaryEdges[iTransect], iTransect] \
            = southernBounderyEdgeSigns[iTransect]

        transectCount = len(southernBoundaryEdges[iTransect])
        transectEdgeGlobalIDs[iTransect, 0:transectCount] \
            = southernBoundaryEdges[iTransect] + 1

        transectVertexMasks[southernBoundaryVertices[iTransect], iTransect] = 1

        transectCount = len(southernBoundaryVertices[iTransect])
        transectVertexGlobalIDs[iTransect, 0:transectCount] \
            = southernBoundaryVertices[iTransect] + 1

    mocMask['transectEdgeMasks'] = \
        (('nEdges', 'nTransects'), transectEdgeMasks)
    mocMask['transectEdgeMaskSigns'] = (('nEdges', 'nTransects'),
                                        transectEdgeMaskSigns)
    mocMask['transectEdgeGlobalIDs'] = (('nTransects', 'maxEdgesInTransect'),
                                        transectEdgeGlobalIDs)

    mocMask['transectVertexMasks'] = (('nVertices', 'nTransects'),
                                      transectVertexMasks)
    mocMask['transectVertexGlobalIDs'] = \
        (('nTransects', 'maxVerticesInTransect'), transectVertexGlobalIDs)

    mocMask['transectNames'] = mocMask.regionNames.rename(
        {'nRegions': 'nTransects'})


if __name__ == "__main__":

    parser = \
        argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-f', '--in_file', dest='in_file',
                        help='Input file with MOC masks', metavar='IN_FILE',
                        required=True)
    parser.add_argument('-m', '--mesh_file', dest='mesh_file',
                        help='Input mesh file', metavar='MESH_FILE',
                        required=True)
    parser.add_argument('-o', '--out_file', dest='out_file',
                        help='Output file for MOC masks and southern-boundary '
                        'transects', metavar='OUT_FILE',
                        required=True)
    args = parser.parse_args()

    mocMask = xarray.open_dataset(args.in_file)
    mesh = xarray.open_dataset(args.mesh_file)

    southernBoundaryEdges, southernBounderyEdgeSigns, \
        southernBoundaryVertices = \
        extractSouthernBounary(mesh, mocMask, latBuffer=3.*numpy.pi/180.)

    addTransectsToMOC(mesh, mocMask, southernBoundaryEdges,
                      southernBounderyEdgeSigns,
                      southernBoundaryVertices)

    mocMask.to_netcdf(args.out_file)
