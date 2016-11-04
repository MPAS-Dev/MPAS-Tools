#!/usr/bin/env python

'''
This script takes a mesh file (-m flag) and a file with MOC regions masks
(-f flag) produce by the MPAS mask creator.  The script produces a copy of
the contents of the MOC mask file, adding transects that mark the southern
boundary of each region in a file indicated with the -o flag.  The transect
is applied only to vertices and edges, not cells, because the need for southern
boundary transect data on cells is not foreseen.

Author: Xylar Asay-Davis
last modified: 11/02/2016
'''


import xarray
import argparse
import numpy


def extractSouthernBounary(mesh, moc, latBuffer):
    # Extrcts the southern boundary of each region mask in moc.  Mesh info
    # is taken from mesh.  latBuffer is a number of radians above the southern-
    # most point that should be considered to definitely be in the southern
    # boundary.

    def getEdgeSequenceOnBoundary(startEdge, isBoundaryEdge):
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

        boundaryEdgesOnEdge = -numpy.ones((nEdges, 2), int)

        boundaryEdges = numpy.arange(nEdges)[isBoundaryEdge]
        nBoundaryEdges = len(boundaryEdges)

        # Find the edges on vertex of the vertices on each boundary edge.
        # Each boundary edge must have valid vertices, so none should be out
        # of bounds.
        edgesOnVerticesOnBoundaryEdge = \
            edgesOnVertex[verticesOnEdge[boundaryEdges, :], :]

        # The (typically 3) edges on each vertex of a boundary edge
        # will be the edge itself, another boundary edge and 1 or more
        # non-boundary edges.  We want only the other boundary edge

        # other edge not be this edge
        mask = numpy.not_equal(edgesOnVerticesOnBoundaryEdge,
                               boundaryEdges.reshape((nBoundaryEdges, 1, 1)))

        # other edge must be in range
        mask = numpy.logical_and(mask, edgesOnVerticesOnBoundaryEdge >= 0)
        mask = numpy.logical_and(mask, edgesOnVerticesOnBoundaryEdge < nEdges)

        # other edge must be a boundary edge
        otherEdgeMask = mask.copy()
        otherEdgeMask[mask] = \
            isBoundaryEdge[edgesOnVerticesOnBoundaryEdge[mask]]

        # otherEdgeMask should have exactly one non-zero entry per vertex
        assert(numpy.all(numpy.equal(numpy.sum(numpy.array(otherEdgeMask, int),
                                               axis=2), 1)))

        (edgeIndices, voeIndices, eovIndices) = numpy.nonzero(otherEdgeMask)

        boundaryEdgesOnEdge = -numpy.ones((nEdges, 2), int)
        boundaryEdgesOnEdge[boundaryEdges[edgeIndices], voeIndices] = \
            edgesOnVerticesOnBoundaryEdge[edgeIndices, voeIndices, eovIndices]

        iEdge = startEdge
        edgeSequence = []
        edgeSigns = []
        vertexSequence = []
        signs = (1, -1)
        vertexOnEdgeIndex = 1
        nextEdge = boundaryEdgesOnEdge[iEdge, vertexOnEdgeIndex]
        while True:
            edgeSequence.append(iEdge)
            edgeSigns.append(signs[vertexOnEdgeIndex])
            vertexSequence.append(verticesOnEdge[iEdge, vertexOnEdgeIndex])

            # a trick to determine which is the next vertex and edge to follow
            vertexOnEdgeIndex = int(boundaryEdgesOnEdge[nextEdge, 0] == iEdge)

            iEdge = nextEdge
            nextEdge = boundaryEdgesOnEdge[nextEdge, vertexOnEdgeIndex]
            if iEdge == startEdge:
                break

        edgeSequence = numpy.array(edgeSequence)
        edgeSigns = numpy.array(edgeSigns)
        vertexSequence = numpy.array(vertexSequence)

        return (edgeSequence, edgeSigns, vertexSequence)

    southernBoundaryEdges = []
    southernBounderyEdgeSigns = []
    southernBoundaryVertices = []
    nCells = mesh.dims['nCells']
    nEdges = mesh.dims['nEdges']

    nRegions = moc.dims['nRegions']
    assert(moc.dims['nCells'] == nCells)

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
        cellMask = moc.variables['regionCellMasks'][:, iRegion].values

        # land cells are outside not in the MOC region
        cellsOnEdgeMask = numpy.zeros(cellsOnEdge.shape, bool)
        # set mask values for cells that are in range (not land)
        cellsOnEdgeMask[cellsOnEdgeInRange] = \
            cellMask[cellsOnEdge[cellsOnEdgeInRange]] == 1

        isMOCBoundaryEdge = (cellsOnEdgeMask[:, 0] != cellsOnEdgeMask[:, 1])
        edgesMOCBoundary = numpy.arange(nEdges)[isMOCBoundaryEdge]

        startEdge = numpy.argmin(latEdge[isMOCBoundaryEdge])
        startEdge = edgesMOCBoundary[startEdge]
        minLat = latEdge[startEdge]

        # follow the boundary from this point to get a loop of edges
        # Note: it is possible but unlikely that the southern-most point is
        # not within bulk region of the MOC mask if the region is not a single
        # shape
        edgeSequence, edgeSigns, vertexSequence = \
            getEdgeSequenceOnBoundary(startEdge, isMOCBoundaryEdge)

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
            southernBounderyEdgeSigns.append(edgeSigns)
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
        southernBounderyEdgeSigns.append(edgeSigns[indices])

        # we want one extra vertex in the vertex sequence
        indices = numpy.arange(endIndices[longest],
                               startIndices[longest] + len(edgeSequence) + 1)
        indices = numpy.mod(indices, len(edgeSequence))

        southernBoundaryVertices.append(vertexSequence[indices])

    return (southernBoundaryEdges, southernBounderyEdgeSigns,
            southernBoundaryVertices)


def addTransectsToMOC(mesh, moc, southernBoundaryEdges,
                      southernBounderyEdgeSigns, southernBoundaryVertices):
    # Creates transect fields in moc from the edges, edge signs and vertices
    # defining the southern boundaries.  Mesh info (nEdges and nVertices) is
    # taken from the mesh file.

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

    moc['transectEdgeMasks'] = (('nEdges', 'nTransects'), transectEdgeMasks)
    moc['transectEdgeMaskSigns'] = (('nEdges', 'nTransects'),
                                    transectEdgeMaskSigns)
    moc['transectEdgeGlobalIDs'] = (('nTransects', 'maxEdgesInTransect'),
                                    transectEdgeGlobalIDs)

    moc['transectVertexMasks'] = (('nVertices', 'nTransects'),
                                  transectVertexMasks)
    moc['transectVertexGlobalIDs'] = (('nTransects', 'maxVerticesInTransect'),
                                      transectVertexGlobalIDs)


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

    moc = xarray.open_dataset(args.in_file)
    mesh = xarray.open_dataset(args.mesh_file)

    southernBoundaryEdges, southernBounderyEdgeSigns, \
        southernBoundaryVertices = \
        extractSouthernBounary(mesh, moc, latBuffer=3.*numpy.pi/180.)

    addTransectsToMOC(mesh, moc, southernBoundaryEdges,
                      southernBounderyEdgeSigns,
                      southernBoundaryVertices)

    moc.to_netcdf(args.out_file)
