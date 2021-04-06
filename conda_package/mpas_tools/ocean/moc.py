from __future__ import absolute_import, division, print_function, \
    unicode_literals

import xarray
import numpy
import logging
import sys
from geometric_features.aggregation.ocean import moc

import mpas_tools.mesh.conversion
from mpas_tools.io import write_netcdf


def make_moc_basins_and_transects(gf, mesh_filename,
                                  mask_and_transect_filename,
                                  geojson_filename=None,
                                  mask_filename=None,
                                  logger=None,
                                  dir=None):
    """
    Builds features defining the ocean basins and southern transects used in
    computing the meridional overturning circulation (MOC)

    Parameters
    ----------
    gf : geometric_features.GeometricFeatures
        An object that knows how to download and read geometric features

    mesh_filename : str
        A file with MPAS mesh information

    mask_and_transect_filename : str
        A file to write the MOC region masks and southern-boundary transects to

    geojson_filename : str, optional
        A file to write MOC regions to

    mask_filename : str, optional
        A file to write MOC region masks to

    logger : ``logging.Logger``, optional
        A logger for the output if not stdout

    dir : str, optional
        A directory in which a temporary directory will be added with files
        produced during conversion and then deleted upon completion.

    Returns
    -------
    fc : geometric_features.FeatureCollection
        The new feature collection
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    fcMOC = moc(gf)

    if geojson_filename is not None:
        fcMOC.to_geojson(geojson_filename)

    dsMesh = xarray.open_dataset(mesh_filename)
    dsMasks = mpas_tools.mesh.conversion.mask(dsMesh=dsMesh, fcMask=fcMOC,
                                              logger=logger, dir=dir)

    if mask_filename is not None:
        write_netcdf(dsMasks, mask_filename, char_dim_name='StrLen')

    dsMasksAndTransects = add_moc_southern_boundary_transects(dsMasks, dsMesh,
                                                              logger=logger)
    write_netcdf(dsMasksAndTransects, mask_and_transect_filename,
                 char_dim_name='StrLen')


def add_moc_southern_boundary_transects(dsMask, dsMesh, logger=None):
    """
    Parameters
    ----------
    dsMask : ``xarray.Dataset``
        Region masks defining MOC basins

    dsMesh : ``xarray.Dataset``, optional
        An MPAS mesh on which the masks should be created

    logger : ``logging.Logger``, optional
        A logger for the output if not stdout

    Returns
    -------
    dsMask : ``xarray.Dataset``
        Region masks defining MOC basins and the corresponding southern-boundary
        transects
    """

    useStdout = logger is None
    if useStdout:
        logger = logging.getLogger()
        logger.addHandler(logging.StreamHandler(sys.stdout))
        logger.setLevel(logging.INFO)

    southernBoundaryEdges, southernBoundaryEdgeSigns, \
        southernBoundaryVertices = \
        _extract_southern_boundary(dsMesh, dsMask, latBuffer=3.*numpy.pi/180.,
                                   logger=logger)

    _add_transects_to_moc(dsMesh, dsMask, southernBoundaryEdges,
                          southernBoundaryEdgeSigns,
                          southernBoundaryVertices)

    if useStdout:
        logger.handlers = []

    return dsMask


def _extract_southern_boundary(mesh, mocMask, latBuffer, logger):
    """
    Extracts the southern boundary of each region mask in mocMask.  Mesh info
    is taken from mesh.  latBuffer is a number of radians above the southern-
    most point that should be considered to definitely be in the southern
    boundary.
    """

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
    southernBoundaryEdgeSigns = []
    southernBoundaryVertices = []

    for iRegion in range(nRegions):
        name = mocMask.regionNames[iRegion].values.astype('U')
        logger.info(name)
        cellMask = mocMask.variables['regionCellMasks'][:, iRegion].values

        # land cells are outside not in the MOC region
        cellsOnEdgeMask = numpy.zeros(cellsOnEdge.shape, bool)
        # set mask values for cells that are in range (not land)
        cellsOnEdgeMask[cellsOnEdgeInRange] = \
            cellMask[cellsOnEdge[cellsOnEdgeInRange]] == 1

        logger.info('  computing edge sign...')
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
        logger.info('  done.')

        startEdge = numpy.argmin(latEdge[isMOCBoundaryEdge])
        startEdge = edgesMOCBoundary[startEdge]
        minLat = latEdge[startEdge]

        logger.info('  getting edge sequence...')
        # follow the boundary from this point to get a loop of edges
        # Note: it is possible but unlikely that the southern-most point is
        # not within bulk region of the MOC mask if the region is not a single
        # shape
        edgeSequence, edgeSequenceSigns, vertexSequence = \
            _get_edge_sequence_on_boundary(startEdge, edgeSign, edgesOnVertex,
                                           verticesOnEdge)

        logger.info('  done: {} edges in transect.'.format(len(edgeSequence)))

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
            southernBoundaryEdgeSigns.append(edgeSequenceSigns)
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
        southernBoundaryEdgeSigns.append(edgeSequenceSigns[indices])

        # we want one extra vertex in the vertex sequence
        indices = numpy.arange(endIndices[longest],
                               startIndices[longest] + len(edgeSequence) + 1)
        indices = numpy.mod(indices, len(edgeSequence))

        southernBoundaryVertices.append(vertexSequence[indices])

    return southernBoundaryEdges, southernBoundaryEdgeSigns, \
        southernBoundaryVertices


def _add_transects_to_moc(mesh, mocMask, southernBoundaryEdges,
                          southernBoiundaryEdgeSigns, southernBoundaryVertices):
    """
    Creates transect fields in mocMask from the edges, edge signs and
    vertices defining the southern boundaries.  Mesh info (nEdges and
    nVertices) is taken from the mesh file.
    """

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
            = southernBoiundaryEdgeSigns[iTransect]

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

    if 'nRegionsInGroup' not in mocMask:
        nRegions = mocMask.sizes['nRegions']
        nRegionGroups = 2
        nRegionsInGroup = nRegions*numpy.ones(nRegionGroups, dtype=int)
        regionsInGroup = numpy.zeros((nRegionGroups, nRegions), dtype=int)
        regionGroupNames = ['MOCBasinRegionsGroup', 'all']
        regionNames = mocMask.regionNames.values
        nChar = 64
        for index in range(nRegionGroups):
            regionsInGroup[index, :] = numpy.arange(1, nRegions+1)

        mocMask['nRegionsInGroup'] = (('nRegionGroups',), nRegionsInGroup)

        mocMask['regionsInGroup'] = (('nRegionGroups', 'maxRegionsInGroup'),
                                     regionsInGroup)

        mocMask['regionGroupNames'] = \
            (('nRegionGroups',), numpy.zeros((nRegionGroups,),
                                             dtype='|S{}'.format(nChar)))

        for index in range(nRegionGroups):
            mocMask['regionGroupNames'][index] = regionGroupNames[index]

        # we need to make sure the region names use the same string length
        mocMask['regionNames'] = \
            (('nRegions',), numpy.zeros((nRegions,),
                                        dtype='|S{}'.format(nChar)))

        for index in range(nRegions):
            mocMask['regionNames'][index] = regionNames[index]

    mocMask['transectNames'] = mocMask.regionNames.rename(
        {'nRegions': 'nTransects'})

    mocMask['nTransectsInGroup'] = mocMask.nRegionsInGroup.rename(
        {'nRegionGroups': 'nTransectGroups'})

    mocMask['transectsInGroup'] = mocMask.regionsInGroup.rename(
        {'nRegionGroups': 'nTransectGroups',
         'maxRegionsInGroup': 'maxTransectsInGroup'})

    mocMask['transectGroupNames'] = mocMask.regionGroupNames.rename(
        {'nRegionGroups': 'nTransectGroups'})


def _get_edge_sequence_on_boundary(startEdge, edgeSign, edgesOnVertex,
                                   verticesOnEdge):
    """
    Follows the boundary from a starting edge to produce a sequence of
    edges that form a closed loop.

    startEdge is an edge on the boundary that will be both the start and
    end of the loop.

    isBoundaryEdge is a mask that indicates which edges are on the
    boundary

    returns lists of edges, edge signs and vertices
    """

    iEdge = startEdge
    edgeSequence = []
    vertexSequence = []
    while True:
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

    return edgeSequence, edgeSequenceSigns, vertexSequence
