from __future__ import absolute_import, division, print_function, \
    unicode_literals

import numpy


def add_critical_land_blockages(dsMask, dsBlockages):
    '''
    Parameters
    ----------
    dsMask : `xarray.Dataset`
        The mask to which critical blockages should be added
    dsBlockage : `xarray.Dataset`
        The transect masks defining critical land regions that should block
        ocean flow (e.g. the Antarctic Peninsula)

    Returns
    -------
    dsMask : `xarray.Dataset`
        The mask with critical blockages included
    '''

    dsMask = dsMask.copy()

    nTransects = dsBlockages.sizes['nTransects']
    for transectIndex in range(nTransects):
        dsMask.regionCellMasks[:, 0] = numpy.maximum(
            dsBlockages.transectCellMasks[:, transectIndex],
            dsMask.regionCellMasks[:, 0])

    return dsMask


def widen_transect_edge_masks(dsMask, dsMesh, latitude_threshold=43.0):
    '''
    Parameters
    ----------
    dsMask : `xarray.Dataset`
        The mask to which critical blockages should be added
    dsMesh : `xarray.Dataset`
        The transect masks defining critical land regions that should block
        ocean flow (e.g. the Antarctic Peninsula)
    latitude_threshold : float
        Minimum latitude, degrees, for transect widening

    Returns
    -------
    dsMask : `xarray.Dataset`
        The mask with critical blockages included
    '''
    latitude_threshold_radians = numpy.deg2rad(latitude_threshold)

    dsMask = dsMask.copy()

    maxEdges = dsMesh.sizes['maxEdges']

    latMask = numpy.abs(dsMesh.latEdge) > latitude_threshold_radians

    edgeMask = numpy.logical_and(
        latMask, dsMask.transectEdgeMasks == 1)
    for iEdge in range(maxEdges):
        eoc = dsMesh.edgesOnCell[:, iEdge]-1
        mask = numpy.logical_and(eoc >= 0,
                                 edgeMask[eoc])
        # cells with a neighboring transect edge should be masked to 1
        dsMask['transectCellMasks'] = dsMask.transectCellMasks.where(
            numpy.logical_not(mask), 1.)

    return dsMask
