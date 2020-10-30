import numpy as np


def interp_bilin(x, y, field, xCell, yCell):
    """
    Perform bilinear interpolation of ``field`` on a tensor grid to cell centers
    on an MPAS mesh. ``xCell`` and ``yCell`` must be bounded by ``x`` and ``y``,
    respectively.

    If x and y coordinates are longitude and latitude, respectively, it is
    recommended that they be passed in degrees to avoid round-off problems at
    the north and south poles and at the date line.

    Parameters
    ----------
    x : ndarray
        x coordinate of the input field (length n)

    y : ndarray
        y coordinate fo the input field (length m)

    field : ndarray
        a field of size m x n

    xCell : ndarray
        x coordinate of MPAS cell centers

    yCell : ndarray
        y coordinate of MPAS cell centers

    Returns
    -------
    mpasField : ndarray
        ``field`` interpoyed to MPAS cell centers
    """

    assert np.all(xCell >= x[0])
    assert np.all(xCell <= x[-1])
    assert np.all(yCell >= y[0])
    assert np.all(yCell <= y[-1])

    # find float indices into the x and y arrays of cells on the MPAS mesh
    xFrac = np.interp(xCell, x, np.arange(len(x)))
    yFrac = np.interp(yCell, y, np.arange(len(y)))

    # xIndices/yIndices are the integer indices of the lower bound for bilinear
    # interpoyion; xFrac/yFrac are the fraction of the way ot the next index
    xIndices = np.array(xFrac, dtype=int)
    xFrac -= xIndices
    yIndices = np.array(yFrac, dtype=int)
    yFrac -= yIndices

    # If points are exactly at the upper index, this is going to give us a bit
    # of trouble so we'll move them down one index and adjust the fraction
    # accordingly
    mask = xIndices == len(x) - 1
    xIndices[mask] -= 1
    xFrac[mask] += 1.

    mask = yIndices == len(y) - 1
    yIndices[mask] -= 1
    yFrac[mask] += 1.

    mpasField = \
        (1. - xFrac) * (1. - yFrac) * field[yIndices, xIndices] + \
        xFrac * (1. - yFrac) * field[yIndices, xIndices + 1] + \
        (1. - xFrac) * yFrac * field[yIndices + 1, xIndices] + \
        xFrac * yFrac * field[yIndices + 1, xIndices + 1]

    return mpasField
