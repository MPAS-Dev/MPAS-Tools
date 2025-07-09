"""
Efficient interpolation from regular (tensor) grids to MPAS mesh cell centers.

This module provides fast, memory-efficient routines for bilinear interpolation
from a regular grid (e.g., longitude/latitude) to the unstructured cell centers
of an MPAS mesh. It is intended as a replacement for `scipy` interpolation
routines, which are slow and memory-intensive for large meshes.

Main function:
- interp_bilin: Bilinear interpolation from a tensor grid to MPAS mesh cell
  centers.

Note:
    - Input coordinates (x, y) and cell center coordinates (xCell, yCell)
      should be in the same units (typically degrees for lon/lat).
    - No extrapolation is performed: all xCell/yCell must be within the bounds
      of x/y.
"""

import numpy as np


def interp_bilin(x, y, field, xCell, yCell):
    """
    Perform bilinear interpolation of ``field`` on a regular (tensor) grid to
    cell centers on an MPAS mesh.

    This function is designed to be much faster and more memory-efficient than
    using `scipy.interpolate` routines for large MPAS meshes.

    Parameters
    ----------
    x : ndarray
        1D array of x coordinates of the input grid (length n).
        For geographic data, this is typically longitude in degrees.

    y : ndarray
        1D array of y coordinates of the input grid (length m).
        For geographic data, this is typically latitude in degrees.

    field : ndarray
        2D array of field values of shape (m, n), defined on the (y, x) grid.

    xCell : ndarray
        1D array of x coordinates of MPAS cell centers (same units as x).

    yCell : ndarray
        1D array of y coordinates of MPAS cell centers (same units as y).

    Returns
    -------
    mpasField : ndarray
        1D array of interpolated field values at MPAS cell centers.

    Notes
    -----
    - All xCell and yCell values must be within the bounds of x and y.
    - No extrapolation is performed.
    - For longitude/latitude grids, it is recommended to use degrees to avoid
      round-off issues at the poles or dateline.
    - This function is intended for use cases where the input grid is regular
      (tensor product), not scattered points.

    Examples
    --------
    >>> import numpy as np
    >>> from mpas_tools.mesh.interpolation import interp_bilin
    >>> x = np.linspace(-180, 180, 361)
    >>> y = np.linspace(-90, 90, 181)
    >>> field = np.random.rand(181, 361)
    >>> xCell = np.array([0.0, 45.0])
    >>> yCell = np.array([0.0, 45.0])
    >>> values = interp_bilin(x, y, field, xCell, yCell)
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
    xFrac[mask] += 1.0

    mask = yIndices == len(y) - 1
    yIndices[mask] -= 1
    yFrac[mask] += 1.0

    mpasField = (
        (1.0 - xFrac) * (1.0 - yFrac) * field[yIndices, xIndices]
        + xFrac * (1.0 - yFrac) * field[yIndices, xIndices + 1]
        + (1.0 - xFrac) * yFrac * field[yIndices + 1, xIndices]
        + xFrac * yFrac * field[yIndices + 1, xIndices + 1]
    )

    return mpasField
