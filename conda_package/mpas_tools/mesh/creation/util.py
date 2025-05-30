import numpy as np


def circumcenter(on_sphere, x1, y1, z1, x2, y2, z2, x3, y3, z3):
    """
    Compute the circumcenter of the triangle (possibly on a sphere)
    with the three given vertices in Cartesian coordinates.

    Parameters
    ----------
    on_sphere : bool
        If True, the circumcenter is computed on a sphere.
        If False, the circumcenter is computed in Cartesian coordinates.

    x1, y1, z1 : float or numpy.ndarray
        Cartesian coordinates of the first vertex

    x2, y2, z2 : float or numpy.ndarray
        Cartesian coordinates of the second vertex

    x3, y3, z3 : float or numpy.ndarray
        Cartesian coordinates of the third vertex

    Returns
    -------
    xv, yv, zv : float or numpy.ndarray
        The circumcenter(s) of the triangle(s)
    """
    x1 = np.asarray(x1)
    y1 = np.asarray(y1)
    z1 = np.asarray(z1)
    x2 = np.asarray(x2)
    y2 = np.asarray(y2)
    z2 = np.asarray(z2)
    x3 = np.asarray(x3)
    y3 = np.asarray(y3)
    z3 = np.asarray(z3)

    if on_sphere:
        a = (x2 - x3) ** 2 + (y2 - y3) ** 2 + (z2 - z3) ** 2
        b = (x3 - x1) ** 2 + (y3 - y1) ** 2 + (z3 - z1) ** 2
        c = (x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2

        pbc = a * (-a + b + c)
        apc = b * (a - b + c)
        abp = c * (a + b - c)

        denom = pbc + apc + abp
        xv = (pbc * x1 + apc * x2 + abp * x3) / denom
        yv = (pbc * y1 + apc * y2 + abp * y3) / denom
        zv = (pbc * z1 + apc * z2 + abp * z3) / denom
    else:
        d = 2 * (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2))

        xv = (
            (x1**2 + y1**2) * (y2 - y3)
            + (x2**2 + y2**2) * (y3 - y1)
            + (x3**2 + y3**2) * (y1 - y2)
        ) / d
        yv = (
            (x1**2 + y1**2) * (x3 - x2)
            + (x2**2 + y2**2) * (x1 - x3)
            + (x3**2 + y3**2) * (x2 - x1)
        ) / d
        zv = np.zeros_like(xv)

    return xv, yv, zv


def lonlat2xyz(lon, lat, R=6378206.4):
    """
    Convert from longitude and latitude to Cartesian coordinates

    Parameters
    ----------
    lon : float or numpy.ndarray
        longitude
    lat : float or numpy.ndarray
        latitude
    R : float, optional
        Earth radius in meters

    Returns
    -------
    x, y, z: float or numpy.array
        Cartesian coordinates
    """

    lon = np.deg2rad(lon)
    lat = np.deg2rad(lat)
    x = R * np.multiply(np.cos(lon), np.cos(lat))
    y = R * np.multiply(np.sin(lon), np.cos(lat))
    z = R * np.sin(lat)

    return x, y, z
