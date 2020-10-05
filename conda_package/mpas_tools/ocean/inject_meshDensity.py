#!/usr/bin/env python
# Simple script to inject mesh density onto a mesh
# example usage:
#   ./inject_meshDensity.py cellWidthVsLatLon.nc base_mesh.nc
# where:
#   cellWidthVsLatLon.nc is a netcdf file with cellWidth
#   base_mesh.nc is the mpas netcdf file where meshDensity is added
# Mark Petersen, 7/24/2018

from __future__ import absolute_import, division, print_function, \
    unicode_literals

import numpy as np
import netCDF4 as nc4
import sys

from mpas_tools.mesh.interpolation import interp_bilin


def inject_meshDensity_from_file(cw_filename, mesh_filename, on_sphere=True):
    """
    Add a ``meshDensity`` field into an MPAS mesh.  The mesh density is defined
    as:

      meshDensity = (minCellWidth / cellWidth)**4

    Parameters
    ----------
    cw_filename : str
        The file name to read ``cellWidth`` and coordinates from

    mesh_filename : str
        The mesh file to add ``meshDensity`` to

    on_sphere : bool, optional
        Whether the mesh is spherical (as opposed to planar)
    """
    print('Read cell width field from nc file regular grid...')
    ds = nc4.Dataset(cw_filename,'r')
    cellWidth = ds.variables['cellWidth'][:]
    if on_sphere:
        lon = ds.variables['lon'][:]
        lat = ds.variables['lat'][:]
        ds.close()
        inject_spherical_meshDensity(cellWidth, lon, lat, mesh_filename)
    else:
        x = ds.variables['x'][:]
        y = ds.variables['y'][:]
        ds.close()
        inject_spherical_meshDensity(cellWidth, x, y, mesh_filename)



def inject_spherical_meshDensity(cellWidth, lon, lat, mesh_filename):
    """
    Add a ``meshDensity`` field into a spherical MPAS mesh.  The mesh density is
    defined as:

      meshDensity = (minCellWidth / cellWidth)**4

    Parameters
    ----------
    cellWidth : ndarray
        m x n array of cell width in km

    lon : ndarray
        longitude in degrees (length n and between -180 and 180)

    lat : ndarray
        longitude in degrees (length m and between -90 and 90)

    mesh_filename : str
        The mesh file to add ``meshDensity`` to
    """

    minCellWidth = cellWidth.min()
    meshDensityVsXY = (minCellWidth / cellWidth)**4
    print('  minimum cell width in grid definition: {0:.0f} km'.format(
        minCellWidth))
    print('  maximum cell width in grid definition: {0:.0f} km'.format(
        cellWidth.max()))

    print('Open unstructured MPAS mesh file...')
    ds = nc4.Dataset(mesh_filename, 'r+')
    meshDensity = ds.variables['meshDensity']
    lonCell = ds.variables['lonCell'][:]
    latCell = ds.variables['latCell'][:]

    lonCell = np.mod(np.rad2deg(lonCell) + 180., 360.) - 180.
    latCell = np.rad2deg(latCell)

    print('Interpolating and writing meshDensity...')
    mpasMeshDensity = interp_bilin(lon, lat, meshDensityVsXY, lonCell, latCell)

    meshDensity[:] = mpasMeshDensity

    ds.close()


def inject_planar_meshDensity(cellWidth, x, y, mesh_filename):
    """
    Add a ``meshDensity`` field into a planar MPAS mesh.  The mesh density is
    defined as:

      meshDensity = (minCellWidth / cellWidth)**4

    Parameters
    ----------
    cellWidth : ndarray
        m x n array of cell width in km

    x, y : ndarray
        Planar coordinates in meters

    mesh_filename : str
        The mesh file to add ``meshDensity`` to
    """
    minCellWidth = cellWidth.min()
    meshDensityVsXY = (minCellWidth / cellWidth)**4
    print('  minimum cell width in grid definition: {0:.0f} km'.format(minCellWidth))
    print('  maximum cell width in grid definition: {0:.0f} km'.format(cellWidth.max()))

    print('Open unstructured MPAS mesh file...')
    ds = nc4.Dataset(mesh_filename, 'r+')
    meshDensity = ds.variables['meshDensity']
    xCell = ds.variables['xCell'][:]
    yCell = ds.variables['xCell'][:]

    print('Interpolating and writing meshDensity...')
    mpasMeshDensity = interp_bilin(x, y, meshDensityVsXY, xCell, yCell)

    meshDensity[:] = mpasMeshDensity

    ds.close()


if __name__ == "__main__":

    inject_meshDensity_from_file(cw_filename=sys.argv[1],
                                 mesh_filename=sys.argv[2])
