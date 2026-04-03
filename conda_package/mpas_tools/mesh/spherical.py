import numpy as np
import xarray as xr

from mpas_tools.transects import lon_lat_to_cartesian


def recompute_angle_edge(ds_mesh):
    """
    Recompute ``angleEdge`` from edge and vertex locations on the sphere.

    Parameters
    ----------
    ds_mesh : xarray.Dataset
        An MPAS spherical mesh dataset containing edge and vertex locations.

    Returns
    -------
    angle_edge : xarray.DataArray
        ``angleEdge`` recomputed from spherical geometry.
    """
    normal_east_north = calc_edge_normal_vector(ds_mesh)
    angle_edge = xr.zeros_like(ds_mesh.angleEdge)
    angle_edge.values = np.atan2(
        normal_east_north[:, 1], normal_east_north[:, 0]
    )
    return angle_edge


def calc_edge_normal_vector(ds_mesh):
    """
    Compute edge-normal vectors projected onto local east/north coordinates.

    Parameters
    ----------
    ds_mesh : xarray.Dataset
        An MPAS spherical mesh dataset containing edge and vertex locations.

    Returns
    -------
    normal_east_north : numpy.ndarray
        A ``(nEdges, 2)`` array of unit normal vectors in local east/north
        coordinates.
    """
    edge_cartesian = np.array(
        lon_lat_to_cartesian(
            ds_mesh.lonEdge, ds_mesh.latEdge, 1.0, degrees=False
        )
    )

    vertex_1 = ds_mesh.verticesOnEdge.isel(TWO=0).values - 1
    vertex_2 = ds_mesh.verticesOnEdge.isel(TWO=1).values - 1

    lon_vertex_1 = ds_mesh.lonVertex.isel(nVertices=vertex_1)
    lat_vertex_1 = ds_mesh.latVertex.isel(nVertices=vertex_1)
    lon_vertex_2 = ds_mesh.lonVertex.isel(nVertices=vertex_2)
    lat_vertex_2 = ds_mesh.latVertex.isel(nVertices=vertex_2)

    vertex_1_cartesian = np.array(
        lon_lat_to_cartesian(lon_vertex_1, lat_vertex_1, 1.0, degrees=False)
    )
    vertex_2_cartesian = np.array(
        lon_lat_to_cartesian(lon_vertex_2, lat_vertex_2, 1.0, degrees=False)
    )

    dvertex_cartesian = vertex_2_cartesian - vertex_1_cartesian
    normal_cartesian = np.cross(dvertex_cartesian, edge_cartesian, axis=0)

    edge_east, edge_north = calc_vector_east_north(
        edge_cartesian[0, :], edge_cartesian[1, :], edge_cartesian[2, :]
    )

    normal_east_north = np.zeros((ds_mesh.sizes['nEdges'], 2))
    normal_east_north[:, 0] = np.sum(edge_east * normal_cartesian, axis=0)
    normal_east_north[:, 1] = np.sum(edge_north * normal_cartesian, axis=0)

    norm = np.linalg.norm(normal_east_north, axis=1)
    nonzero = norm > 0.0
    normal_east_north[nonzero, :] /= norm[nonzero, np.newaxis]

    return normal_east_north


def calc_vector_east_north(x, y, z):
    """
    Compute local east and north unit vectors on the sphere.

    Parameters
    ----------
    x, y, z : numpy.ndarray
        Cartesian coordinates of points on the unit sphere.

    Returns
    -------
    east, north : tuple of numpy.ndarray
        Local east and north unit vectors, each with shape ``(3, nPoints)``.
    """
    axis = np.array([0.0, 0.0, 1.0])
    xyz = np.stack((x, y, z), axis=1)
    east = np.cross(axis, np.transpose(xyz), axis=0)
    north = np.cross(np.transpose(xyz), east, axis=0)

    east /= np.linalg.norm(east, axis=0)
    north /= np.linalg.norm(north, axis=0)

    return east, north
