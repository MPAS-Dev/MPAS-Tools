#!/usr/bin/env python

import numpy as np

from mpas_tools.cime.constants import constants
from mpas_tools.transects import (
    angular_distance,
    cartesian_to_great_circle_distance,
    lon_lat_to_cartesian,
    subdivide_great_circle,
    subdivide_planar,
)
from mpas_tools.vector import Vector


def test_subdivide_great_circle():
    _, _, x, y, z, earth_radius = _get_transect()
    max_res = 10e3
    x_div, y_div, z_div, d, d_div = subdivide_great_circle(
        x, y, z, max_res, earth_radius
    )
    assert np.amax(d_div <= max_res)


def test_cartesian_to_great_circle_distance():
    _, _, x, y, z, earth_radius = _get_transect()
    distance = cartesian_to_great_circle_distance(x, y, z, earth_radius)
    assert len(distance) == len(x)


def test_subdivide_planar():
    x = 1e3 * np.array([-10.0, 0.0, 10.0])
    y = 1e3 * np.array([-10.0, 0.0, 10.0])
    max_res = 1e3
    x_div, y_div, d, d_div = subdivide_planar(y, x, max_res)
    assert np.amax(d_div <= max_res)


def test_angular_distance_xyz():
    _, _, x, y, z, earth_radius = _get_transect()
    ang_dist = angular_distance(x=x, y=y, z=z)
    assert len(ang_dist) == len(x) - 1


def test_angular_distance_first_second():
    _, _, x, y, z, earth_radius = _get_transect()
    ang_dist1 = angular_distance(x=x, y=y, z=z)
    first = Vector(x[0:-1], y[0:-1], z[0:-1])
    second = Vector(x[1:], y[1:], z[1:])
    ang_dist2 = first.angular_distance(second)
    assert np.all(ang_dist1 == ang_dist2)


def test_intersects_scalar():
    a1 = _lon_lat_to_vector(-10.0, -10.0)
    a2 = _lon_lat_to_vector(10.0, 10.0)
    b1 = _lon_lat_to_vector(-10.0, 10.0)
    b2 = _lon_lat_to_vector(10.0, -10.0)
    assert Vector.intersects(a1, a2, b1, b2)


def test_intersects_array():
    a1 = _lon_lat_to_vector(np.array([-10.0, -5.0]), np.array([-10.0, -5.0]))
    a2 = _lon_lat_to_vector(np.array([-5.0, 10.0]), np.array([-5.0, 10.0]))
    b1 = _lon_lat_to_vector(-10.0, 10.0)
    b2 = _lon_lat_to_vector(10.0, -10.0)
    result = Vector.intersects(a1, a2, b1, b2)
    assert np.all(result == np.array([False, True]))


def test_intersection_scalar():
    a1 = _lon_lat_to_vector(-10.0, -10.0)
    a2 = _lon_lat_to_vector(10.0, 10.0)
    b1 = _lon_lat_to_vector(-10.0, 10.0)
    b2 = _lon_lat_to_vector(10.0, -10.0)
    earth_radius = constants['SHR_CONST_REARTH']
    cross = _lon_lat_to_vector(0.0, 0.0)
    point = Vector.intersection(a1, a2, b1, b2)
    assert np.all(
        np.isclose(
            [
                earth_radius * point.x,
                earth_radius * point.y,
                earth_radius * point.z,
            ],
            [cross.x, cross.y, cross.z],
        )
    )


def test_intersection_array():
    a1 = _lon_lat_to_vector(np.array([-10.0]), np.array([-10.0]))
    a2 = _lon_lat_to_vector(np.array([10.0]), np.array([10.0]))
    b1 = _lon_lat_to_vector(-10.0, 10.0)
    b2 = _lon_lat_to_vector(10.0, -10.0)
    earth_radius = constants['SHR_CONST_REARTH']
    cross = _lon_lat_to_vector(0.0, 0.0)
    point = Vector.intersection(a1, a2, b1, b2)
    assert np.all(
        np.isclose(
            [
                earth_radius * point.x[0],
                earth_radius * point.y[0],
                earth_radius * point.z[0],
            ],
            [cross.x, cross.y, cross.z],
        )
    )


def _get_transect():
    lon = np.array([-10.0, 0.0, 10.0])
    lat = np.array([-10.0, 0.0, 10.0])
    earth_radius = constants['SHR_CONST_REARTH']
    x, y, z = lon_lat_to_cartesian(lon, lat, earth_radius, degrees=True)
    return lon, lat, x, y, z, earth_radius


def _lon_lat_to_vector(lon, lat):
    earth_radius = constants['SHR_CONST_REARTH']
    x, y, z = lon_lat_to_cartesian(lon, lat, earth_radius, degrees=True)
    return Vector(x=x, y=y, z=z)


if __name__ == '__main__':
    test_subdivide_great_circle()
    test_cartesian_to_great_circle_distance()
    test_subdivide_planar()
    test_angular_distance_xyz()
    test_angular_distance_first_second()
    test_intersects_scalar()
    test_intersects_array()
    test_intersection_scalar()
    test_intersection_array()
