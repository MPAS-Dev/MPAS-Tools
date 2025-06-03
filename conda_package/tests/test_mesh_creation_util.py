import collections

import numpy as np

from mpas_tools.mesh.creation.util import circumcenter

point = collections.namedtuple('Point', ['x', 'y', 'z'])


def test_circumcenter_planar():
    # List of test cases: (cell0, cell1, cell2, expected)
    test_cases = [
        # From a regular hex mesh
        (
            (0.0, 17320.50807569, 0.0),
            (10000.0, 17320.50807569, 0.0),
            (5000.0, 25980.76211353, 0.0),
            (5000.0, 20207.25942164, 0.0),
        ),
        # From an Antarctic mesh
        (
            (2194652.05105505, 1144816.29348609, 0.0),
            (2198726.62725077, 1145337.65538494, 0.0),
            (2195474.22967839, 1148983.95766333, 0.0),
            (2196492.12915069, 1146618.22089598, 0.0),
        ),
        # From an Antarctic mesh
        (
            (2195474.22967839, 1148983.95766333, 0.0),
            (2198726.62725077, 1145337.65538494, 0.0),
            (2199279.33328257, 1149240.68069275, 0.0),
            (2197485.28573298, 1147504.08822421, 0.0),
        ),
        # From an Antarctic mesh
        (
            (2150629.6663189, 1441451.4899963, 0.0),
            (2149349.39244275, 1445285.59533973, 0.0),
            (2146912.55160507, 1442187.732174, 0.0),
            (2149013.33940856, 1443042.57600619, 0.0),
        ),
        # From an Antarctic mesh
        (
            (2570868.91779026, -97880.09977304, 0.0),
            (2575332.71614832, -97980.66509542, 0.0),
            (2573585.71564672, -95138.63592864, 0.0),
            (2573113.05578246, -97387.13758518, 0.0),
        ),
    ]

    for i, (cell0, cell1, cell2, expected) in enumerate(test_cases):
        x1, y1, z1 = cell0
        x2, y2, z2 = cell1
        x3, y3, z3 = cell2

        xv, yv, zv = circumcenter(False, x1, y1, z1, x2, y2, z2, x3, y3, z3)
        assert np.allclose([xv, yv, zv], expected, atol=1e-8), (
            f'Planar circumcenter failed for case {i}: got {[xv, yv, zv]}, '
            f'expected {expected}'
        )

        xv_old, yv_old, zv_old = _old_circumcenter(
            False, x1, y1, z1, x2, y2, z2, x3, y3, z3
        )
        print(f'Planar circumcenter (case {i}):')
        print('  circumcenter:      ', (xv, yv, zv))
        print('  _old_circumcenter: ', (xv_old, yv_old, zv_old))
        print('  difference:        ', (xv - xv_old, yv - yv_old, zv - zv_old))
        assert np.array_equal([xv, yv, zv], [xv_old, yv_old, zv_old]), (
            f'circumcenter and _old_circumcenter results differ in case {i}: '
            f'({xv}, {yv}, {zv}) != ({xv_old}, {yv_old}, {zv_old})'
        )


def test_circumcenter_sphere():
    # Use provided cell coordinates and vertex as expected output
    x1, y1, z1 = -3590.95704026, -684.52527287, -6371218.95125671
    x2, y2, z2 = -9926.94769501, 18734.26792821, -6371184.72274307
    x3, y3, z3 = 12921.92507847, 11340.69768405, -6371196.8028643
    expected = (296.83807146, 11327.05862805, -6371209.92418473)

    xv, yv, zv = circumcenter(True, x1, y1, z1, x2, y2, z2, x3, y3, z3)
    assert np.allclose([xv, yv, zv], expected, atol=1e-8), (
        f'Spherical circumcenter failed: got {[xv, yv, zv]}, expected '
        f'{expected}'
    )

    xv_old, yv_old, zv_old = _old_circumcenter(
        True, x1, y1, z1, x2, y2, z2, x3, y3, z3
    )
    # Print both sets and their difference
    print('Spherical circumcenter:')
    print('  circumcenter:      ', (xv, yv, zv))
    print('  _old_circumcenter: ', (xv_old, yv_old, zv_old))
    print('  difference:        ', (xv - xv_old, yv - yv_old, zv - zv_old))
    # Check for exact equality
    assert np.array_equal([xv, yv, zv], [xv_old, yv_old, zv_old]), (
        f'circumcenter and _old_circumcenter results differ: '
        f'({xv}, {yv}, {zv}) != ({xv_old}, {yv_old}, {zv_old})'
    )


def _old_circumcenter(on_sphere, x1, y1, z1, x2, y2, z2, x3, y3, z3):
    """
    Compute the circumcenter of the triangle (possibly on a sphere)
    with the three given vertices in Cartesian coordinates.

    Returns
    -------
    center : point
        The circumcenter of the triangle with x, y and z attributes

    """
    p1 = point(x1, y1, z1)
    p2 = point(x2, y2, z2)
    p3 = point(x3, y3, z3)
    if on_sphere:
        a = (p2.x - p3.x) ** 2 + (p2.y - p3.y) ** 2 + (p2.z - p3.z) ** 2
        b = (p3.x - p1.x) ** 2 + (p3.y - p1.y) ** 2 + (p3.z - p1.z) ** 2
        c = (p1.x - p2.x) ** 2 + (p1.y - p2.y) ** 2 + (p1.z - p2.z) ** 2

        pbc = a * (-a + b + c)
        apc = b * (a - b + c)
        abp = c * (a + b - c)

        xv = (pbc * p1.x + apc * p2.x + abp * p3.x) / (pbc + apc + abp)
        yv = (pbc * p1.y + apc * p2.y + abp * p3.y) / (pbc + apc + abp)
        zv = (pbc * p1.z + apc * p2.z + abp * p3.z) / (pbc + apc + abp)
    else:
        d = 2 * (
            p1.x * (p2.y - p3.y) + p2.x * (p3.y - p1.y) + p3.x * (p1.y - p2.y)
        )

        xv = (
            (p1.x**2 + p1.y**2) * (p2.y - p3.y)
            + (p2.x**2 + p2.y**2) * (p3.y - p1.y)
            + (p3.x**2 + p3.y**2) * (p1.y - p2.y)
        ) / d
        yv = (
            (p1.x**2 + p1.y**2) * (p3.x - p2.x)
            + (p2.x**2 + p2.y**2) * (p1.x - p3.x)
            + (p3.x**2 + p3.y**2) * (p2.x - p1.x)
        ) / d
        zv = 0.0

        # Optional method to use barycenter instead.
        # xv = p1.x + p2.x + p3.x
        # xv = xv / 3.0
        # yv = p1.y + p2.y + p3.y
        # yv = yv / 3.0
    return point(xv, yv, zv)
