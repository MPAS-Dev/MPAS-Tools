import collections

point = collections.namedtuple('Point', ['x', 'y', 'z'])


def circumcenter(on_sphere, x1, y1, z1, x2, y2, z2, x3, y3, z3):  # {{{
    p1 = point(x1, y1, z1)
    p2 = point(x2, y2, z2)
    p3 = point(x3, y3, z3)
    if on_sphere:
        a = (p2.x - p3.x)**2 + (p2.y - p3.y)**2 + (p2.z - p3.z)**2
        b = (p3.x - p1.x)**2 + (p3.y - p1.y)**2 + (p3.z - p1.z)**2
        c = (p1.x - p2.x)**2 + (p1.y - p2.y)**2 + (p1.z - p2.z)**2

        pbc = a * (-a + b + c)
        apc = b * (a - b + c)
        abp = c * (a + b - c)

        xv = (pbc * p1.x + apc * p2.x + abp * p3.x) / (pbc + apc + abp)
        yv = (pbc * p1.y + apc * p2.y + abp * p3.y) / (pbc + apc + abp)
        zv = (pbc * p1.z + apc * p2.z + abp * p3.z) / (pbc + apc + abp)
    else:
        d = 2 * (p1.x * (p2.y - p3.y) + p2.x *
                 (p3.y - p1.y) + p3.x * (p1.y - p2.y))

        xv = ((p1.x**2 + p1.y**2) * (p2.y - p3.y) + (p2.x**2 + p2.y**2)
              * (p3.y - p1.y) + (p3.x**2 + p3.y**2) * (p1.y - p2.y)) / d
        yv = ((p1.x**2 + p1.y**2) * (p3.x - p2.x) + (p2.x**2 + p2.y**2)
              * (p1.x - p3.x) + (p3.x**2 + p3.y**2) * (p2.x - p1.x)) / d
        zv = 0.0

        # Optional method to use barycenter instead.
        # xv = p1.x + p2.x + p3.x
        # xv = xv / 3.0
        # yv = p1.y + p2.y + p3.y
        # yv = yv / 3.0
    return point(xv, yv, zv)

# }}}