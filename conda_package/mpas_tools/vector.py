import numpy as np


class Vector:
    """
    A class for representing Cartesian vectors with ``x``, ``y`` and ``z``
    components that are either ``float`` or ``numpy.array`` objects of
    identical size.

    Attributes
    ----------
    x : float or numpy.ndarray
        The x component(s)

    y : float or numpy.ndarray
        The y component(s)

    z : float or numpy.ndarray
        The z component(s)
    """
    def __init__(self, x, y, z):
        """
        A class for representing Cartesian vectors with ``x``, ``y`` and ``z``
        components that are either ``float`` or ``numpy.array`` objects of
        identical size.

        Parameters
        ----------
        x : float or numpy.ndarray
            The x component(s)

        y : float or numpy.ndarray
            The y component(s)

        z : float or numpy.ndarray
            The z component(s)
        """
        self.x = x
        self.y = y
        self.z = z

    def angular_distance(self, other):
        """
        Compute angular distance between points on the sphere, following:
        https://en.wikipedia.org/wiki/Great-circle_distance

        Parameters
        ----------
        other : mpas_tools.vector.Vector
            The vector to compute the angular distance to

        Returns
        -------
        angularDistance : numpy.ndarray
            The angular distance (in radians) between segments of the transect.
        """
        angular_distance = np.arctan2(self.cross(other).mag(), self.dot(other))
        return angular_distance

    @staticmethod
    def intersects(a1, a2, b1, b2):
        """
        Based on https://stackoverflow.com/a/26669130/7728169
        Determine if the great circle arc from ``a1`` to ``a2`` intersects that
        from ``b1`` to ``b2``.

        Parameters
        ----------
        a1 : mpas_tools.vector.Vector
            Cartesian coordinates of the end point of a great circle arc.
            The types of the attributes ``x``, ``y``, and ``z`` must either be
            ``numpy.arrays`` of identical size for all 4 vectors (in which case
            intersections are found element-wise), or scalars for
            at least one of either ``a1`` and ``a2`` or ``b1`` and ``b2``.

        a2 : mpas_tools.vector.Vector
            Cartesian coordinates of the other end point of a great circle arc.

        b1 : mpas_tools.vector.Vector
            Cartesian coordinates of an end point of a second great circle arc.

        b2 : mpas_tools.vector.Vector
            Cartesian coordinates of the other end point of the second great
            circle arc.

        Returns
        -------
        intersect : numpy.ndarray
            A boolean array of the same size as ``a1`` and ``a2`` or ``b1`` and
            ``b2``, whichever is greater, indicating if the particular pair of
            arcs intersects
        """
        return np.logical_and(Vector.straddles(a1, a2, b1, b2),
                              Vector.straddles(b1, b2, a1, a2))

    @staticmethod
    def intersection(a1, a2, b1, b2):
        """
        Based on https://stackoverflow.com/a/26669130/7728169
        Find the intersection point as a unit vector between great circle arc
        from ``a1`` to ``a2`` and from ``b1`` to ``b2``.  The arcs should have
        already have been found to intersect by calling ``intersects()``

        Parameters
        ----------
        a1 : mpas_tools.vector.Vector
            Cartesian coordinates of the end point of a great circle arc.
            The types of the attributes ``x``, ``y``, and ``z`` must either be
            ``numpy.arrays`` of identical size for all 4 vectors (in which case
            intersections are found element-wise), or scalars for
            at least one of either ``a1`` and ``a2`` or ``b1`` and ``b2``.

        a2 : mpas_tools.vector.Vector
            Cartesian coordinates of the other end point of a great circle arc.

        b1 : mpas_tools.vector.Vector
            Cartesian coordinates of an end point of a second great circle arc.

        b2 : mpas_tools.vector.Vector
            Cartesian coordinates of the other end point of the second great
            circle arc.

        Returns
        -------
        points : mpas_tools.vector.Vector
            An array of Cartesian points *on the unit sphere* indicating where
            the arcs intersect
        """
        points = (a1.cross(a2)).cross(b1.cross(b2))
        s = np.sign(Vector.det(a1, b1, b2))/points.mag()
        points = Vector(s*points.x,  s*points.y, s*points.z)
        return points

    @staticmethod
    def straddles(a1, a2, b1, b2):
        """
        Based on https://stackoverflow.com/a/26669130/7728169
        Determines if the great circle segment determined by (a1, a2)
        straddles the great circle determined by (b1, b2)

        Parameters
        ----------
        a1: mpas_tools.vector.Vector
            Cartesian coordinates of first end point of first great circle arc.
            The types of the attributes ``x``, ``y``, and ``z`` must either be
            ``numpy.arrays`` of identical size for all 4 vectors (in which case
            intersections are found element-wise), or scalars for
            at least one of either the ``a``s or the ``b``s.

        a2 : mpas_tools.vector.Vector
            Second end point of first great circle arc.

        b1 : mpas_tools.vector.Vector
            First end point of second great circle arc.

        b2 : mpas_tools.vector.Vector
            Second end point of second great circle arc.

        Returns
        -------
        straddle : numpy.ndarray
            A boolean array of the same size as the ``a``s or the ``b``s, whichever
            is greater, indicating if the great circle segment determined by
            (a1, a2) straddles the great circle determined by (b1, b2)
        """
        return Vector.det(a1, b1, b2) * Vector.det(a2, b1, b2) < 0

    def dot(self, other):
        """
        Compute the dot product between this vector and ``other``.

        Parameters
        ----------
        other : mpas_tools.vector.Vector
            The other vector

        Returns
        -------
        dot_product : numpy.ndarray
            The dot product
        """
        return self.x * other.x + self.y * other.y + self.z * other.z

    def cross(self, other):
        """
        Compute the dot product between this vector and ``other``.

        Parameters
        ----------
        other : mpas_tools.vector.Vector
            The other vector

        Returns
        -------
        cross_product : mpas_tools.vector.Vector
            The cross product
        """
        return Vector(self.y * other.z - self.z * other.y,
                      self.z * other.x - self.x * other.z,
                      self.x * other.y - self.y * other.x)

    @staticmethod
    def det(v1, v2, v3):
        """
        The determinant of the matrix defined by the three ``Vector`` objects

        Parameters
        ----------
        v1 : mpas_tools.vector.Vector
            First row of the matrix

        v2 : mpas_tools.vector.Vector
            Second row

        v3 : mpas_tools.vector.Vector
            Third row

        Returns
        -------
        determinant : numpy.ndarray
            The determinant of the matrix
        """
        return v1.dot(v2.cross(v3))

    def mag(self):
        """
        The magnitude of the vector

        Returns
        -------
        magnitude : numpy.ndarray
            The magnitude of the vector
        """
        return np.sqrt(self.dot(self))
