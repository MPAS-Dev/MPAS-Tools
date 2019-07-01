from __future__ import absolute_import, division, print_function

import os
import sys
import types
import numpy as np

from limited_area.points import PointsParser
from limited_area.mesh import latlon_to_xyz, xyz_to_latlon
from limited_area.mesh import rotate_about_vector
import numpy as np

EARTH_RADIUS = 6371229 # Meters

""" region_spec.py - Provide a number of operations for defining a 
region. """

if sys.version_info[0] > 2:
    create_bound_method = types.MethodType
else:
    def create_bound_method(func, obj):
        return types.MethodType(func, obj, obj.__class__)


def normalize_cords(lat, lon):
    """ Returned lat, and lon to be in radians and the same 
    range as MPAS - Lat: -pi/2 to pi/2 - Lon: 0 to 2*pi
       
    Lat - Latitude in degrees
    Lon - Longitude in degrees

    """
    lat *= np.pi / 180.0
    lon *= np.pi / 180.0

    return lat, lon


class RegionSpec:
    """ RegionSpec - Method for generating regional specifications
    Region spec works upon a contract. It will need the following 
    information when it is called:

    filename  - Path to the file that is to be read 

    And will then return, the following:

    filename         - Output filename (if desired and specified in the specification file)
    points array     - A 1-dimensional list of lat lon cords specifying boundary 
                       in counter clockwise order
    in-point         - pair of points that inside the boundary
    algorithm choice - The desired algorithm for choosing boundary points and 
                       relaxation layers (if specified within the specification file)
    """
    # TODO: Update __init__ with fileName
    def __init__(self, *args, **kwargs):
        """ init for region Spec

        Keyword Arguments: 
            DEBUG - Debug value for verbose output - Default 0
        """
        # Keyword Args
        self._DEBUG_ = kwargs.get('DEBUG', 0)
        self._gen_spec = create_bound_method(PointsParser, self)

    def gen_spec(self, fileName, *args, **kwargs):
        """ Generate the specifications and return, name, in point and a list of points.

        Call the method we bound above, and then do any processing here
        to do things like convert coordinates to radians, or anything else
        we need to get return contract variables.
        
        fileName - The file that specifies the region

        Return values:
            name     - The name of the region
            in_point - A point that is within the region
            points   - A 1-Dimensional list of latitude and longitude points
                       (in degrees) that list the boundary points of the region
                       in counter-clockwise. ie: [lat1, lon1, lat2, lon2, ... , latN, lonN]
        """

        self._gen_spec(fileName, *args, **kwargs)

        if self.type == 'custom':
            if self._DEBUG_ > 0:
                print("DEBUG: Using a custom polygon for generating a region")

            self.points = np.array(self.points)

            # Convert the points to radians and set them to be between
            # Lon: 0 to 2*Pi and Lat: -pi to +pi
            for cord in range(0, len(self.points), 2):
                 self.points[cord], self.points[cord+1] = normalize_cords(
                                                          self.points[cord], 
                                                          self.points[cord+1])

            self.in_point[0], self.in_point[1] = normalize_cords(
                                                    self.in_point[0],
                                                    self.in_point[1])

            return self.name, self.in_point, [self.points]
        elif self.type == 'circle':
            if self._DEBUG_ > 0:
                print("DEBUG: Using the circle method for region generation")

            self.in_point[0], self.in_point[1] = normalize_cords(
                                                    self.in_point[0],
                                                    self.in_point[1])

            # Convert to meters, then divide by radius to get radius upon sphere w/ r = 1
            self.radius = (self.radius * 1000) / EARTH_RADIUS
            self.points = self.circle(self.in_point[0], self.in_point[1], self.radius)
            return self.name, self.in_point, [self.points.flatten()]
        elif self.type == 'ellipse':
            if self._DEBUG_ > 0:
                print("DEBUG: Using the ellipse method for region generation")

            # Convert ellipse center point from degrees to radians
            self.in_point[0], self.in_point[1] = normalize_cords(
                                                    self.in_point[0],
                                                    self.in_point[1])

            self.points = self.ellipse(self.in_point[0], self.in_point[1],
                                       self.semimajor, self.semiminor,
                                       self.orientation)
            return self.name, self.in_point, [self.points.flatten()]
        elif self.type == 'channel':
            if self._DEBUG_ > 0:
                print("DEBUG: Using the channel method for region generation")

            if self.ulat == self.llat:
                print("ERROR: Upper and lower latitude for channel specification")
                print("ERROR: cannot be equal")
                print("ERROR: Upper-lat: ", self.ulat)
                print("ERROR: Lower-lat: ", self.llat)
                sys.exit(-1)

            self.ulat, self.llat = normalize_cords(self.ulat, self.llat)

            self.boundaries = []

            upperBdy = np.empty([100, 2])
            lowerBdy = np.empty([100, 2])

            upperBdy[:,0] = self.ulat
            upperBdy[:,1] = np.linspace(0.0, 2.0 * np.pi, 100)

            lowerBdy[:,0] = self.llat
            lowerBdy[:,1] = np.linspace(0.0, 2.0 * np.pi, 100)

            self.in_point = np.array([(self.ulat + self.llat) / 2, 0])

            self.boundaries.append(upperBdy.flatten())
            self.boundaries.append(lowerBdy.flatten())

            return self.name, self.in_point, self.boundaries

    def circle(self, center_lat, center_lon, radius):
        """ Return a list of latitude and longitude points in degrees that
        area radius away from (center_lat, center_lon)

        center_lat - Circle center latitude in radians
        center_lon - Circle center longitude in radians
        radius     - Radius of desire circle in radians upon the unit sphere

        """

        P = []

        if self._DEBUG_ > 1:
            print("DEBUG: center_lat: ", center_lat,
                        " center_lon: ", center_lon,
                        " radius: ", radius)


        C = latlon_to_xyz(center_lat, center_lon, 1.0)

        # Find a point not equal to C or -C
        K = np.zeros(3)
        K[0] = 1
        K[1] = 0
        K[2] = 0

        if abs(np.dot(K,C)) >= 0.9:
            K[0] = 0
            K[1] = 1
            K[2] = 0


        # S is then a vector orthogonal to C
        S = np.cross(C, K)
        S = S / np.linalg.norm(S)

        P0 = rotate_about_vector(C, S, radius)

        for r in np.linspace(0.0, 2.0*np.pi, 100):
            P.append(rotate_about_vector(P0, C, r))

        ll = []
        # TODO: The efficiency here can be improved for memory
        # and probably comp time
        for i in range(len(P)):
            ll.append(xyz_to_latlon(P[i])) # Convert back to latlon

        return np.array(ll)

    def ellipse(self, center_lat, center_lon, semi_major, semi_minor, orientation):
        """ Return a list of points that form an ellipse around [center_lat, center_lon] 

        center_lat - The center latitude of the ellipse in degrees
        center_lon - The center longitude of the ellipse in degrees
        semi_major - The length of the semi_major axies in meters
        semi_minor - The legnth of the semi_minor axies in meters
        orientaiton - The orientation of the desired ellipse, rotated clockwise from north.

        Given the center_lat and center_lon of an ellipse, create a region that forms an ellipse
        with the semi major axies length being == `semi_major` and the semi minor axies being ==
        `semi_minor` and rotate the ellipse clockwise by `orientation` from due North.

        """

        P = []

        # Convert ellipse center from (lat,lon) to Cartesian
        C = latlon_to_xyz(center_lat, center_lon, 1.0)

        # Convert semi-major and semi-minor axis lengths from meters to radians
        semi_major = semi_major / EARTH_RADIUS
        semi_minor = semi_minor / EARTH_RADIUS

        # Convert orientation angle to radians
        orientation = orientation * np.pi / 180.0

        # Find a point not equal to C or -C
        K = np.zeros(3)
        K[0] = 0
        K[1] = 0
        K[2] = 1

        if abs(np.dot(K,C)) >= 0.9:
            K[0] = 1
            K[1] = 0
            K[2] = 0

        # S is then a vector orthogonal to C
        S = np.cross(C, K)
        S = S / np.linalg.norm(S)

        for r in np.linspace(0.0, 2.0*np.pi, 100):
            radius = np.sqrt((semi_major * np.cos(r))**2 + (semi_minor * np.sin(r))**2)
            P0 = rotate_about_vector(C, S, radius)
            ang = np.arctan2(semi_minor * np.sin(r), semi_major * np.cos(r))
            P.append(rotate_about_vector(P0, C, ang-orientation))

        ll = []
        # TODO: The efficency here can be improved for memory
        # and probably comp time
        for i in range(len(P)):
            ll.append(xyz_to_latlon(P[i])) # Convert back to latlon

        return np.array(ll)
