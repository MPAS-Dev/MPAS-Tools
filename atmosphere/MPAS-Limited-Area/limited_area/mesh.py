from __future__ import absolute_import, division, print_function
import sys
import os

import numpy as np
from netCDF4 import Dataset

""" mesh.py - Handle NetCDF file operations as well as helpful
calculations upon on MPAS grid. """

class MeshHandler:
    """ Handle the operations related to NetCDF/MPAS grids. """

    def __init__(self, fname, mode, format='NETCDF3_64BIT_OFFSET', *args, **kwargs):
        """ Open fname with mode, for either reading, or creating
        
        fname - A valid netCDF4 file OR, if `mode==w` then a name of a 
                desired netCDF that will be created for writing.
        mode  - Mode for opening the file, options are 'r' and 'w' for read
                and write respectively.
        format - NetCDF file format for the regional mesh. This is only used
                if `mode==w`. For more information on NetCDF versions see the
                netCDF4 Python documentation found here:

                http://unidata.github.io/netcdf4-python/netCDF4/index.html#netCDF4.Dataset.__init__
        """
        self._DEBUG_ = kwargs.get('DEBUG', 0)
        self.fname = fname

        if mode == 'r':
            if self.check_file(fname):
                self._load_vars()
                return
            else:
                sys.exit(-1)
        elif mode == 'w':
            self.create_file(fname, mode, format)


    def create_file(self, fname, mode, format):
        """ Create and open a new NetCDF file with name fname and access mode mode
        with the NetCDF file version being format. """
        try:
            self.mesh = Dataset(fname, mode, format=format)
            return
        except:
            print("ERROR: There was a problem creating the file ", fname)
            sys.exit(-1)


    def check_file(self, fname):
        """ Check to see that fname exists and it is a valid NetCDF file """
        if os.path.isfile(fname):
            try:
                mesh = open(fname, 'rb')
                nc_bytes = mesh.read()
                mesh.close()
                
                self.mesh = Dataset(fname, 'r', memory=nc_bytes)
                return True
            except OSError as E: 
                print("ERROR: ", E)
                print("ERROR: This file was not a valid NetCDF file")
                sys.exit(-1)
        else:
            print("ERROR: This file did not exist!")
            return False

    def _load_vars(self):
        """ Pre-load variables to avoid multiple, unnecessary IO calls 
            
            Pulling variables from a netCDF4 interface like the following:
            ```self.mesh.variables['lonCell']{[:]```, will read it from disk
            each time, thus we can pre-load the variables into memory to reduce
            I/O calls.
        """
        if self._DEBUG_ > 2:
            print("DEBUG: In Load Vars")

        # Dimensions
        self.nCells = self.mesh.dimensions['nCells'].size
        self.nEdges = self.mesh.dimensions['nEdges'].size
        self.maxEdges = self.mesh.dimensions['maxEdges'].size
        self.nVertices = self.mesh.dimensions['nVertices'].size
        self.vertexDegree = self.mesh.dimensions['vertexDegree'].size

        # Variables
        self.latCells = self.mesh.variables['latCell'][:]
        self.lonCells = self.mesh.variables['lonCell'][:]

        self.nEdgesOnCell = self.mesh.variables['nEdgesOnCell'][:]
        self.cellsOnCell = self.mesh.variables['cellsOnCell'][:]
        self.cellsOnEdge = self.mesh.variables['cellsOnEdge'][:]
        self.cellsOnVertex = self.mesh.variables['cellsOnVertex'][:]

        self.indexToCellIDs = self.mesh.variables['indexToCellID'][:]
        self.indexToEdgeIDs = self.mesh.variables['indexToEdgeID'][:]
        self.indexToVertexIDs = self.mesh.variables['indexToVertexID'][:]

        # Attributes
        self.sphere_radius = self.mesh.sphere_radius

        self.variables = { 'latCells' : self.latCells,
                           'lonCells' : self.lonCells,
                           'nEdgesOnCell' : self.nEdgesOnCell,
                           'cellsOnCell' : self.cellsOnCell,
                           'cellsOnEdge' : self.cellsOnEdge,
                           'cellsOnVertex' : self.cellsOnVertex,
                           'indexToCellID' : self.indexToCellIDs,
                           'indexToEdgeID' : self.indexToEdgeIDs,
                           'indexToVertexID' : self.indexToVertexIDs
                         }


    def nearest_cell(self, lat, lon):
        """ Find the nearest cell of this mesh to lat and lon

        lat - Latitude - Radians
        lon - Longitude - Radians
        """
        nearest_cell = 0 # Start at the first cell
        current_cell = -1

        while (nearest_cell != current_cell):
            current_cell = nearest_cell
            current_distance = sphere_distance(self.latCells[current_cell],
                                               self.lonCells[current_cell],
                                               lat,
                                               lon,
                                               self.sphere_radius)
            
            nearest_cell = current_cell
            nearest_distance = current_distance
            
            for edges in range(self.nEdgesOnCell[current_cell]):
                iCell = self.cellsOnCell[current_cell, edges] - 1
                if (iCell <= self.nCells):
                    iDistance = sphere_distance(self.latCells[iCell],
                                                self.lonCells[iCell],
                                                lat,
                                                lon,
                                                self.sphere_radius)

                    if (iDistance <= nearest_distance):
                        nearest_cell = iCell
                        nearest_distance = iDistance

    
        if self._DEBUG_ > 3:
            print("DEBUG: nearest_cell latLon: ", nearest_cell, '\t',
                                                  self.latCells[nearest_cell] * (180.0/np.pi),
                                                  self.lonCells[nearest_cell] * (180.0/np.pi),
                  ' Given lat lon: ', lat * (180.0/np.pi), lon * (180.0/np.pi))


        return nearest_cell

    def create_graph_file(self, graphFname):
        """ Create a graph.info file for the current mesh """

        # In the limited_area program, this function will always create
        # a graph.info file for a regional mesh. Thus, the variables here
        # are not preloaded and are being read from disk
        nCells = self.mesh.dimensions['nCells'].size
        nEdges = self.mesh.dimensions['nEdges'].size

        nEdgesOnCell = self.mesh.variables['nEdgesOnCell'][:]
        cellsOnCell = self.mesh.variables['cellsOnCell'][:]
        cellsOnEdge = self.mesh.variables['cellsOnEdge'][:]

        lines = []
        line = ''

        nEdgesInterior = 0

        for i in range(nEdges):
            if (cellsOnEdge[i,0] > 0 and cellsOnEdge[i,1] > 0):
                nEdgesInterior = nEdgesInterior + 1

        with open(graphFname, 'w') as f:
            f.write(repr(nCells)+' '+repr(nEdgesInterior)+'\n')
            for i in range(nCells):
                for j in range(nEdgesOnCell[i]):
                    if (cellsOnCell[i,j] > 0):
                        f.write(repr(cellsOnCell[i,j])+' ')
                f.write('\n')

        return graphFname

    def subset_fields(self, 
                      regionalFname, 
                      bdyMaskCell,
                      bdyMaskEdge,
                      bdyMaskVertex,
                      inside,
                      unmarked,
                      format='NETCDF3_64BIT_OFFSET',
                      *args, 
                      **kwargs):
        """ Subset the current mesh and return a new regional mesh with
        subsetted fields 
        
        regionalFname -- Desired filename for the regional subset
        bdyMaskCell   -- Global mesh mask denoting regional cells
        bdyMaskEdge   -- Global mesh mask denoting regional edges
        bdyMaskVertex -- Global mesh mask denoting regional vertices
        inside        -- The integer value that was used to mark the 
                         cells, edges, vertices as being 'inside' the 
                         regional within the bdyMasks
        unmarked      -- The integer value that was used to mark cells,
                         edges, vertices as being 'outside' of the regional
                         mesh.
        """

        # Don't pass on DEBUG to the regional mess - tone down output
        kwargs.pop('DEBUG')

        indexingFields = {}
        indexingFields['indexToCellID'] = bdyMaskCell
        indexingFields['indexToEdgeID'] = bdyMaskEdge
        indexingFields['indexToVertexID'] = bdyMaskVertex
        indexingFields['cellsOnEdge'] = bdyMaskCell
        indexingFields['edgesOnCell'] = bdyMaskEdge
        indexingFields['edgesOnEdge'] = bdyMaskEdge
        indexingFields['cellsOnCell'] = bdyMaskCell
        indexingFields['verticesOnCell'] = bdyMaskVertex
        indexingFields['verticesOnEdge'] = bdyMaskVertex
        indexingFields['edgesOnVertex'] = bdyMaskEdge
        indexingFields['cellsOnVertex'] = bdyMaskCell

        glbBdyCellIDs = self.indexToCellIDs[np.where(bdyMaskCell != unmarked)] - 1
        glbBdyEdgeIDs = self.indexToEdgeIDs[np.where(bdyMaskEdge != unmarked)] - 1
        glbBdyVertexIDs = self.indexToVertexIDs[np.where(bdyMaskVertex != unmarked)] - 1


        if self._DEBUG_ > 0:
            print("DEBUG: nCells of new region: ", len(glbBdyCellIDs))
            print("DEBUG: nEdges of new region: ", len(glbBdyEdgeIDs))
            print("DEBUG: nVertex of new region: ", len(glbBdyVertexIDs))

        # Check to see the user didn't mess specifying the region. If 
        # len(bdyIndexToCellIDs) == nCells, then the specification was probably not
        # specified correctly
        force = False
        if len(glbBdyCellIDs) == self.nCells and not force:
            print("ERROR: The number of Cells in the specified region ",
                  "(", len(glbBdyCellIDs), ")")
            print("ERROR: appears to be equal number of cells in the global mesh",
                  "(", self.nCells, ")")
            print("ERROR: which means there was perhaps a problem in specifying the")
            print("ERROR: region. Please insure your region specification is correct")
            sys.exit(-1)

        # Create a new grid
        region = MeshHandler(regionalFname, 'w', format=format, *args, **kwargs)

        # Dimensions - Create dimensions
        for dim in self.mesh.dimensions:
            if dim == 'nCells':
                region.mesh.createDimension(dim, 
                                            len(glbBdyCellIDs))
            elif dim == 'nEdges':
                region.mesh.createDimension(dim, 
                                            len(glbBdyEdgeIDs))
            elif dim == 'nVertices':
                region.mesh.createDimension(dim, 
                                            len(glbBdyVertexIDs))
            else:
                if self.mesh.dimensions[dim].isunlimited():
                    region.mesh.createDimension(dim,
                                                None)
                else:
                    region.mesh.createDimension(dim,
                                                self.mesh.dimensions[dim].size)

        # Make boundary Mask's between 0 and the number of specified relaxation
        # layers
        region.mesh.createVariable('bdyMaskCell', 'i4', ('nCells',))
        region.mesh.createVariable('bdyMaskEdge', 'i4', ('nEdges',))
        region.mesh.createVariable('bdyMaskVertex', 'i4', ('nVertices')) 

        region.mesh.variables['bdyMaskCell'][:] = bdyMaskCell[bdyMaskCell != 0] - 1 
        region.mesh.variables['bdyMaskEdge'][:] = bdyMaskEdge[bdyMaskEdge != 0] - 1
        region.mesh.variables['bdyMaskVertex'][:] = bdyMaskVertex[bdyMaskVertex != 0] - 1

        scan(bdyMaskCell)
        scan(bdyMaskEdge)
        scan(bdyMaskVertex)

        # Variables - Create Variables
        for var in self.mesh.variables:
            # If we're subsetting a static file, don't copy variables for bdyMaskCell,
            # bdyMaskEdge or bdyMaskVertex if they exist in the static file as they have
            # been created by this program
            if var != 'bdyMaskCell' and var != 'bdyMaskEdge' and var != 'bdyMaskVertex':
                region.mesh.createVariable(var, self.mesh.variables[var].dtype,
                                                self.mesh.variables[var].dimensions)
                try:
                    region.mesh.variables[var].units = self.mesh.variables[var].units
                    region.mesh.variables[var].long_name = self.mesh.variables[var].long_name
                except:
                    pass

        # Subset global variables into the regional mesh and write them
        # to the regional mesh - re-indexing if necessary
        for var in self.mesh.variables:
            # If subsetting a static file, don't copy any bdyMask variables, as we
            # created them above
            if var == 'bdyMaskCell' or var == 'bdyMaskEdge' or var == 'bdyMaskVertex':
                continue

            print("Copying variable ", var, "...", end=' ', sep=''); sys.stdout.flush()
            if var in self.variables:
                arrTemp = self.variables[var] # Use the pre-loaded variable if possible
            else:
                arrTemp = self.mesh.variables[var][:] # Else, read it from disk

            if 'nCells' in self.mesh.variables[var].dimensions:
                if var in indexingFields:
                    region.mesh.variables[var][:] = reindex_field(arrTemp[glbBdyCellIDs], 
                                                                  indexingFields[var])
                    print('Done!')
                else:
                    print('')
                    region.mesh.variables[var][:] = arrTemp[glbBdyCellIDs]
            elif 'nEdges' in self.mesh.variables[var].dimensions:
                if var in indexingFields:
                    region.mesh.variables[var][:] = reindex_field(arrTemp[glbBdyEdgeIDs], 
                                                                  indexingFields[var])
                    print('Done!')
                else:
                    print('')
                    region.mesh.variables[var][:] = arrTemp[glbBdyEdgeIDs]
            elif 'nVertices' in self.mesh.variables[var].dimensions:
                if var in indexingFields:
                    region.mesh.variables[var][:] = reindex_field(arrTemp[glbBdyVertexIDs], 
                                                                  indexingFields[var])
                    print('Done!')
                else:
                    print('')
                    region.mesh.variables[var][:] = arrTemp[glbBdyVertexIDs]
            else:
                print('')
                region.mesh.variables[var][:] = arrTemp


        return region

    def copy_global_attributes(self, region):
        """ Copy the global attributes into the regional mesh, but not 'np' """
        region.mesh.on_a_sphere = self.mesh.on_a_sphere
        region.mesh.sphere_radius = self.mesh.sphere_radius


def scan(arr):
    """ For values within bdyMaskCell/Vertex/Edge etc. assign them
    values of their edges, cells etc."""
    arr[arr > 0] = np.arange(1, len(arr[arr > 0])+1)


def reindex_field(field, mmap):
    """ Re-index fields to be in range of their dimensions """
    print('reindexing field ...', end=' '); sys.stdout.flush()
    return mmap[field[:]-1]



def latlon_to_xyz(lat, lon, radius):
    """ Calculate and return x, y, z coordinations of lat, lon on the sphere that has
    radius, radius.
    lat - Latitude
    lon - Longitude
    radius - Radius of sphere
    """
    z = radius * np.sin(lat)
    x = radius * np.cos(lon) * np.cos(lat)
    y = radius * np.sin(lon) * np.cos(lat)

    return np.array([x, y, z])


def xyz_to_latlon(point):
    """ Convert a Cartesian coordinate point into a latitude,
        longitude point in radians """
    x = point[0]
    y = point[1]
    z = point[2]

    eps = float(1.0e-10)
    lat = np.arcsin(z)

    if np.fabs(x) > eps:
        if np.fabs(y) > eps:
            lon = np.arctan(np.fabs(y/x))

            if x <= 0.0 and y >= 0.0:
                lon = np.pi - lon
            elif x <= 0.0 and y <= 0.0:
                lon = lon + np.pi
            elif x >= 0.0 and y <= 0.0:
                lon = 2.0 * np.pi - lon
        else:
            if x > 0.0:
                lon = 0.0
            else:
                lon = np.pi

    elif np.fabs(y) > eps:
        if y > 0.0:
            lon = 0.5 * np.pi
        else:
            lon = 1.5 * np.pi
    else:
        lon = 0.0

    return lat, lon


def sphere_distance(lat1, lon1, lat2, lon2, radius, **kwargs):
    """ Calculate the sphere distance between point1 and point2. 

    lat1 - Float - Radians - -pi:pi
    lon1 - Float - Radians - 0:2*pi
    lat2 - Float - Radians - -pi:pi
    lon2 - Float - Radians - 0:2*pi
    radius - Radius of the earth (or sphere) - Units can be ignored

    """ 
    return (2 * radius * np.arcsin(
                         np.sqrt(
                         np.sin(0.5 * (lat2 - lat1))**2
                       + np.cos(lat1) 
                       * np.cos(lat2) 
                       * np.sin(0.5 * (lon2 - lon1))**2)))


def rotate_about_vector(X, U, theta):
   """  Rotates the point X through an angle theta about the vector U

   X - [x, y, z] - The point to be rotated
   U - [u, v, w] - The point to rotate X around
   theta - The angle to rotate X around U

   Reference: https://sites.google.com/site/glennmurray/Home/rotation-matrices-and-formulas/rotation-about-an-arbitrary-axis-in-3-dimensions
   """

   x = X[0]
   y = X[1]
   z = X[2]

   u = U[0]
   v = U[1]
   w = U[2]

   vw2 = v*v + w*w
   uw2 = u*u + w*w
   uv2 = u*u + v*v
   m = np.sqrt(u*u + v*v + w*w)

   xp = (u*(u*x+v*y+w*z) + (x*vw2+u*(-v*y-w*z))*np.cos(theta) + m*(-w*y+v*z)*np.sin(theta))/(m*m)
   yp = (v*(u*x+v*y+w*z) + (y*uw2+v*(-u*x-w*z))*np.cos(theta) + m*( w*x-u*z)*np.sin(theta))/(m*m)
   zp = (w*(u*x+v*y+w*z) + (z*uv2+w*(-u*x-v*y))*np.cos(theta) + m*(-v*x+u*y)*np.sin(theta))/(m*m)

   return np.array([xp, yp, zp])
