from netCDF4 import Dataset
import numpy as np
import math

degreesToRadians = math.pi / 180.0
radiansToDegrees = 180.0 / math.pi

#---------------------------------------------------------------------

def write_scrip_file(scripFilename, title, nCells, maxEdges, latCell, lonCell, corner_lat, corner_lon, mask=None):
    """
    A low-level function for writing a SCRIP file for the given MPAS-Seaice
    mesh

    Parameters
    ----------
    scripFilename : str
        The name of the resulting SCRIP file

    title : str
        A string to include as the ``title`` attribute in the SCRIP file

    nCells : int
        The number of cells in the mesh

    maxEdges : int
        The maximum number of edges/vertices on a cell in the mesh

    latCell : numpy.ndarray
        The latitude (in radians) of cell centers

    lonCell : numpy.ndarray
        The longitude (in radians) of cell centers

    corner_lat : numpy.ndarray
        The latitude (in radians) of cell vertices

    corner_lon : numpy.ndarray
        The longitude (in radians) of cell vertices

    mask  : numpy.ndarray, optional
        A mask of where cells are valid
    """

    # output mesh file
    scripFile = Dataset(scripFilename, "w", format="NETCDF3_CLASSIC")

    scripFile.title = title.strip()

    scripFile.createDimension("grid_size", nCells)
    scripFile.createDimension("grid_corners", maxEdges)
    scripFile.createDimension("grid_rank", 1)

    grid_dims       = scripFile.createVariable("grid_dims",       "i", dimensions=["grid_rank"])
    grid_center_lat = scripFile.createVariable("grid_center_lat", "d", dimensions=["grid_size"])
    grid_center_lon = scripFile.createVariable("grid_center_lon", "d", dimensions=["grid_size"])
    grid_imask      = scripFile.createVariable("grid_imask",      "i", dimensions=["grid_size"])
    grid_corner_lat = scripFile.createVariable("grid_corner_lat", "d", dimensions=["grid_size","grid_corners"])
    grid_corner_lon = scripFile.createVariable("grid_corner_lon", "d", dimensions=["grid_size","grid_corners"])

    grid_center_lat.units = "radians"
    grid_center_lon.units = "radians"
    grid_imask.units      = "unitless"
    grid_corner_lat.units = "radians"
    grid_corner_lon.units = "radians"

    grid_dims      [:]   = [nCells]
    grid_center_lat[:]   = latCell
    grid_center_lon[:]   = lonCell
    if (mask is None):
        grid_imask [:]   = np.ones(nCells,dtype="i")
    else:
        grid_imask [:]   = mask[:]
    grid_corner_lat[:,:] = corner_lat[:,:]
    grid_corner_lon[:,:] = corner_lon[:,:]

    scripFile.close()

#------------------------------------------------------------------------------------

def write_2D_scripfile(filenameScripOut, scripTitle, nColumns, nRows, latsCentre, lonsCentre, latsVertex, lonsVertex, degrees=False):
    """
    Write a SCRIP file for the given 2D grid

    Parameters
    ----------
    filenameScripOut : str
        The name of the resulting SCRIP file

    scripTitle : str
        A string to include as the ``title`` attribute in the SCRIP file

    nColumns : int
        The number of columns in the grid

    nRows : int
        The number of rows in the grid

    latsCentre : numpy.ndarray
        The latitude (in radians) of cell centers

    lonsCentre : numpy.ndarray
        The longitude (in radians) of cell centers

    latsVertex : numpy.ndarray
        The latitude (in radians) of cell vertices

    lonsVertex : numpy.ndarray
        The longitude (in radians) of cell vertices

    degrees  : bool, optional
        Whether the latitude and longitude variables are in degrees (as opposed
        to radians)
    """

    if (degrees):
        scaling = np.pi/180.
    else:
        scaling = 1.0

    fileOut = Dataset(filenameScripOut,"w",format="NETCDF3_CLASSIC")

    fileOut.title = scripTitle

    grid_size    = nRows*nColumns
    grid_corners = 4
    grid_rank    = 2

    fileOut.createDimension("grid_size",    grid_size)
    fileOut.createDimension("grid_corners", grid_corners)
    fileOut.createDimension("grid_rank",    grid_rank)

    grid_dimsVar       = fileOut.createVariable("grid_dims"      , "i", dimensions=("grid_rank"))
    grid_center_latVar = fileOut.createVariable("grid_center_lat", "d", dimensions=("grid_size"))
    grid_center_lonVar = fileOut.createVariable("grid_center_lon", "d", dimensions=("grid_size"))
    grid_imaskVar      = fileOut.createVariable("grid_imask"     , "i", dimensions=("grid_size"))
    grid_corner_latVar = fileOut.createVariable("grid_corner_lat", "d", dimensions=("grid_size","grid_corners"))
    grid_corner_lonVar = fileOut.createVariable("grid_corner_lon", "d", dimensions=("grid_size","grid_corners"))

    grid_dims       = np.zeros(grid_rank, dtype="i")
    grid_center_lat = np.zeros(grid_size, dtype="d")
    grid_center_lon = np.zeros(grid_size, dtype="d")
    grid_imask      = np.zeros(grid_size, dtype="i")
    grid_corner_lat = np.zeros((grid_size,grid_corners), dtype="d")
    grid_corner_lon = np.zeros((grid_size,grid_corners), dtype="d")

    grid_dims[0] = nColumns
    grid_dims[1] = nRows

    for iRow in range(0,nRows):
        for iColumn in range(0,nColumns):

            i = iColumn + iRow*nColumns

            grid_center_lat[i] = latsCentre[iRow,iColumn] * scaling
            grid_center_lon[i] = lonsCentre[iRow,iColumn] * scaling

            grid_imask     [i] = 1

            for iVertex in range(0,4):
                grid_corner_lat[i,iVertex] = latsVertex[iRow,iColumn,iVertex] * scaling
                grid_corner_lon[i,iVertex] = lonsVertex[iRow,iColumn,iVertex] * scaling

    grid_dimsVar      [:] = grid_dims
    grid_center_latVar[:] = grid_center_lat
    grid_center_lonVar[:] = grid_center_lon
    grid_imaskVar     [:] = grid_imask
    grid_corner_latVar[:] = grid_corner_lat
    grid_corner_lonVar[:] = grid_corner_lon

    grid_center_latVar.units = "radians"
    grid_center_lonVar.units = "radians"
    grid_imaskVar.units = "unitless"
    grid_corner_latVar.units = "radians"
    grid_corner_lonVar.units = "radians"

    fileOut.close()

#---------------------------------------------------------------------

def make_mpas_scripfile_on_cells(meshFilename, scripFilename, title):
    """
    Write a SCRIP file for cel quantities on the given MPAS-Seaice mesh

    Parameters
    ----------
    meshFilename : str
        The name of a file containing the MPAS-Seaice mesh

    scripFilename : str
        The name of the resulting SCRIP file

    title : str
        A string to include as the ``title`` attribute in the SCRIP file
    """

    # input mesh data
    meshFile = Dataset(meshFilename, "r")

    nCells   = len(meshFile.dimensions["nCells"])
    maxEdges = len(meshFile.dimensions["maxEdges"])

    latCell = meshFile.variables["latCell"][:]
    lonCell = meshFile.variables["lonCell"][:]

    latVertex = meshFile.variables["latVertex"][:]
    lonVertex = meshFile.variables["lonVertex"][:]

    nEdgesOnCell   = meshFile.variables["nEdgesOnCell"][:]
    verticesOnCell = meshFile.variables["verticesOnCell"][:]

    meshFile.close()

    # make corner arrays
    corner_lat = np.zeros((nCells,maxEdges), dtype="d")
    corner_lon = np.zeros((nCells,maxEdges), dtype="d")

    for iCell in range(0, nCells):

        for iVertexOnCell in range(0, nEdgesOnCell[iCell]):

            iVertex = verticesOnCell[iCell,iVertexOnCell] - 1

            corner_lat[iCell,iVertexOnCell] = latVertex[iVertex]
            corner_lon[iCell,iVertexOnCell] = lonVertex[iVertex]

        for iVertexOnCell in range(nEdgesOnCell[iCell], maxEdges):

            corner_lat[iCell,iVertexOnCell] = corner_lat[iCell,nEdgesOnCell[iCell]-1]
            corner_lon[iCell,iVertexOnCell] = corner_lon[iCell,nEdgesOnCell[iCell]-1]

    # create the scrip file
    write_scrip_file(scripFilename, title, nCells, maxEdges, latCell, lonCell, corner_lat, corner_lon)

#---------------------------------------------------------------------

def make_mpas_scripfile_on_vertices(meshFilename, scripFilename, title):
    """
    Write a SCRIP file for vertex quantities on the given MPAS-Seaice mesh

    Parameters
    ----------
    meshFilename : str
        The name of a file containing the MPAS-Seaice mesh

    scripFilename : str
        The name of the resulting SCRIP file

    title : str
        A string to include as the ``title`` attribute in the SCRIP file
    """

    # input mesh data
    meshFile = Dataset(meshFilename, "r")

    nVertices    = len(meshFile.dimensions["nVertices"])
    vertexDegree = len(meshFile.dimensions["vertexDegree"])

    latCell = meshFile.variables["latCell"][:]
    lonCell = meshFile.variables["lonCell"][:]

    latVertex = meshFile.variables["latVertex"][:]
    lonVertex = meshFile.variables["lonVertex"][:]

    cellsOnVertex  = meshFile.variables["cellsOnVertex"][:]

    meshFile.close()

    # make corner arrays
    corner_lat = np.zeros((nVertices,vertexDegree), dtype="d")
    corner_lon = np.zeros((nVertices,vertexDegree), dtype="d")

    for iVertex in range(0, nVertices):

        for iCellOnVertex in range(0, vertexDegree):

            iCell = cellsOnVertex[iVertex,iCellOnVertex] - 1

            if (iCell != -1):

                corner_lat[iVertex,iCellOnVertex] = latCell[iCell]
                corner_lon[iVertex,iCellOnVertex] = lonCell[iCell]

            else:

                corner_lat[iVertex,iCellOnVertex], corner_lon[iVertex,iCellOnVertex] = _estimate_missing_cell_latlon(
                        iVertex, iCellOnVertex, vertexDegree, cellsOnVertex, latCell, lonCell, latVertex, lonVertex)

    # create the scrip file
    write_scrip_file(scripFilename, title, nVertices, vertexDegree, latVertex, lonVertex, corner_lat, corner_lon)


#---------------------------------------------------------------------
# Private functions
#---------------------------------------------------------------------

def _rotate_about_vector(vectorToRotate, vectorAxis, rotationAngle):
    """
    Rotate the given Cartesian vector around the given axis by the given
    angle.

    Parameters
    ----------
    vectorToRotate : numpy.ndarray
        The 3-element vector to rotate

    vectorAxis : numpy.ndarray
        The 3-element axis of rotation

    rotationAngle : float
        The angle of rotation (in radians)

    Returns
    -------
    rotatedVector : numpy.ndarray
        The rotated vector

    """

    # http://ksuweb.kennesaw.edu/~plaval//math4490/rotgen.pdf

    r = vectorAxis / np.linalg.norm(vectorAxis)

    C = math.cos(rotationAngle)
    S = math.sin(rotationAngle)
    t = 1.0 - C

    rot = np.zeros((3,3))

    rot[0,0] = t * r[0] * r[0] + C
    rot[0,1] = t * r[0] * r[1] - S * r[2]
    rot[0,2] = t * r[0] * r[2] + S * r[1]

    rot[1,0] = t * r[0] * r[1] + S * r[2]
    rot[1,1] = t * r[1] * r[1] + C
    rot[1,2] = t * r[1] * r[2] - S * r[0]

    rot[2,0] = t * r[0] * r[2] - S * r[1]
    rot[2,1] = t * r[1] * r[2] + S * r[0]
    rot[2,2] = t * r[2] * r[2] + C

    # get relative position of missing cell
    rotatedVector = np.matmul(rot, vectorToRotate)

    return rotatedVector

#---------------------------------------------------------------------

def _wrap_index(i, n):

    return i % n

#---------------------------------------------------------------------

def _rotate_about_vertex(latCell, lonCell, latVertex, lonVertex, angle):

    vertexPosition = _get_position_from_lat_lon(latVertex, lonVertex)
    cellPosition   = _get_position_from_lat_lon(latCell, lonCell)

    vertexToCellVector = cellPosition - vertexPosition

    vertexToCellVectorMagnitude = _vector_magnitude(vertexToCellVector)

    vertexToCellVectorNormalized = _normalize_vector(vertexToCellVector)

    perpendicularVector = _cross_product(vertexPosition, vertexToCellVector)

    perpendicularVectorNormalized = _normalize_vector(perpendicularVector)

    rotatedVectorNormalized = math.cos(angle * degreesToRadians) * vertexToCellVectorNormalized + \
                              math.sin(angle * degreesToRadians) * perpendicularVectorNormalized

    rotatedVector = rotatedVectorNormalized * vertexToCellVectorMagnitude

    newPositionVector = vertexPosition + rotatedVector

    latRotated, lonRotated = _get_lat_lon_from_position(newPositionVector)

    return latRotated, lonRotated

#---------------------------------------------------------------------

def _estimate_missing_cell_latlon(iVertex, iCellOnVertexNeeded, vertexDegree, cellsOnVertex, latCell, lonCell, latVertex, lonVertex):
    """
    Estimate the latitude and longitude of a "missing" neighbor cell that has
    been culled from the MPAS-Seaice mesh

    Parameters
    ----------
    iVertex : int
        The index of the vertex with the missing neighbor

    iCellOnVertexNeeded : int
        The local index on the vertex of the missing cell

    vertexDegree : int
        The maximum number of cells neighboring each vertex

    cellsOnVertex : numpy.ndarray
        The indices of cells neighboring each vertex in the mesh

    latCell : numpy.ndarray
        The latitude (in radians) of cell centers

    lonCell : numpy.ndarray
        The longitude (in radians) of cell centers

    latVertex : numpy.ndarray
        The latitude (in radians) of vertices

    lonVertex : numpy.ndarray
        The longitude (in radians) of vertices

    Returns
    -------
    lat : float
        The approximate latitude of the missing cell center

    lon : float
        The approximate longitude of the missing cell center

    Raises
    ------
    ValueError
        If no reasonable location for the missing cell can be found
    """

    # firstly we find a cell on vertex that exists
    nRotations = 0

    for iCellOnVertex in range(0,vertexDegree):

        iCell = cellsOnVertex[iVertex, _wrap_index(iCellOnVertexNeeded - iCellOnVertex - 1, vertexDegree)] - 1
        nRotations = nRotations + 1

        if (iCell != -1):

            # find relative vector to cell we have which needs to be rotated
            cellPosition   = _get_position_from_lat_lon(latCell  [iCell],   lonCell  [iCell])
            vertexPosition = _get_position_from_lat_lon(latVertex[iVertex], lonVertex[iVertex])

            cellVector = cellPosition - vertexPosition

            # rotate cell vector around vertex vector
            theta = nRotations * ((2.0 * math.pi) / float(vertexDegree))
            missingCellVector = _rotate_about_vector(cellVector, vertexPosition, theta)

            # get absolute position
            missingCellPosition = vertexPosition + missingCellVector

            # get missing cell lat, lon
            lat, lon = _get_lat_lon_from_position(missingCellPosition)

            return lat, lon

    raise ValueError("Can't find position!")

#---------------------------------------------------------------------

def _get_position_from_lat_lon(lat, lon):

    position = np.zeros(3)

    position[0] = math.cos(lat) * math.cos(lon)
    position[1] = math.cos(lat) * math.sin(lon)
    position[2] = math.sin(lat)

    return position

#---------------------------------------------------------------------

def _get_lat_lon_from_position(position):

    lat = math.asin(position[2])
    lon = math.atan2(position[1], position[0])

    return lat, lon

#---------------------------------------------------------------------

def _vector_magnitude(vector):

    magnitude = math.sqrt(math.pow(vector[0],2) + \
                          math.pow(vector[1],2) + \
                          math.pow(vector[2],2))

    return magnitude

#---------------------------------------------------------------------

def _normalize_vector(vector):

    magnitude = _vector_magnitude(vector)

    normalizedVector = vector / magnitude

    return normalizedVector

#---------------------------------------------------------------------

def _cross_product(u, v):

    cp = np.zeros(3)

    cp[0] = u[1] * v[2] - u[2] * v[1]
    cp[1] = u[2] * v[0] - u[0] * v[2]
    cp[2] = u[0] * v[1] - u[1] * v[0]

    return cp
