# Create a SCRIP file from an MPAS mesh.
# See for details: http://www.earthsystemmodeling.org/esmf_releases/public/ESMF_5_2_0rp1/ESMF_refdoc/node3.html#SECTION03024000000000000000

import sys
import netCDF4
import numpy as np

from optparse import OptionParser
from mpas_tools.cime.constants import constants

def scrip_from_mpas(mpasFile, scripFile, useLandIceMask=False):
    """
    Create a SCRIP file from an MPAS mesh

    Parameters
    ----------
    mpasFile : str
        The path to a NetCDF file with the MPAS mesh

    scripFile : str
        The path to the output SCRIP file

    useLandIceMask : bool
        Whether to use the landIceMask field for masking
    """
    if useLandIceMask:
        print(" -- Landice Masks are enabled")
    else:
        print(" -- Landice Masks are disabled")

    # make a space in stdout before further output
    print('')
    fin = netCDF4.Dataset(mpasFile, 'r')
    # This will clobber existing files
    fout = netCDF4.Dataset(scripFile, 'w')

    # Get info from input file
    latCell = fin.variables['latCell'][:]
    lonCell = fin.variables['lonCell'][:]
    latVertex = fin.variables['latVertex'][:]
    lonVertex = fin.variables['lonVertex'][:]
    verticesOnCell = fin.variables['verticesOnCell'][:] - 1
    nEdgesOnCell = fin.variables['nEdgesOnCell'][:]
    nCells = len(fin.dimensions['nCells'])
    maxVertices = len(fin.dimensions['maxEdges'])
    areaCell = fin.variables['areaCell'][:]
    sphereRadius = float(fin.sphere_radius)
    on_a_sphere = str(fin.on_a_sphere)

    # check the longitude convention to use positive values [0 2pi]
    if np.any(np.logical_or(lonCell < 0, lonCell > 2.0 * np.pi)):
       raise ValueError("lonCell is not in the desired range (0, 2pi)")

    if np.any(np.logical_or(lonVertex < 0, lonVertex > 2.0 * np.pi)):
       raise ValueError("lonVertex is not in the desired range (0, 2pi)")

    if sphereRadius <= 0:
       sphereRadius =  constants['SHR_CONST_REARTH']
       print(f" -- WARNING: sphereRadius<0 so setting sphereRadius = "
             f"{constants['SHR_CONST_REARTH']}")

    if on_a_sphere == "NO":
        print(" -- WARNING: 'on_a_sphere' attribute is 'NO', which means that "
              "there may be some disagreement regarding area between the "
              "planar (source) and spherical (target) mesh")

    if useLandIceMask:
        landIceMask = fin.variables['landIceMask'][:]
    else:
        landIceMask = None

    # Write to output file
    # Dimensions
    fout.createDimension("grid_size", nCells)
    fout.createDimension("grid_corners", maxVertices)
    fout.createDimension("grid_rank", 1)

    # Variables
    grid_center_lat = fout.createVariable('grid_center_lat', 'f8',
                                          ('grid_size',))
    grid_center_lat.units = 'radians'
    grid_center_lon = fout.createVariable('grid_center_lon', 'f8',
                                          ('grid_size',))
    grid_center_lon.units = 'radians'
    grid_corner_lat = fout.createVariable('grid_corner_lat', 'f8',
                                          ('grid_size', 'grid_corners'))
    grid_corner_lat.units = 'radians'
    grid_corner_lon = fout.createVariable('grid_corner_lon', 'f8',
                                          ('grid_size', 'grid_corners'))
    grid_corner_lon.units = 'radians'
    grid_area = fout.createVariable('grid_area', 'f8', ('grid_size',))
    grid_area.units = 'radian^2'
    grid_imask = fout.createVariable('grid_imask', 'i4', ('grid_size',))
    grid_imask.units = 'unitless'
    grid_dims = fout.createVariable('grid_dims', 'i4', ('grid_rank',))

    grid_center_lat[:] = latCell[:]
    grid_center_lon[:] = lonCell[:]
    # SCRIP uses square radians
    grid_area[:] = areaCell[:]/(sphereRadius**2)
    grid_dims[:] = nCells

    # grid corners:
    grid_corner_lon_local = np.zeros((nCells, maxVertices))
    grid_corner_lat_local = np.zeros((nCells, maxVertices))
    cellIndices = np.arange(nCells)
    lastValidVertex = verticesOnCell[cellIndices, nEdgesOnCell-1]
    for iVertex in range(maxVertices):
        mask = iVertex < nEdgesOnCell
        grid_corner_lat_local[mask, iVertex] = \
            latVertex[verticesOnCell[mask, iVertex]]
        grid_corner_lon_local[mask, iVertex] = \
            lonVertex[verticesOnCell[mask, iVertex]]

        mask = iVertex >= nEdgesOnCell
        grid_corner_lat_local[mask, iVertex] = latVertex[lastValidVertex[mask]]
        grid_corner_lon_local[mask, iVertex] = lonVertex[lastValidVertex[mask]]

    if useLandIceMask:
        # If useLandIceMask are enabled, mask out ocean under land ice.
        grid_imask[:] = 1 - landIceMask[0, :]
    else:
        # If landiceMasks are not enabled, don't mask anything out.
        grid_imask[:] = 1

    grid_corner_lat[:] = grid_corner_lat_local[:]
    grid_corner_lon[:] = grid_corner_lon_local[:]

    print("Input latCell min/max values (radians): {}, {}".format(
        latCell[:].min(), latCell[:].max()))
    print("Input lonCell min/max values (radians): {}, {}".format(
        lonCell[:].min(), lonCell[:].max()))
    print("Calculated grid_center_lat min/max values (radians): {}, {}".format(
         grid_center_lat[:].min(), grid_center_lat[:].max()))
    print("Calculated grid_center_lon min/max values (radians): {}, {}".format(
        grid_center_lon[:].min(), grid_center_lon[:].max()))
    print("Calculated grid_area min/max values (sq radians): {}, {}".format(
        grid_area[:].min(), grid_area[:].max()))

    fin.close()
    fout.close()

    print("Creation of SCRIP file is complete.")


def main():
    print("== Gathering information.  (Invoke with --help for more details. All"
          " arguments are optional)")
    parser = OptionParser()
    parser.description = "This script takes an MPAS grid file and generates " \
                         "a SCRIP grid file."
    parser.add_option("-m", "--mpas", dest="mpasFile",
                      help="MPAS grid file name used as input.",
                      default="grid.nc", metavar="FILENAME")
    parser.add_option("-s", "--scrip", dest="scripFile",
                      help="SCRIP grid file to output.", default="scrip.nc",
                      metavar="FILENAME")
    parser.add_option("-l", "--landice", dest="landiceMasks",
                      help="If flag is on, landice masks will be computed "
                           "and used.",
                      action="store_true")
    for option in parser.option_list:
        if option.default != ("NO", "DEFAULT"):
            option.help += (" " if option.help else "") + "[default: %default]"
    options, args = parser.parse_args()

    if not options.mpasFile:
        sys.exit('Error: MPAS input grid file is required.  Specify with -m '
                 'command line argument.')
    if not options.scripFile:
        sys.exit('Error: SCRIP output grid file is required.  Specify with -s '
                 'command line argument.')

    if not options.landiceMasks:
        options.landiceMasks = False

    scrip_from_mpas(options.mpasFile, options.scripFile, options.landiceMasks)
