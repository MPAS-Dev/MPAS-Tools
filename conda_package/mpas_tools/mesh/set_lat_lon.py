import sys
import numpy as np
import netCDF4
from pyproj import Transformer, CRS
from optparse import OptionParser
from datetime import datetime

from mpas_tools.landice.projections import projections


def main():
    """
    Take MPAS planar grid and populate the lat/lon fields based on a specified
    projection.
    """

    print("== Gathering information.  (Invoke with --help for more details. "
          "All arguments are optional)")
    parser = OptionParser()
    parser.description = "This script populates the MPAS lat and lon fields " \
                         "based on the projection specified by the -p option."
    parser.add_option(
        "-f", "--file", dest="fileInName",
        help="MPAS land ice file name.", default="landice_grid.nc",
        metavar="FILENAME")
    parser.add_option(
        "-p", "--proj", dest="projection",
        help=f"projection used for the data. Valid options are: "
             f"{list(projections.keys())}",
        metavar="PROJ")
    for option in parser.option_list:
        if option.default != ("NO", "DEFAULT"):
            option.help += (" " if option.help else "") + "[default: %default]"
    options, args = parser.parse_args()

    if not options.projection:
        raise ValueError(
            f'Data projection required with -p or --proj command line '
            f'argument. Valid options are: {list(projections.keys())}')

    if not options.fileInName:
        print("No filename specified, so using 'landice_grid.nc'.")
        options.fileInName = 'landice_grid.nc'
    # make a space in stdout before further output
    print('')

    # =================================================

    print(f"Using {options.projection} projection, defined as: "
          f"{projections[options.projection]}")

    # get needed fields
    f = netCDF4.Dataset(options.fileInName, 'r+')
    xCell = f.variables['xCell']
    yCell = f.variables['yCell']
    xVertex = f.variables['xVertex']
    yVertex = f.variables['yVertex']
    xEdge = f.variables['xEdge']
    yEdge = f.variables['yEdge']

    latCellVar = f.variables['latCell']
    lonCellVar = f.variables['lonCell']
    latVertexVar = f.variables['latVertex']
    lonVertexVar = f.variables['lonVertex']
    latEdgeVar = f.variables['latEdge']
    lonEdgeVar = f.variables['lonEdge']

    print("Input file xCell min/max values:", xCell[:].min(), xCell[:].max())
    print("Input file yCell min/max values:", yCell[:].min(), yCell[:].max())

    # make a CRS (coordinate reference system) for projections from Proj
    # string:
    crs_in = CRS.from_proj4(projections[options.projection])
    crs_out = CRS.from_proj4(projections['latlon'])

    # define a transformer
    t = Transformer.from_crs(crs_in, crs_out)

    # populate x,y fields
    # MPAS uses lat/lon in radians, so have pyproj return fields in radians.
    lonCell, latCell = t.transform(xCell[:], yCell[:], radians=True)
    lonVertex, latVertex = t.transform(xVertex[:], yVertex[:], radians=True)
    lonEdge, latEdge = t.transform(xEdge[:], yEdge[:], radians=True)

    # change the longitude convention to use positive values [0 2pi]
    lonCell = np.mod(lonCell, 2.0*np.pi)
    lonVertex = np.mod(lonVertex, 2.0*np.pi)
    lonEdge = np.mod(lonEdge, 2.0*np.pi)

    print(f"Calculated latCell min/max values (radians): "
          f"{latCell.min()}, {latCell.max()}")
    print(f"Calculated lonCell min/max values (radians): "
          f"{lonCell.min()}, {lonCell.max()}")

    latCellVar[:] = latCell
    lonCellVar[:] = lonCell
    latVertexVar[:] = latVertex
    lonVertexVar[:] = lonVertex
    latEdgeVar[:] = latEdge
    lonEdgeVar[:] = lonEdge

    # Update history attribute of netCDF file
    thiscommand = datetime.now().strftime("%a %b %d %H:%M:%S %Y") + ": " + \
        " ".join(sys.argv[:])
    if hasattr(f, 'history'):
        newhist = '\n'.join([thiscommand, getattr(f, 'history')])
    else:
        newhist = thiscommand
        setattr(f, 'history', newhist)

    f.close()

    print("Lat/lon calculations completed.  File has been written.")
