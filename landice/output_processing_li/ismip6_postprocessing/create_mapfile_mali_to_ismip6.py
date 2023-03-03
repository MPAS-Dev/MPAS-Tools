from mpas_tools.scrip.from_mpas import scrip_from_mpas
from subprocess import check_call
import os
import netCDF4
import xarray as xr

def build_mapping_file(mali_mesh_file,
                       mapping_file, res_ismip6_grid,
                       ismip6_grid_file=None,
                       method_remap=None):
    """
    Build a mapping file if it does not exist.
    Mapping file is then used to remap the MALI source file to the
    ISMIP6 ppolarstero grid

    Parameters
    ----------

    mali_mesh_file : str
        mali file

    mapping_file : str
        weights for interpolation from mali_mesh_file to ismip6_grid_file

    res_ismip6_grid: str
        resolution of the ismip6 grid in kilometers
    ismip6_grid_file : str, optional
        The ISMIP6 file if mapping file does not exist

    method_remap : str, optional
        Remapping method used in building a mapping file
    """

    if os.path.exists(mapping_file):
        print(f"Mapping file exists. Not building a new one.")
        return

    if ismip6_grid_file is None:
        raise ValueError("Mapping file does not exist. To build one, ISMIP6 "
                         "grid file with '-f' should be provided. "
                         "Type --help for info")

    if method_remap is None:
        method_remap = "conserve"

    ismip6_scripfile = f"temp_ismip6_{res_ismip6_grid}km_scrip.nc"
    mali_scripfile = "temp_mali_scrip.nc"
    ismip6_projection = "ais-bedmap2"

    # create the ismip6 scripfile if mapping file does not exist
    # this is the projection of ismip6 data for Antarctica
    print(f"Mapping file does not exist. Building one based on "
          f"the input/ouptut meshes")
    print(f"Creating temporary scripfiles "
          f"for ismip6 grid and mali mesh...")

    # create a scripfile for ismip6 grid
    create_atm_scrip(ismip6_grid_file, ismip6_scripfile)

    # create a MALI mesh scripfile
    scrip_from_mpas(mali_mesh_file, mali_scripfile)

    # create a mapping file using ESMF weight gen
    print(f"Creating a mapping file... "
          f"Mapping method used: {method_remap}")

    # generate a mapping file using the scrip files

    # On compute node, ESMG regridder needs to be called with srun
    # On head nodes or local machines it does not.
    # Here, assuming a compute node has hostname starting with nid,
    # which is the case on Cori.  Modify as needed for other machines.
    # Also assuming we can use multiple cores on a compute node.
    hostname = os.uname()[1]
    if hostname.startswith('nid'):
        args = (["srun", "-n", "12", "ESMF_RegridWeightGen",
             "-s", mali_scripfile,
             "-d", ismip6_scripfile,
             "-w", mapping_file,
             "-m", method_remap,
             "-i", "-64bit_offset",
             "--dst_regional", "--src_regional"])
    else:
        args = (["ESMF_RegridWeightGen",
             "-s", mali_scripfile,
             "-d", ismip6_scripfile,
             "-w", mapping_file,
             "-m", method_remap,
             "-i", "-64bit_offset",
             "--dst_regional", "--src_regional"])

    check_call(args)

    # remove the temporary scripfiles once the mapping file is generated
    print(f"Removing the temporary mesh and scripfiles...")
    os.remove(ismip6_scripfile)
    os.remove(mali_scripfile)


def create_atm_scrip(source_grid_file, source_grid_scripfile):
    """
    create a scripfile for the ismip6 atmospheric forcing data.
    Note: the atmospheric forcing data do not have 'x' and 'y' coordinates and
    only have dimensions of them. This function uses 'lat' and 'lon'
    coordinates to generate a scripfile.
    Parameters
    ----------
    source_grid_file : str
        input smb grid file
    source_grid_scripfile : str
        output scrip file of the input smb data
    """

    ds = xr.open_dataset(source_grid_file)
    out_file = netCDF4.Dataset(source_grid_scripfile, 'w')

    if "rlon" in ds and "rlat" in ds:  # this is for RACMO's rotated-pole grid
        nx = ds.sizes["rlon"]
        ny = ds.sizes["rlat"]
    else:
        nx = ds.sizes["x"]
        ny = ds.sizes["y"]
    units = "degrees"

    grid_size = nx * ny

    out_file.createDimension("grid_size", grid_size)
    out_file.createDimension("grid_corners", 4)
    out_file.createDimension("grid_rank", 2)

    # Variables
    grid_center_lat = out_file.createVariable("grid_center_lat", "f8",
                                              ("grid_size",))
    grid_center_lat.units = units
    grid_center_lon = out_file.createVariable("grid_center_lon", "f8",
                                              ("grid_size",))
    grid_center_lon.units = units
    grid_corner_lat = out_file.createVariable("grid_corner_lat", "f8",
                                              ("grid_size", "grid_corners"))
    grid_corner_lat.units = units
    grid_corner_lon = out_file.createVariable("grid_corner_lon", "f8",
                                              ("grid_size", "grid_corners"))
    grid_corner_lon.units = units
    grid_imask = out_file.createVariable("grid_imask", "i4", ("grid_size",))
    grid_imask.units = "unitless"
    out_file.createVariable("grid_dims", "i4", ("grid_rank",))

    out_file.variables["grid_center_lat"][:] = ds.lat.values.flat
    out_file.variables["grid_center_lon"][:] = ds.lon.values.flat
    out_file.variables["grid_dims"][:] = [nx, ny]
    out_file.variables["grid_imask"][:] = 1

    if "lat_bnds" in ds and "lon_bnds" in ds:
        lat_corner = ds.lat_bnds
        if "time" in lat_corner.dims:
            lat_corner = lat_corner.isel(time=0)

        lon_corner = ds.lon_bnds
        if "time" in lon_corner.dims:
            lon_corner = lon_corner.isel(time=0)

        lat_corner = lat_corner.values
        lon_corner = lon_corner.values
    else:
        # this part is used for RACMO as it does not have lat_bnds & lon_bnds
        lat_corner = _unwrap_corners(interp_extrap_corners_2d(ds.lat.values))
        lon_corner = _unwrap_corners(interp_extrap_corners_2d(ds.lon.values))

    grid_corner_lat = lat_corner.reshape((grid_size, 4))
    grid_corner_lon = lon_corner.reshape((grid_size, 4))

    out_file.variables["grid_corner_lat"][:] = grid_corner_lat
    out_file.variables["grid_corner_lon"][:] = grid_corner_lon

    out_file.close()


def _unwrap_corners(in_field):
    """
    Turn a 2D array of corners into an array of rectangular mesh elements
    """
    out_field = np.zeros(((in_field.shape[0] - 1) *
                          (in_field.shape[1] - 1), 4))
    # corners are counterclockwise
    out_field[:, 0] = in_field[0:-1, 0:-1].flat
    out_field[:, 1] = in_field[0:-1, 1:].flat
    out_field[:, 2] = in_field[1:, 1:].flat
    out_field[:, 3] = in_field[1:, 0:-1].flat

    return out_field
