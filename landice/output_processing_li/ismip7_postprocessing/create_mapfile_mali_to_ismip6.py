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
    args = ["create_SCRIP_file_from_planar_rectangular_grid.py",
            "--input", ismip6_grid_file,
            "--scrip", ismip6_scripfile,
            "--proj", ismip6_projection,
            "--rank", "2"]

    check_call(args)

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
