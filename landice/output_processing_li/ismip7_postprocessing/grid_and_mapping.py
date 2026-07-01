from mpas_tools.scrip.from_mpas import scrip_from_mpas
from subprocess import check_call
import os
import numpy as np
import netCDF4
from netCDF4 import num2date
import xarray as xr


VALID_EXPERIMENTS = [f"C{i:03d}" for i in range(1, 12)]
VALID_RESOLUTIONS = {
    'AIS': [2, 4, 8, 16],
    'GrIS': [1, 2, 4, 8, 16],
}
CRS_DICT = {
    'AIS': 'epsg:3031',
    'GrIS': 'epsg:3413',
}

# AIS: polar stereographic EPSG:3031, standard parallel 71S, central meridian 0
# Cell centres span -3,040,000 m to 3,040,000 m in both x and y.
AIS_X_MIN = -3040000
AIS_X_MAX = 3040000
AIS_Y_MIN = -3040000
AIS_Y_MAX = 3040000

# GrIS: polar stereographic EPSG:3413, standard parallel 70N, central meridian 315E
# Cell centres span x: -720,000 to 960,000 m; y: -3,450,000 to -570,000 m.
GRIS_X_MIN = -720000
GRIS_X_MAX = 960000
GRIS_Y_MIN = -3450000
GRIS_Y_MAX = -570000


def check_res(res, icesheet):
    """
    Validate the ISMIP7 grid resolution for the given ice sheet.

    Parameters
    ----------
    res : str or int
        Resolution in kilometres to validate.
    icesheet : str
        Ice sheet domain ('AIS' or 'GrIS').

    Raises
    ------
    ValueError
        If the resolution is not valid for the specified ice sheet.
    """
    valid = VALID_RESOLUTIONS[icesheet]
    try:
        res_int = int(res)
    except (ValueError, TypeError):
        raise ValueError(
            f"Resolution '{res}' is not a valid integer. "
            f"Valid resolutions for {icesheet} (km) are: {valid}")
    if res_int not in valid:
        raise ValueError(
            f"Invalid resolution '{res_int}' km for {icesheet}. "
            f"Valid resolutions (km) are: {valid}")
    print(f"Resolution {res_int} km is valid for {icesheet}.")


def check_exp_name(exp):
    """
    Validate the ISMIP7 experiment name.

    Parameters
    ----------
    exp : str
        Experiment name to validate (e.g. 'C001').

    Raises
    ------
    ValueError
        If the experiment name is not in the list of valid experiments.
    """
    if exp not in VALID_EXPERIMENTS:
        raise ValueError(
            f"Invalid experiment name '{exp}'. "
            f"Valid experiments are: {', '.join(VALID_EXPERIMENTS)}")
    print(f"Experiment name '{exp}' is valid.")


def get_time_range(globalstats_files):
    """
    Derive a 'YYYY-YYYY' time-range string from globalStats (1D) files.

    The start and end years are computed from the actual daysSinceStart values
    in the data (not from simulationStartTime, which may be from a historical
    restart). Uses netCDF calendar helpers to properly account for leap years.

    Parameters
    ----------
    globalstats_files : list of str
        Sorted list of MALI globalStats output file paths.

    Returns
    -------
    str
        Year range string in format 'YYYY-YYYY' (e.g., '2000-2014').
    """
    # Get simulationStartTime and calendar from first file
    with xr.open_dataset(globalstats_files[0], decode_cf=False) as ds:
        sim_start = (
            ds['simulationStartTime'].values
            .tobytes().decode('utf-8').strip().strip('\x00')
        )
        # Read and validate calendar type from global attributes
        calendar = ds.attrs.get('config_calendar_type')
        if calendar not in ['noleap', 'gregorian']:
            raise ValueError(
                f"config_calendar_type must be 'noleap' or 'gregorian', "
                f"got '{calendar}'"
            )
        if calendar != 'gregorian':
            print(f"Warning: config_calendar_type is '{calendar}', not 'gregorian'")
    
    # Parse the full datetime (format: YYYY-MM-DD_HH:MM:SS)
    sim_start_date = sim_start.split('_')[0]  # YYYY-MM-DD

    # Get the first and last daysSinceStart values from all files
    with xr.open_mfdataset(
            globalstats_files, combine='nested', concat_dim='Time',
            decode_cf=False, data_vars='minimal',
            coords='minimal', compat='override') as ds:
        days_first = ds['daysSinceStart'].values[0]
        days_last = ds['daysSinceStart'].values[-1]

    # Use netCDF calendar helpers to convert days to proper dates
    units = f'days since {sim_start_date}'
    
    # Convert first and last daysSinceStart to datetime objects
    date_first = num2date(days_first, units=units, calendar=calendar)
    date_last = num2date(days_last, units=units, calendar=calendar)

    # Extract years from the dates
    start_year = date_first.year
    # End year is one less than the calendar year of the final timestamp
    # (assuming final timestamp is at midnight Jan 1 of subsequent year)
    end_year = date_last.year - 1

    return f"{start_year}-{end_year}"


def create_ismip7_grid_file(icesheet, res_km, output_file):
    """
    Create a minimal ISMIP7 standard grid file containing only x and y
    projected coordinate variables (metres).

    Parameters
    ----------
    icesheet : str
        Ice sheet domain: 'AIS' or 'GrIS'.
    res_km : int
        Grid resolution in kilometres.
    output_file : str
        Path for the output NetCDF file.
    """
    res_m = int(res_km) * 1000
    if icesheet == 'AIS':
        x = np.arange(AIS_X_MIN, AIS_X_MAX + res_m, res_m, dtype=float)
        y = np.arange(AIS_Y_MIN, AIS_Y_MAX + res_m, res_m, dtype=float)
    elif icesheet == 'GrIS':
        x = np.arange(GRIS_X_MIN, GRIS_X_MAX + res_m, res_m, dtype=float)
        y = np.arange(GRIS_Y_MIN, GRIS_Y_MAX + res_m, res_m, dtype=float)
    else:
        raise ValueError(f"Unknown icesheet '{icesheet}'.")

    ds = netCDF4.Dataset(output_file, 'w')
    ds.createDimension('x', len(x))
    ds.createDimension('y', len(y))
    xv = ds.createVariable('x', 'f8', ('x',))
    yv = ds.createVariable('y', 'f8', ('y',))
    xv[:] = x
    yv[:] = y
    xv.units = 'm'
    xv.standard_name = 'projection_x_coordinate'
    xv.long_name = 'x'
    yv.units = 'm'
    yv.standard_name = 'projection_y_coordinate'
    yv.long_name = 'y'
    ds.close()
    print(f"Created ISMIP7 grid file: {output_file} "
          f"({len(x)} x {len(y)} cells at {res_km} km)")


def check_ismip7_grid_file(ismip7_grid_file_path, res_ismip7_grid):
    """
    Ensure the ISMIP7 grid file has 'x' and 'y' coordinate variables and that
    its corners match the ISMIP7-required extents.

    Parameters
    ----------
    ismip7_grid_file_path : str
        Path to the ISMIP7 grid file supplied by the user.
    res_ismip7_grid : str
        Resolution of the ISMIP7 grid in kilometres (e.g. '8').

    Returns
    -------
    None
    """
    print("\n---Checking the coordinate variables of the ismip7 grid file---")
    data_ismip7 = netCDF4.Dataset(ismip7_grid_file_path, "r")

    if 'x' in data_ismip7.variables and 'y' in data_ismip7.variables:
        ismip7_grid_file = ismip7_grid_file_path
        print("'x' and 'y' coordinates exist in the file.")
    else:
        data_ismip7.close()
        raise ValueError(
            f"'x' and/or 'y' coordinate variables are missing from the "
            f"ISMIP7 grid file '{ismip7_grid_file_path}'. Please provide a "
            f"grid file that includes 'x' and 'y' coordinate variables.")

    data_ismip7.close()

    print("Checking the grid corners...")
    check_ds = netCDF4.Dataset(ismip7_grid_file, "r")
    x = check_ds.variables["x"]
    y = check_ds.variables["y"]
    if not x[0] == -3040000 or not y[0] == -3040000:
        raise ValueError(
            f"The lower left corner values must be at "
            f"-3040000m and -3040000m. But the values are at "
            f"{x[0]}m and {y[0]}m. Check the value you "
            f"provided for '--res' matches with the resolution of "
            f"the MALI output files. ")
    elif not x[-1] == 3040000 or not y[-1] == 3040000:
        raise ValueError(
            f"The upper right corner values must be at "
            f"3040000m and 3040000m. But the values are at "
            f"{x[-1]}m and {y[-1]}m. Check the value you "
            f"provided for '--res' matches with the resolution of "
            f"the MALI output files. ")
    else:
        print(f"Grid corners are as ismip7-required: "
              f"lower left corner values at {x[0]}m and {y[0]}m, and "
              f"upper right corner values at {x[-1]}m and {y[-1]}m")
    check_ds.close()


def build_mapping_file(mali_mesh_file,
                       mapping_file, res_ismip7_grid,
                       icesheet=None,
                       ismip7_grid_file=None,
                       method_remap=None):
    """
    Build a mapping file if it does not exist.
    Mapping file is then used to remap the MALI source file to the
    ISMIP7 ppolarstero grid

    Parameters
    ----------

    mali_mesh_file : str
        mali file

    mapping_file : str
        weights for interpolation from mali_mesh_file to ismip7_grid_file

    res_ismip7_grid: str
        resolution of the ismip7 grid in kilometers

    icesheet : str, optional
        Ice sheet domain ('AIS' or 'GrIS') for determining projection

    ismip7_grid_file : str, optional
        The ISMIP7 file if mapping file does not exist

    method_remap : str, optional
        Remapping method used in building a mapping file
    """

    if os.path.exists(mapping_file):
        print("Mapping file exists. Not building a new one.")
        return

    if ismip7_grid_file is None:
        raise ValueError("Mapping file does not exist. To build one, ISMIP7 "
                         "grid file with '-f' should be provided. "
                         "Type --help for info")

    if method_remap is None:
        method_remap = "conserve"

    ismip7_scripfile = f"temp_ismip7_{res_ismip7_grid}km_scrip.nc"
    mali_scripfile = "temp_mali_scrip.nc"
    if icesheet == 'GrIS':
        ismip7_projection = "gis-gimp"
    elif icesheet == 'AIS':
        ismip7_projection = "ais-bedmap2"
    else:
        raise ValueError(
            f"Unknown icesheet '{icesheet}'. Must be 'AIS' or 'GrIS'.")

    # create the ismip7 scripfile if mapping file does not exist
    # this is the projection of ismip7 data for Antarctica
    print("Mapping file does not exist. Building one based on "
          "the input/ouptut meshes")
    print("Creating temporary scripfiles "
          "for ismip7 grid and mali mesh...")

    # create a scripfile for ismip7 grid
    script_path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "..", "..", "..",
        "conda_package", "mesh_tools", "create_SCRIP_files",
        "create_SCRIP_file_from_planar_rectangular_grid.py"
    )
    script_path = os.path.abspath(script_path)
    args = ["python", script_path,
            "--input", ismip7_grid_file,
            "--scrip", ismip7_scripfile,
            "--proj", ismip7_projection,
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
                 "-d", ismip7_scripfile,
                 "-w", mapping_file,
                 "-m", method_remap,
                 "-i", "-64bit_offset",
                 "--dst_regional", "--src_regional"])
    else:
        args = (["ESMF_RegridWeightGen",
                 "-s", mali_scripfile,
                 "-d", ismip7_scripfile,
                 "-w", mapping_file,
                 "-m", method_remap,
                 "-i", "-64bit_offset",
                 "--dst_regional", "--src_regional"])

    print(f"Running remapping command: {' '.join(args)}")
    check_call(args)

    # remove the temporary scripfiles once the mapping file is generated
    print("Removing the temporary mesh and scripfiles...")
    os.remove(ismip7_scripfile)
    os.remove(mali_scripfile)
