"""
This script has functions that are needed to post-process and write state
output variables from ISMIP7 simulations.
"""

from netCDF4 import Dataset, num2date, date2num, default_fillvals
import xarray as xr
import numpy as np
import os
from subprocess import check_call
from validate import validate_mali_files
from output_naming import build_output_filename
from time_utils import get_calendar, days_to_ismip7_time
from netcdf_utils import to_masked, zero_where_invalid


# ISMIP7 variables for which missing values are not permitted: where there is
# no ice the value is 0.
ZERO_WHERE_NO_ICE_VARIABLES = ('lithk', 'sftgif', 'sftgrf', 'sftflf')

EXPECTED_STATE_VARIABLES = [
    'daysSinceStart', 'simulationStartTime',
    'cellMask', 'basalTemperature', 'betaSolve',
    'uReconstructX', 'uReconstructY', 'upperSurface',
]

# Keep only variables needed downstream by generate_output_2d_state_vars
# and write_netcdf_2d_state_vars.
REQUIRED_STATE_REMAPPING_VARIABLES = [
    'daysSinceStart', 'simulationStartTime',
    'thickness', 'upperSurface', 'lowerSurface', 'bedTopography',
    'uReconstructX_sfc', 'uReconstructY_sfc',
    'uReconstructX_base', 'uReconstructY_base',
    'xvelmean', 'yvelmean', 'surfaceTemperature',
    'litempbotgr', 'litempbotfl', 'strbasemag',
    'sftgif', 'sftgrf', 'sftflf',
]

# Written only if present in the MALI output (runs with a GIA model such as
# FastIsostasy or the 1D sea-level model that outputs geoidChange).
OPTIONAL_STATE_REMAPPING_VARIABLES = ['geoidChange']


def check_state_files(files):
    """
    Validate a list of MALI state output files before processing.

    See validate.validate_mali_files for full details of checks performed.
    """
    validate_mali_files(files, EXPECTED_STATE_VARIABLES, label='state')


def process_state_vars(files, tmp_file):
    """
    files: list of MALI state output file paths
    tmp_file: temporary output file name
    """

    inputfile_state_vars = xr.open_mfdataset(files, combine='nested',
                                             concat_dim='Time',
                                             engine='netcdf4',
                                             decode_cf=False,
                                             data_vars='minimal',
                                             coords='minimal',
                                             compat='override')
    # Delete the 'units' attr on daysSinceStart to prevent xarray from
    # encoding it as a timedelta type (corrupts values after ~250 years).
    if 'units' in inputfile_state_vars['daysSinceStart'].attrs:
        del inputfile_state_vars['daysSinceStart'].attrs['units']

    # Skip the first time level if it is at the simulation start time
    days_first = inputfile_state_vars['daysSinceStart'].values[0]
    if inputfile_state_vars.sizes['Time'] > 0 and np.isclose(days_first, 0.0, rtol=0.0):
        print("Skipping state data at simulation start time.")
        inputfile_state_vars = inputfile_state_vars.isel(Time=slice(1, None))

    # get the mesh description data
    # number of vertical interfaces (uReconstructX/Y, which are required, are
    # on nVertInterfaces; the files may have no nVertLevels dimension)
    nInterface = inputfile_state_vars.sizes['nVertInterfaces']
    cellMask = inputfile_state_vars['cellMask'][:, :]
    basalTemperature = inputfile_state_vars['basalTemperature'][:, :]
    betaSolve = inputfile_state_vars['betaSolve'][:, :]

    inputfile_state_vars['litempbotfl'] = basalTemperature * \
        (cellMask[:, :] & 4) / 4
    inputfile_state_vars['litempbotgr'] = basalTemperature * \
        (1 - (cellMask[:, :] & 4) / 4)

    uxsurf = inputfile_state_vars['uReconstructX'][:, :, 0]
    uysurf = inputfile_state_vars['uReconstructY'][:, :, 0]
    uxbase = inputfile_state_vars['uReconstructX'][:, :, nInterface - 1]
    uybase = inputfile_state_vars['uReconstructY'][:, :, nInterface - 1]
    inputfile_state_vars['uReconstructX_sfc'] = uxsurf
    inputfile_state_vars['uReconstructY_sfc'] = uysurf
    inputfile_state_vars['uReconstructX_base'] = uxbase
    inputfile_state_vars['uReconstructY_base'] = uybase

    inputfile_state_vars['upperSurface'] = np.maximum(
        0.0, inputfile_state_vars['upperSurface'])

    floating_mask = (cellMask[:, :] & 4) / 4
    dynamic_mask = (cellMask[:, :] & 2) / 2
    grounded_mask = (cellMask[:, :] * 0 + 1) - floating_mask

    # floating and dynamic
    inputfile_state_vars['sftflf'] = floating_mask * dynamic_mask
    # grounded: not-floating and dynamic
    inputfile_state_vars['sftgrf'] = grounded_mask * dynamic_mask
    # dynamic ice
    inputfile_state_vars['sftgif'] = dynamic_mask

    speed_base = (uxbase[:, :] ** 2 + uybase[:, :] ** 2) ** 0.5
    seconds_per_year = 3600.0 * 24.0 * 365.0
    inputfile_state_vars['strbasemag'] = (
        betaSolve[:, :] * speed_base * seconds_per_year * grounded_mask *
        dynamic_mask
    )
    # Compliance checker requires this limit.
    # Larger values are not meaningful so no qualms about it.
    inputfile_state_vars['strbasemag'] = inputfile_state_vars['strbasemag'].clip(max=300000.0)

    # Compliance checker requires this limit. 
    # Vinson Massif is 4892 m, so presumably this limit will be increased
    # in the checker.  MJH put in an issue for that to them.
    inputfile_state_vars['bedTopography'] = inputfile_state_vars['bedTopography'].clip(max=4000.0)

    missing = [
        var for var in REQUIRED_STATE_REMAPPING_VARIABLES
        if var not in inputfile_state_vars
    ]
    if missing:
        raise ValueError(
            "Processed state dataset is missing required output variables: "
            f"{missing}"
        )

    keep = list(REQUIRED_STATE_REMAPPING_VARIABLES)
    for var in OPTIONAL_STATE_REMAPPING_VARIABLES:
        if var in inputfile_state_vars:
            keep.append(var)
        else:
            print(f"Optional variable {var} not found in state files; "
                  "the corresponding ISMIP7 variable will not be written.")
    output_state_vars = inputfile_state_vars[keep]
    output_state_vars.to_netcdf(tmp_file)
    output_state_vars.close()
    inputfile_state_vars.close()


def _finalize(values, ismip7_var_name):
    """Fill value for missing data, or 0 for variables that must be 0 where
    there is no ice."""
    if ismip7_var_name in ZERO_WHERE_NO_ICE_VARIABLES:
        return zero_where_invalid(values)
    return to_masked(values)


def write_netcdf_2d_state_vars(
        mali_var_name,
        ismip7_var_name,
        var_std_name,
        var_units,
        var_varname,
        remapped_mali_outputfile,
        ismip7_grid_file,
        output_path,
        metadata,
        extra_var_attrs=None,
        calendar='gregorian'):
    """
    mali_var_name: variable name on MALI side
    ismip7_var_name: variable name required by ISMIP7
    var_std_name: standard variable name
    var_units: variable units
    var_varname: variable variable name
    remapped_mali_outputfile: mali state file remapped on the ISMIP7 grid
    ismip7_grid_file: original ISMIP7 file
    exp: experiment name
    output_path: output path to which the output files will be saved
    extra_var_attrs: optional dict of additional attributes to set on the
        output variable (e.g., reference_geoid for deltag)
    """

    data_ismip7 = Dataset(ismip7_grid_file, 'r')
    var_x = data_ismip7.variables['x'][:]
    var_y = data_ismip7.variables['y'][:]

    data = Dataset(remapped_mali_outputfile, 'r')
    data.set_auto_mask(False)
    simulationStartTime = data.variables['simulationStartTime'][:].tobytes(
    ).decode('utf-8').strip().strip('\x00')
    simulationStartDate = simulationStartTime.split("_")[0]
    # Convert from the MALI model-calendar time axis to Gregorian days since
    # 1850-01-01
    daysSinceStart = data.variables['daysSinceStart'][:]
    daysSinceStart = days_to_ismip7_time(daysSinceStart, simulationStartDate,
                                         calendar)
    var_sftgif = data.variables['sftgif'][:, :, :]
    var_sftgrf = data.variables['sftgrf'][:, :, :]
    var_sftflf = data.variables['sftflf'][:, :, :]
    var_mali = data.variables[mali_var_name][:, :, :]
    var_mali[np.where(abs(var_mali + 1e34) < 1e33)] = np.nan
    timeSteps, latN, lonN = np.shape(var_mali)

    output_filename = build_output_filename(
        output_path,
        ismip7_var_name,
        metadata,
    )
    dataOut = Dataset(output_filename, 'w', format='NETCDF4_CLASSIC')
    dataOut.createDimension('time', None)  # unlimited dimension (record)
    dataOut.createDimension('x', lonN)
    dataOut.createDimension('y', latN)
    dataValues = dataOut.createVariable(ismip7_var_name, 'f4',
                                        ('time', 'y', 'x'),
                                        fill_value=default_fillvals['f4'])
    xValues = dataOut.createVariable('x', 'f4', ('x',),
                                     fill_value=default_fillvals['f4'])
    yValues = dataOut.createVariable('y', 'f4', ('y',),
                                     fill_value=default_fillvals['f4'])
    timeValues = dataOut.createVariable('time', 'f4', ('time',),
                                        fill_value=default_fillvals['f4'])
    timeValues[:] = daysSinceStart
    for i in range(timeSteps):
        if ismip7_var_name == 'sftgif':
            dataValues[i, :, :] = _finalize(var_mali[i, :, :], ismip7_var_name)
        else:
            if ismip7_var_name == 'litempbotgr':
                # todo: get masking division to work right for litempbotgr and litempbotfl
                # for now, only calculate these var for full grounded or full floating cells
                mask = (var_sftgrf[i, :, :])
                with np.errstate(divide='ignore', invalid='ignore'):
                    var_mali[i, :, :] = var_mali[i, :, :] / mask
                mask = (var_sftgrf[i, :, :] > 0.99)
            elif ismip7_var_name == 'litempbotfl':
                mask = (var_sftflf[i, :, :])
                with np.errstate(divide='ignore', invalid='ignore'):
                    var_mali[i, :, :] = var_mali[i, :, :] / mask
                mask = (var_sftflf[i, :, :] > 0.99)
            elif ismip7_var_name in ('topg', 'deltag'):
                # don't mask topg or the geoid anomaly (defined everywhere)
                mask = np.ones(var_mali.shape[1:])
            else:
                mask = var_sftgif[i, :, :]
            tmp = var_mali[i, :, :]
            tmp[mask == 0] = np.nan
            dataValues[i, :, :] = _finalize(tmp, ismip7_var_name)

    xValues[:] = var_x
    yValues[:] = var_y

    dataValues.standard_name = var_std_name
    dataValues.units = var_units
    if extra_var_attrs:
        for key, value in extra_var_attrs.items():
            setattr(dataValues, key, value)
    timeValues.units = 'days since 1850-01-01'
    timeValues.calendar = 'standard'
    timeValues.standard_name = 'time'
    timeValues.long_name = 'time'
    xValues.units = 'm'
    xValues.standard_name = 'x'
    xValues.long_name = 'x'
    yValues.units = 'm'
    yValues.standard_name = 'y'
    yValues.long_name = 'y'
    # Write all metadata as global attributes
    for key, value in metadata.items():
        if key != 'time_range':  # time_range is added by caller
            setattr(dataOut, key, value)
    dataOut.VARIABLE = var_varname
    dataOut.close()


def process_state_pipeline(state_files, mapping_file, ismip7_grid_file,
                           output_path, metadata,
                           reference_geoid='EIGEN-6C4'):
    """
    Full state-variable processing pipeline: validate, adjust, remap, write.

    Parameters
    ----------
    state_files : list of str
        Sorted list of MALI state output file paths.
    mapping_file : str
        Path to the ESMF mapping/weights file.
    ismip7_grid_file : str
        Path to the ISMIP7 grid file.
    output_path : str
        Directory for output files.
    metadata : dict
        Submission metadata (exp, icesheet, authors, group, model, date, ...).
        time_range will be used as already computed by the main driver.
    reference_geoid : str
        Value of the reference_geoid attribute written on deltag:
        'EIGEN-6C4' (BedMachine), 'GL04C' (Bedmap), or
        'custom, defined by ref_geoid variable'.
    """
    metadata = metadata.copy()

    check_state_files(state_files)

    with xr.open_dataset(state_files[0], decode_cf=False) as ds_first:
        calendar = get_calendar(ds_first.attrs)
    print(f"Model calendar: {calendar}")

    print("Calculating needed state file adjustments.")
    tmp_file = "tmp_state.nc"
    process_state_vars(state_files, tmp_file)

    processed_and_remapped_file_state = 'processed_and_remapped_state.nc'
    print("Remapping state file.")
    check_call(["ncremap",
                "-i", tmp_file,
                "-o", processed_and_remapped_file_state,
                "-m", mapping_file,
                "-P", "mpas"])

    print("Writing processed and remapped state fields to ISMIP7 file format.")
    generate_output_2d_state_vars(processed_and_remapped_file_state,
                                  ismip7_grid_file, output_path, metadata,
                                  reference_geoid=reference_geoid,
                                  calendar=calendar)

    os.remove(tmp_file)
    os.remove(processed_and_remapped_file_state)


def generate_output_2d_state_vars(file_remapped_mali_state,
                                  ismip7_grid_file, output_path, metadata,
                                  reference_geoid='EIGEN-6C4',
                                  calendar='gregorian'):
    """
    file_remapped_mali_state: output files on mali mesh remapped
    on the ismip7 grid
    ismip7_grid_file: ismip7 original file
    exp: ISMIP7 experiment name
    output_path: path to which the final output files are saved
    """

    # ----------- lithk ------------------
    print("* Writing lithk")
    write_netcdf_2d_state_vars('thickness', 'lithk', 'land_ice_thickness',
                               'm', 'Ice thickness',
                               file_remapped_mali_state,
                               ismip7_grid_file, output_path, metadata,
                               calendar=calendar)

    # ----------- orog ------------------
    print("* Writing orog")
    write_netcdf_2d_state_vars('upperSurface', 'orog', 'surface_altitude', 'm',
                               'Surface elevation',
                               file_remapped_mali_state,
                               ismip7_grid_file, output_path, metadata,
                               calendar=calendar)

    # ----------- base ------------------
    print("* Writing base")
    write_netcdf_2d_state_vars('lowerSurface', 'base', 'base_altitude', 'm',
                               'Base elevation',
                               file_remapped_mali_state,
                               ismip7_grid_file, output_path, metadata,
                               calendar=calendar)

    # ----------- topg ------------------
    print("* Writing topg")
    write_netcdf_2d_state_vars(
        'bedTopography',
        'topg',
        'bedrock_altitude',
        'm',
        'Bedrock elevation',
        file_remapped_mali_state,
        ismip7_grid_file,
        output_path,
        metadata,
        calendar=calendar)

    # ----------- deltag ------------------
    # Anomaly in geoid height from the reference geoid (ISMIP7 GIA protocol,
    # Table 1).  Only written if the MALI output includes geoidChange (a GIA
    # model that computes geoid changes, e.g., FastIsostasy with SSH feedback
    # or the 1D sea-level model).  MALI initializes geoidChange to zero at the
    # start of the simulation, so it is the change relative to the initial
    # geoid, which is the reference geoid of the bed dataset used (BedMachine:
    # EIGEN-6C4).  refgeoid is not written because it is only required if the
    # reference geoid is not EIGEN-6C4 or GL04C.
    with Dataset(file_remapped_mali_state, 'r') as ds_check:
        have_geoid_change = 'geoidChange' in ds_check.variables
    if have_geoid_change:
        print("* Writing deltag")
        write_netcdf_2d_state_vars(
            'geoidChange',
            'deltag',
            'geopotential_height_anomaly',
            'm',
            'Anomaly in geopotential height from reference geoid',
            file_remapped_mali_state,
            ismip7_grid_file,
            output_path,
            metadata,
            extra_var_attrs={'reference_geoid': reference_geoid},
            calendar=calendar)
    else:
        print("* Skipping deltag (geoidChange not in MALI output)")

    # ----------- hfgeoubed------------------
    # Note: even though this is a flux variable, we are taking a snapshot of it
    # as it does not change with time.
    # Uncomment once basalHeatFlux is outputted in the output stream.
    # write_netcdf_2d_state_vars('basalHeatFlux', 'hfgeoubed',
    #                            'upward_geothermal_heat_flux_in_land_ice',
    #                            'W m-2', 'Geothermal heat flux',
    #                            file_remapped_mali_state, ismip7_grid_file,
    #                            exp, output_path)

    # ----------- xvelsurf ------------------
    print("* Writing xvelsurf")
    write_netcdf_2d_state_vars('uReconstructX_sfc', 'xvelsurf',
                               'land_ice_surface_x_velocity',
                               'm s-1', 'Surface velocity in x',
                               file_remapped_mali_state,
                               ismip7_grid_file, output_path, metadata,
                               calendar=calendar)

    # -----------yxvelsurf ------------------
    print("* Writing yvelsurf")
    write_netcdf_2d_state_vars('uReconstructY_sfc', 'yvelsurf',
                               'land_ice_surface_y_velocity',
                               'm s-1', 'Surface velocity in x',
                               file_remapped_mali_state,
                               ismip7_grid_file, output_path, metadata,
                               calendar=calendar)

    # ----------- xvelbase ------------------
    print
    write_netcdf_2d_state_vars('uReconstructX_base', 'xvelbase',
                               'land_ice_basal_x_velocity',
                               'm s-1', 'Basal velocity in x',
                               file_remapped_mali_state,
                               ismip7_grid_file, output_path, metadata,
                               calendar=calendar)

    # ----------- yvelbase ------------------
    print("* Writing yvelbase")
    write_netcdf_2d_state_vars('uReconstructY_base', 'yvelbase',
                               'land_ice_basal_y_velocity',
                               'm s-1', 'Basal velocity in y',
                               file_remapped_mali_state,
                               ismip7_grid_file, output_path, metadata,
                               calendar=calendar)

    # ----------- zvelsurf & zvelbase ------------------
    # ISMIP7 requires these variables, but MALI does not output them.
    # So, we are not processing/writing these variables out.

    # ----------- xvelmean ------------------
    print("* Writing xvelmean")
    write_netcdf_2d_state_vars('xvelmean', 'xvelmean',
                               'land_ice_vertical_mean_x_velocity',
                               'm s-1', 'Mean velocity in x',
                               file_remapped_mali_state,
                               ismip7_grid_file, output_path, metadata,
                               calendar=calendar)

    # ----------- yvelmean ------------------
    print("* Writing yvelmean")
    write_netcdf_2d_state_vars('yvelmean', 'yvelmean',
                               'land_ice_vertical_mean_y_velocity',
                               'm s-1', 'Mean velocity in y',
                               file_remapped_mali_state,
                               ismip7_grid_file, output_path, metadata,
                               calendar=calendar)

    # ----------- litemptop ------------------
    print("* Writing litemptop")
    write_netcdf_2d_state_vars('surfaceTemperature', 'litemptop',
                               'temperature_at_top_of_ice_sheet_model', 'K',
                               'Surface temperature',
                               file_remapped_mali_state,
                               ismip7_grid_file, output_path, metadata,
                               calendar=calendar)

    # ----------- litempbotgr ------------------
    print("* Writing litempbotgr")
    write_netcdf_2d_state_vars('litempbotgr', 'litempbotgr',
                               'temperature_at_base_of_ice_sheet_model', 'K',
                               'Basal temperature beneath grounded ice sheet',
                               file_remapped_mali_state,
                               ismip7_grid_file, output_path, metadata,
                               calendar=calendar)

    # ----------- litempbotfl ------------------
    print("* Writing litempbotfl")
    write_netcdf_2d_state_vars('litempbotfl', 'litempbotfl',
                               'temperature_at_base_of_ice_sheet_model', 'K',
                               'Basal temperature beneath floating ice shelf',
                               file_remapped_mali_state,
                               ismip7_grid_file, output_path, metadata,
                               calendar=calendar)

    # ----------- strbasemag ------------------
    print("* Writing strbasemag")
    write_netcdf_2d_state_vars('strbasemag', 'strbasemag',
                               'land_ice_basal_drag', 'Pa',
                               'Basal drag',
                               file_remapped_mali_state,
                               ismip7_grid_file, output_path, metadata,
                               calendar=calendar)

    # ----------- sftgif ------------------
    print("* Writing sftgif")
    write_netcdf_2d_state_vars('sftgif', 'sftgif',
                               'land_ice_area_fraction', '1',
                               'Land ice area fraction',
                               file_remapped_mali_state,
                               ismip7_grid_file, output_path, metadata,
                               calendar=calendar)

    # ----------- sftgrf ------------------
    print("* Writing sftgrf")
    write_netcdf_2d_state_vars('sftgrf', 'sftgrf',
                               'grounded_ice_sheet_area_fraction', '1',
                               'Grounded ice sheet area fraction',
                               file_remapped_mali_state,
                               ismip7_grid_file, output_path, metadata,
                               calendar=calendar)

    # ----------- sftflf ------------------
    print("* Writing sftflf")
    write_netcdf_2d_state_vars('sftflf', 'sftflf',
                               'floating_ice_shelf_area_fraction', '1',
                               'Floating ice shelf area fraction',
                               file_remapped_mali_state,
                               ismip7_grid_file, output_path, metadata,
                               calendar=calendar)
