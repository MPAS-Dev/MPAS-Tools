"""
This script has functions that are needed to post-process and write state
output variables from ISMIP7 simulations.
"""

from netCDF4 import Dataset, num2date, date2num
import xarray as xr
import numpy as np
import os
from subprocess import check_call
from validate import validate_mali_files
from output_naming import build_output_filename


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
    if inputfile_state_vars.sizes['Time'] > 0 and np.isclose(days_first, 0.0):
        print("Skipping state data at simulation start time.")
        inputfile_state_vars = inputfile_state_vars.isel(Time=slice(1, None))

    # get the mesh description data
    nLayer = inputfile_state_vars.sizes['nVertLevels']
    nInterface = nLayer + 1  # inputfile_state_vars.sizes['nVertInterfaces']
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

    missing = [
        var for var in REQUIRED_STATE_REMAPPING_VARIABLES
        if var not in inputfile_state_vars
    ]
    if missing:
        raise ValueError(
            "Processed state dataset is missing required output variables: "
            f"{missing}"
        )

    output_state_vars = inputfile_state_vars[REQUIRED_STATE_REMAPPING_VARIABLES]
    output_state_vars.to_netcdf(tmp_file)
    output_state_vars.close()
    inputfile_state_vars.close()


def write_netcdf_2d_state_vars(
        mali_var_name,
        ismip7_var_name,
        var_std_name,
        var_units,
        var_varname,
        remapped_mali_outputfile,
        ismip7_grid_file,
        output_path,
        metadata):
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
    """

    data_ismip7 = Dataset(ismip7_grid_file, 'r')
    var_x = data_ismip7.variables['x'][:]
    var_y = data_ismip7.variables['y'][:]

    data = Dataset(remapped_mali_outputfile, 'r')
    data.set_auto_mask(False)
    simulationStartTime = data.variables['simulationStartTime'][:].tobytes(
    ).decode('utf-8').strip().strip('\x00')
    simulationStartDate = simulationStartTime.split("_")[0]
    # Convert from MALI noleap time axis to Gregorian days since 1850-01-01
    daysSinceStart = data.variables['daysSinceStart'][:]
    daysSinceStart = date2num(
        num2date(
            daysSinceStart,
            units=f'days since {simulationStartDate}',
            calendar='noleap',
        ),
        units='days since 1850-01-01',
        calendar='standard',
    )
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
    dataValues = dataOut.createVariable(ismip7_var_name, 'd',
                                        ('time', 'y', 'x'), fill_value=np.nan)
    xValues = dataOut.createVariable('x', 'd', ('x'))
    yValues = dataOut.createVariable('y', 'd', ('y'))
    timeValues = dataOut.createVariable('time', 'd', ('time'))
    timeValues[:] = daysSinceStart
    for i in range(timeSteps):
        if ismip7_var_name == 'sftgif':
            dataValues[i, :, :] = var_mali[i, :, :]
        else:
            if ismip7_var_name == 'litempbotgr':
                mask = var_sftgrf[i, :, :]
            elif ismip7_var_name == 'litempbotfl':
                mask = var_sftflf[i, :, :]
            elif ismip7_var_name == 'topg':
                mask = np.ones(var_mali.shape[1:])  # don't mask topg
            else:
                mask = var_sftgif[i, :, :]
            tmp = var_mali[i, :, :]
            tmp[mask == 0] = np.nan
            dataValues[i, :, :] = tmp

    for i in range(latN):
        xValues[i] = var_x[i]

    for i in range(lonN):
        yValues[i] = var_y[i]

    dataValues.standard_name = var_std_name
    dataValues.units = var_units
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
                           output_path, metadata):
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
    """
    metadata = metadata.copy()

    check_state_files(state_files)

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
                                  ismip7_grid_file, output_path, metadata)

    os.remove(tmp_file)
    os.remove(processed_and_remapped_file_state)


def generate_output_2d_state_vars(file_remapped_mali_state,
                                  ismip7_grid_file, output_path, metadata):
    """
    file_remapped_mali_state: output files on mali mesh remapped
    on the ismip7 grid
    ismip7_grid_file: ismip7 original file
    exp: ISMIP7 experiment name
    output_path: path to which the final output files are saved
    """

    # ----------- lithk ------------------
    write_netcdf_2d_state_vars('thickness', 'lithk', 'land_ice_thickness',
                               'm', 'Ice thickness',
                               file_remapped_mali_state,
                               ismip7_grid_file, output_path, metadata)

    # ----------- orog ------------------
    write_netcdf_2d_state_vars('upperSurface', 'orog', 'surface_altitude', 'm',
                               'Surface elevation',
                               file_remapped_mali_state,
                               ismip7_grid_file, output_path, metadata)

    # ----------- base ------------------
    write_netcdf_2d_state_vars('lowerSurface', 'base', 'base_altitude', 'm',
                               'Base elevation',
                               file_remapped_mali_state,
                               ismip7_grid_file, output_path, metadata)

    # ----------- topg ------------------
    write_netcdf_2d_state_vars(
        'bedTopography',
        'topg',
        'bedrock_altitude',
        'm',
        'Bedrock elevation',
        file_remapped_mali_state,
        ismip7_grid_file,
        output_path,
        metadata)

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
    write_netcdf_2d_state_vars('uReconstructX_sfc', 'xvelsurf',
                               'land_ice_surface_x_velocity',
                               'm s-1', 'Surface velocity in x',
                               file_remapped_mali_state,
                               ismip7_grid_file, output_path, metadata)

    # -----------yxvelsurf ------------------
    write_netcdf_2d_state_vars('uReconstructY_sfc', 'yvelsurf',
                               'land_ice_surface_y_velocity',
                               'm s-1', 'Surface velocity in x',
                               file_remapped_mali_state,
                               ismip7_grid_file, output_path, metadata)

    # ----------- xvelbase ------------------
    write_netcdf_2d_state_vars('uReconstructX_base', 'xvelbase',
                               'land_ice_basal_x_velocity',
                               'm s-1', 'Basal velocity in x',
                               file_remapped_mali_state,
                               ismip7_grid_file, output_path, metadata)

    # ----------- yvelbase ------------------
    write_netcdf_2d_state_vars('uReconstructY_base', 'yvelbase',
                               'land_ice_basal_y_velocity',
                               'm s-1', 'Basal velocity in y',
                               file_remapped_mali_state,
                               ismip7_grid_file, output_path, metadata)

    # ----------- zvelsurf & zvelbase ------------------
    # ISMIP7 requires these variables, but MALI does not output them.
    # So, we are not processing/writing these variables out.

    # ----------- xvelmean ------------------
    write_netcdf_2d_state_vars('xvelmean', 'xvelmean',
                               'land_ice_vertical_mean_x_velocity',
                               'm s-1', 'Mean velocity in x',
                               file_remapped_mali_state,
                               ismip7_grid_file, output_path, metadata)

    # ----------- yvelmean ------------------
    write_netcdf_2d_state_vars('yvelmean', 'yvelmean',
                               'land_ice_vertical_mean_y_velocity',
                               'm s-1', 'Mean velocity in y',
                               file_remapped_mali_state,
                               ismip7_grid_file, output_path, metadata)

    # ----------- litemptop ------------------
    write_netcdf_2d_state_vars('surfaceTemperature', 'litemptop',
                               'temperature_at_top_of_ice_sheet_model', 'K',
                               'Surface temperature',
                               file_remapped_mali_state,
                               ismip7_grid_file, output_path, metadata)

    # ----------- litempbotgr ------------------
    write_netcdf_2d_state_vars('litempbotgr', 'litempbotgr',
                               'temperature_at_base_of_ice_sheet_model', 'K',
                               'Basal temperature beneath grounded ice sheet',
                               file_remapped_mali_state,
                               ismip7_grid_file, output_path, metadata)

    # ----------- litempbotfl ------------------
    write_netcdf_2d_state_vars('litempbotfl', 'litempbotfl',
                               'temperature_at_base_of_ice_sheet_model', 'K',
                               'Basal temperature beneath floating ice shelf',
                               file_remapped_mali_state,
                               ismip7_grid_file, output_path, metadata)

    # ----------- strbasemag ------------------
    write_netcdf_2d_state_vars('strbasemag', 'strbasemag',
                               'land_ice_basal_drag ', 'Pa',
                               'Basal drag',
                               file_remapped_mali_state,
                               ismip7_grid_file, output_path, metadata)

    # ----------- sftgif ------------------
    write_netcdf_2d_state_vars('sftgif', 'sftgif',
                               'land_ice_area_fraction', '1',
                               'Land ice area fraction',
                               file_remapped_mali_state,
                               ismip7_grid_file, output_path, metadata)

    # ----------- sftgrf ------------------
    write_netcdf_2d_state_vars('sftgrf', 'sftgrf',
                               'grounded_ice_sheet_area_fraction', '1',
                               'Grounded ice sheet area fraction',
                               file_remapped_mali_state,
                               ismip7_grid_file, output_path, metadata)

    # ----------- sftflf ------------------
    write_netcdf_2d_state_vars('sftflf', 'sftflf',
                               'floating_ice_shelf_area_fraction', '1',
                               'Floating ice shelf area fraction',
                               file_remapped_mali_state,
                               ismip7_grid_file, output_path, metadata)
