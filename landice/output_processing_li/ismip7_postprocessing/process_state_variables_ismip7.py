"""
This script has functions that are needed to post-process and write state
output variables from ISMIP7 simulations.
The input files (i.e., MALI output files) need to have been
concatenated to have yearly data, which can be done using 'ncrcat' command
before using this script.
"""

from netCDF4 import Dataset
import xarray as xr
import numpy as np
from datetime import date
import shutil
import os, sys


def process_state_vars(inputfile_state, tmp_file):
    """
    inputfile_state: output file copy from MALI simulations
    tmp_file: temporary file name
    inputfile_temperature: output temperature file from MALI simulations
    """

    inputfile_state_vars = xr.open_dataset(inputfile_state, engine="netcdf4", decode_cf=False)
    del inputfile_state_vars.daysSinceStart.attrs['units'] # need this line to prevent xarray from reading daysSinceStart as a timedelta type and corrupting values after about 250 years

    # get the mesh description data
    nCells = inputfile_state_vars.dims['nCells']
    nTime = inputfile_state_vars.dims['Time']
    nLayer = inputfile_state_vars.dims['nVertLevels']
    nInterface = nLayer + 1  # inputfile_state_vars.dims['nVertInterfaces']
    cellMask = inputfile_state_vars['cellMask'][:, :]
    basalTemperature = inputfile_state_vars['basalTemperature'][:, :]
    betaSolve = inputfile_state_vars['betaSolve'][:, :]

    inputfile_state_vars['litempbotfl'] = basalTemperature * (cellMask[:, :] & 4) / 4
    inputfile_state_vars['litempbotgr'] = basalTemperature * (1 - (cellMask[:, :] & 4) / 4)

    uxsurf = inputfile_state_vars['uReconstructX'][:, :, 0]
    uysurf = inputfile_state_vars['uReconstructY'][:, :, 0]
    uxbase = inputfile_state_vars['uReconstructX'][:, :, nInterface - 1]
    uybase = inputfile_state_vars['uReconstructY'][:, :, nInterface - 1]
    inputfile_state_vars['uReconstructX_sfc'] = uxsurf
    inputfile_state_vars['uReconstructY_sfc'] = uysurf
    inputfile_state_vars['uReconstructX_base'] = uxbase
    inputfile_state_vars['uReconstructY_base'] = uybase

    inputfile_state_vars['upperSurface'] = np.maximum(0.0, inputfile_state_vars['upperSurface'])

    inputfile_state_vars['sftflf'] = (cellMask[:, :] & 4) / 4 * (cellMask[:, :] & 2) / 2  # floating and dynamic
    inputfile_state_vars['sftgrf'] = ((cellMask[:, :] * 0 + 1) - (cellMask[:, :] & 4) / 4) * (cellMask[:, :] & 2) / 2  # grounded: not-floating & dynamic
    inputfile_state_vars['sftgif'] = (cellMask[:, :] & 2) / 2  # grounded: dynamic ice
    inputfile_state_vars['strbasemag'] = betaSolve[:, :] * ((uxbase[:, :]) ** 2 + (uybase[:, :]) ** 2) **0.5 \
        * (3600.0 * 24.0 * 365.0) \
        * (cellMask[:, :] * 0 + 1 - (cellMask[:, :] & 4) / 4) * (cellMask[:, :] & 2) / 2

    inputfile_state_vars.to_netcdf(tmp_file)
    inputfile_state_vars.close()


def write_netcdf_2d_state_vars(mali_var_name, ismip7_var_name, var_std_name,
                               var_units, var_varname, remapped_mali_outputfile,
                               ismip7_grid_file, exp, output_path):
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
    simulationStartTime = data.variables['simulationStartTime'][:].tostring().decode('utf-8').strip().strip('\x00')
    simulationStartDate = simulationStartTime.split("_")[0]
    daysSinceStart = data.variables['daysSinceStart'][:]
    var_sftgif = data.variables['sftgif'][:, :, :]
    var_sftgrf = data.variables['sftgrf'][:, :, :]
    var_sftflf = data.variables['sftflf'][:, :, :]
    var_mali = data.variables[mali_var_name][:,:,:]
    var_mali[np.where(abs(var_mali + 1e34) < 1e33)] = np.NAN
    timeSteps, latN, lonN = np.shape(var_mali)

    dataOut = Dataset(f'{output_path}/{ismip7_var_name}_{icesheet}_DOE_MALI_{exp}.nc',
                      'w', format='NETCDF4_CLASSIC')
    dataOut.createDimension('time', timeSteps)
    dataOut.createDimension('x', lonN)
    dataOut.createDimension('y', latN)
    dataValues = dataOut.createVariable(ismip7_var_name, 'd',
                                       ('time', 'y', 'x'), fill_value=np.NAN)
    xValues = dataOut.createVariable('x', 'd', ('x'))
    yValues = dataOut.createVariable('y', 'd', ('y'))
    timeValues = dataOut.createVariable('time', 'd', ('time'))
    timeValues[:] = daysSinceStart
    AUTHOR_STR = 'Matthew Hoffman, Trevor Hillebrand, Holly Kyeore Han'
    DATE_STR = date.today().strftime("%d-%b-%Y")

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
            tmp[mask == 0] = np.NAN
            dataValues[i, :, :] = tmp

    for i in range(latN):
        xValues[i] = var_x[i]

    for i in range(lonN):
        yValues[i] = var_y[i]

    dataValues.standard_name = var_std_name
    dataValues.units = var_units
    timeValues.units = f'days since {simulationStartDate}'
    timeValues.calendar = 'noleap'
    timeValues.standard_name = 'time'
    timeValues.long_name = 'time'
    xValues.units = 'm'
    xValues.standard_name = 'x'
    xValues.long_name = 'x'
    yValues.units = 'm'
    yValues.standard_name = 'y'
    yValues.long_name = 'y'
    dataOut.AUTHORS = AUTHOR_STR
    dataOut.MODEL = 'MALI (MPAS-Albany Land Ice model)'
    dataOut.GROUP = 'Los Alamos National Laboratory, Department of Energy'
    dataOut.VARIABLE = var_varname
    dataOut.DATE = DATE_STR
    dataOut.close()


def generate_output_2d_state_vars(file_remapped_mali_state,
                                  ismip7_grid_file, exp, output_path,
                                  icesheet):
    """
    file_remapped_mali_state: output files on mali mesh remapped
    on the ismip7 grid
    ismip7_grid_file: ismip7 original file
    exp: ISMIP7 experiment name
    output_path: path to which the final output files are saved
    """


    # ----------- lithk ------------------
    write_netcdf_2d_state_vars('thickness','lithk', 'land_ice_thickness',
                               'm', 'Ice thickness',
                               file_remapped_mali_state,
                               ismip7_grid_file, exp, output_path)

    # ----------- orog ------------------
    write_netcdf_2d_state_vars('upperSurface','orog', 'surface_altitude', 'm',
                               'Surface elevation',
                               file_remapped_mali_state,
                               ismip7_grid_file,exp, output_path)

    # ----------- base ------------------
    write_netcdf_2d_state_vars('lowerSurface','base', 'base_altitude', 'm',
                               'Base elevation',
                               file_remapped_mali_state,
                               ismip7_grid_file, exp, output_path)

    # ----------- topg ------------------
    write_netcdf_2d_state_vars('bedTopography','topg', 'bedrock_altitude', 'm',
                               'Bedrock elevation',
                               file_remapped_mali_state,
                               ismip7_grid_file, exp, output_path)

    # ----------- hfgeoubed------------------
    # Note: even though this is a flux variable, we are taking a snapshot of it
    # as it does not change with time #Uncomment the function all once basalHeatFlux is outputted in the output stream
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
                               ismip7_grid_file, exp, output_path)

    # -----------yxvelsurf ------------------
    write_netcdf_2d_state_vars('uReconstructY_sfc', 'yvelsurf',
                               'land_ice_surface_y_velocity',
                               'm s-1', 'Surface velocity in x',
                               file_remapped_mali_state,
                               ismip7_grid_file, exp, output_path)

    # ----------- xvelbase ------------------
    write_netcdf_2d_state_vars('uReconstructX_base', 'xvelbase',
                               'land_ice_basal_x_velocity',
                               'm s-1', 'Basal velocity in x',
                               file_remapped_mali_state,
                               ismip7_grid_file, exp, output_path)

    # ----------- yvelbase ------------------
    write_netcdf_2d_state_vars('uReconstructY_base', 'yvelbase',
                               'land_ice_basal_y_velocity',
                               'm s-1', 'Basal velocity in y',
                               file_remapped_mali_state,
                               ismip7_grid_file, exp, output_path)

    # ----------- zvelsurf & zvelbase ------------------
        # ISMIP7 requires these variables, but MALI does not output them.
        # So, we are not processing/writing these variables out.

    # ----------- xvelmean ------------------
    write_netcdf_2d_state_vars('xvelmean', 'xvelmean',
                               'land_ice_vertical_mean_x_velocity',
                               'm s-1', 'Mean velocity in x',
                               file_remapped_mali_state,
                               ismip7_grid_file, exp, output_path)

    # ----------- yvelmean ------------------
    write_netcdf_2d_state_vars('yvelmean', 'yvelmean',
                               'land_ice_vertical_mean_y_velocity',
                               'm s-1', 'Mean velocity in y',
                               file_remapped_mali_state,
                               ismip7_grid_file, exp, output_path)

    # ----------- litemptop ------------------
    write_netcdf_2d_state_vars('surfaceTemperature', 'litemptop',
                               'temperature_at_top_of_ice_sheet_model', 'K',
                               'Surface temperature',
                               file_remapped_mali_state,
                               ismip7_grid_file, exp, output_path)

    # ----------- litempbotgr ------------------
    write_netcdf_2d_state_vars('litempbotgr', 'litempbotgr',
                               'temperature_at_base_of_ice_sheet_model', 'K',
                               'Basal temperature beneath grounded ice sheet',
                               file_remapped_mali_state,
                               ismip7_grid_file,exp, output_path)

    # ----------- litempbotfl ------------------
    write_netcdf_2d_state_vars('litempbotfl', 'litempbotfl',
                               'temperature_at_base_of_ice_sheet_model', 'K',
                               'Basal temperature beneath floating ice shelf',
                               file_remapped_mali_state,
                               ismip7_grid_file, exp, output_path)

    # ----------- strbasemag ------------------
    write_netcdf_2d_state_vars('strbasemag', 'strbasemag',
                               'land_ice_basal_drag ', 'Pa',
                               'Basal drag',
                               file_remapped_mali_state,
                               ismip7_grid_file, exp, output_path)

    # ----------- sftgif ------------------
    write_netcdf_2d_state_vars('sftgif','sftgif',
                               'land_ice_area_fraction', '1',
                               'Land ice area fraction',
                               file_remapped_mali_state,
                               ismip7_grid_file, exp, output_path)

    # ----------- sftgrf ------------------
    write_netcdf_2d_state_vars('sftgrf', 'sftgrf',
                               'grounded_ice_sheet_area_fraction', '1',
                               'Grounded ice sheet area fraction',
                               file_remapped_mali_state,
                               ismip7_grid_file, exp, output_path)

    # ----------- sftflf ------------------
    write_netcdf_2d_state_vars('sftflf','sftflf',
                               'floating_ice_shelf_area_fraction', '1',
                               'Floating ice shelf area fraction',
                               file_remapped_mali_state,
                               ismip7_grid_file, exp, output_path)
