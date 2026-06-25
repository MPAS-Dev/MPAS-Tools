"""
This script has functions that are needed to post-process and write flux
output variables from ISMIP7 simulations.
"""

from netCDF4 import Dataset
import xarray as xr
import numpy as np
from datetime import date
from subprocess import check_call
import os, sys
import warnings


def write_netcdf_2d_flux_vars(mali_var_name, ismip7_var_name, var_std_name,
                              var_units, var_varname, remapped_mali_flux_file,
                              ismip7_grid_file, exp, output_path):

    """
    mali_var_name: variable name on MALI side
    ismip7_var_name: variable name required by ISMIP7
    var_std_name: standard variable name
    var_units: variable units
    var_varname: variable variable name
    remapped_mali_flux_file: mali flux file remapped on the ISMIP7 grid
    ismip7_grid_file: original ISMIP7 file
    exp: experiment name
    output_path: output path to which the output files will be saved
    """

    data_ismip7 = Dataset(ismip7_grid_file, 'r')
    var_x = data_ismip7.variables['x'][:]
    var_y = data_ismip7.variables['y'][:]

    data = Dataset(remapped_mali_flux_file, 'r')
    data.set_auto_mask(False)
    iceMask = data.variables['iceMask'][:, :, :]
    simulationStartTime = data.variables['simulationStartTime'][:].tostring().decode('utf-8').strip().strip('\x00')
    simulationStartDate = simulationStartTime.split("_")[0]
    timeBndsMin = data.variables['timeBndsMin'][:]
    timeBndsMax = data.variables['timeBndsMax'][:]
    if not mali_var_name in data.variables:
        print(f"WARNING: {mali_var_name} not present.  Skipping.")
        data.close()
        return
    var_mali = data.variables[mali_var_name][:,:,:]
    var_mali[np.where(abs(var_mali + 1e34) < 1e33)] = np.NAN
    timeSteps, latN, lonN = np.shape(var_mali)

    dataOut = Dataset(f'{output_path}/{ismip7_var_name}_AIS_DOE_MALI_{exp}.nc',
                      'w', format='NETCDF4_CLASSIC')
    dataOut.createDimension('time', timeSteps)
    dataOut.createDimension('bnds', 2)
    timebndsValues = dataOut.createVariable('time_bnds', 'd', ('time', 'bnds'))
    dataOut.createDimension('x', lonN)
    dataOut.createDimension('y', latN)
    dataValues = dataOut.createVariable(ismip7_var_name, 'd',
                                       ('time', 'y', 'x'), fill_value=np.NAN)
    xValues = dataOut.createVariable('x', 'd', ('x'))
    yValues = dataOut.createVariable('y', 'd', ('y'))
    timeValues = dataOut.createVariable('time', 'd', ('time'))

    AUTHOR_STR = 'Matthew Hoffman, Trevor Hillebrand, Holly Kyeore Han'
    DATE_STR = date.today().strftime("%d-%b-%Y")

    for i in range(timeSteps):
        mask = iceMask[i, :, :]
        tmp = var_mali[i, :, :]
        tmp[mask == 0] = np.NAN
        dataValues[i, :, :] = tmp
        timeValues[i] = (timeBndsMin[i] + timeBndsMax[i]) / 2.0
        timebndsValues[i, 0] = timeBndsMin[i]
        timebndsValues[i, 1] = timeBndsMax[i]

    for i in range(latN):
        xValues[i] = var_x[i]

    for i in range(lonN):
        yValues[i] = var_y[i]

    dataValues.standard_name = var_std_name
    dataValues.units = var_units
    timeValues.bounds = 'time_bnds'
    timeValues.units = f'days since {simulationStartDate}'
    timeValues.calendar = 'noleap'
    timeValues.standard_name = 'time'
    timeValues.long_name = 'time'
    timebndsValues.units = f'days since {simulationStartDate}'
    timebndsValues.calendar = 'noleap'
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
    data.close()


def generate_output_2d_flux_vars(file_remapped_mali_flux,
                                 ismip7_grid_file, exp, output_path):
    """
    file_remapped_mali_flux: flux output file on mali mesh remapped
    onto the ismip7 grid
    ismip7 grid
    ismip7_grid_file: ismip7 original file
    exp: ISMIP7 experiment name
    output_path: path to which the final output files are saved
    """

    print("Writing 2d flux variables")
    # ----------- acabf ------------------
    write_netcdf_2d_flux_vars('avgSMBFlux', 'acabf',
                              'land_ice_surface_specific_mass_balance_flux',
                              'kg m-2 s-1', 'Surface mass balance flux',
                              file_remapped_mali_flux,
                              ismip7_grid_file, exp, output_path)

    # ----------- libmassbffl ------------------
    write_netcdf_2d_flux_vars('avgFloatingBMBFlux', 'libmassbffl',
                              'land_ice_basal_specific_mass_balance_flux',
                              'kg m-2 s-1',
                              'Basal mass balance flux beneath floating ice',
                              file_remapped_mali_flux,
                              ismip7_grid_file, exp, output_path)

    # ----------- libmassbfgr ------------------
    write_netcdf_2d_flux_vars('avgGroundedBMBFlux', 'libmassbfgr',
                              'land_ice_basal_specific_mass_balance_flux',
                              'kg m-2 s-1',
                              'Basal mass balance flux beneath grounded ice',
                              file_remapped_mali_flux,
                              ismip7_grid_file, exp, output_path)

    # ----------- dlithkdt ------------------
    write_netcdf_2d_flux_vars('avgDhdt', 'dlithkdt',
                              'tendency_of_land_ice_thickness',
                              'm s-1',
                              'Ice thickness imbalance',
                              file_remapped_mali_flux,
                              ismip7_grid_file, exp, output_path)

    # ----------- licalvf ------------------
    write_netcdf_2d_flux_vars('avgCalvingFlux', 'licalvf',
                              'land_ice_specific_mass_flux_due_to_calving',
                              'kg m-2 s-1',
                              'Calving flux',
                              file_remapped_mali_flux,
                              ismip7_grid_file, exp, output_path)

    # ----------- lifmassbf ------------------
    # Note: facemelting and calving flux are combined above
    write_netcdf_2d_flux_vars('avgFaceMeltFlux', 'lifmassbf',
                              'TBD by ISMIP7',
                              'kg m-2 s-1',
                              'Ice front melt flux',
                              file_remapped_mali_flux,
                              ismip7_grid_file, exp, output_path)

    # ----------- ligroundf ------------------
    write_netcdf_2d_flux_vars('avgGroundingLineFlux', 'ligroundf',
                              'TBD by ISMIP7',
                              'kg m-2 s-1',
                              'Grounding line flux',
                              file_remapped_mali_flux,
                              ismip7_grid_file, exp, output_path)
