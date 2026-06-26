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
from validate import validate_mali_files


EXPECTED_FLUX_VARIABLES = [
    'daysSinceStart', 'simulationStartTime',
    'timeBndsMin', 'timeBndsMax', 'iceMask',
    'avgSMBFlux', 'avgFloatingBMBFlux', 'avgGroundedBMBFlux',
    'avgDhdt', 'avgCalvingFlux', 'avgGroundingLineFlux',
]


def check_flux_files(files):
    """
    Validate a list of MALI flux output files before processing.

    See validate.validate_mali_files for full details of checks performed.
    """
    validate_mali_files(files, EXPECTED_FLUX_VARIABLES, label='flux')


def process_flux_vars(files, tmp_file):
    """
    Concatenate/prepare flux files into a temporary file for remapping.

    Parameters
    ----------
    files : list of str
        Sorted list of MALI flux output file paths.
    tmp_file : str
        Temporary output file name.
    """
    ds_flux = xr.open_mfdataset(files, combine='nested',
                                concat_dim='Time',
                                engine='netcdf4',
                                decode_cf=False,
                                data_vars='minimal',
                                coords='minimal',
                                compat='override')
    if 'daysSinceStart' in ds_flux and 'units' in ds_flux['daysSinceStart'].attrs:
        del ds_flux['daysSinceStart'].attrs['units']
    ds_flux.to_netcdf(tmp_file)
    ds_flux.close()


def write_netcdf_2d_flux_vars(mali_var_name, ismip7_var_name, var_std_name,
                              var_units, var_varname, remapped_mali_flux_file,
                              ismip7_grid_file, output_path, metadata):

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

    dataOut = Dataset(f'{output_path}/{ismip7_var_name}_{metadata["icesheet"]}_{metadata["group_nickname"]}_MALI_{metadata["exp"]}.nc',
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
    DATE_STR = metadata['date']

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
    dataOut.AUTHORS = metadata['authors']
    dataOut.MODEL = metadata['model']
    dataOut.GROUP = metadata['group']
    dataOut.VARIABLE = var_varname
    dataOut.DATE = metadata['date']
    dataOut.close()
    data.close()


def generate_output_2d_flux_vars(file_remapped_mali_flux,
                                 ismip7_grid_file, output_path, metadata):
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
                              ismip7_grid_file, output_path, metadata)

    # ----------- libmassbffl ------------------
    write_netcdf_2d_flux_vars('avgFloatingBMBFlux', 'libmassbffl',
                              'land_ice_basal_specific_mass_balance_flux',
                              'kg m-2 s-1',
                              'Basal mass balance flux beneath floating ice',
                              file_remapped_mali_flux,
                              ismip7_grid_file, output_path, metadata)

    # ----------- libmassbfgr ------------------
    write_netcdf_2d_flux_vars('avgGroundedBMBFlux', 'libmassbfgr',
                              'land_ice_basal_specific_mass_balance_flux',
                              'kg m-2 s-1',
                              'Basal mass balance flux beneath grounded ice',
                              file_remapped_mali_flux,
                              ismip7_grid_file, output_path, metadata)

    # ----------- dlithkdt ------------------
    write_netcdf_2d_flux_vars('avgDhdt', 'dlithkdt',
                              'tendency_of_land_ice_thickness',
                              'm s-1',
                              'Ice thickness imbalance',
                              file_remapped_mali_flux,
                              ismip7_grid_file, output_path, metadata)

    # ----------- licalvf ------------------
    write_netcdf_2d_flux_vars('avgCalvingFlux', 'licalvf',
                              'land_ice_specific_mass_flux_due_to_calving',
                              'kg m-2 s-1',
                              'Calving flux',
                              file_remapped_mali_flux,
                              ismip7_grid_file, output_path, metadata)

    # ----------- lifmassbf ------------------
    # Note: facemelting and calving flux are combined above
    write_netcdf_2d_flux_vars('avgFaceMeltFlux', 'lifmassbf',
                              'TBD by ISMIP7',
                              'kg m-2 s-1',
                              'Ice front melt flux',
                              file_remapped_mali_flux,
                              ismip7_grid_file, output_path, metadata)

    # ----------- ligroundf ------------------
    write_netcdf_2d_flux_vars('avgGroundingLineFlux', 'ligroundf',
                              'TBD by ISMIP7',
                              'kg m-2 s-1',
                              'Grounding line flux',
                              file_remapped_mali_flux,
                              ismip7_grid_file, output_path, metadata)


def process_flux_pipeline(flux_files, mapping_file, ismip7_grid_file,
                          output_path, metadata):
    """
    Full flux-variable processing pipeline: validate, concatenate, remap, write.

    Parameters
    ----------
    flux_files : list of str
        Sorted list of MALI flux output file paths.
    mapping_file : str
        Path to the ESMF mapping/weights file.
    ismip7_grid_file : str
        Path to the ISMIP7 grid file.
    output_path : str
        Directory for output files.
    metadata : dict
        Submission metadata dict.
    """
    check_flux_files(flux_files)

    print("Preparing concatenated flux file for remapping.")
    tmp_flux_file = 'tmp_flux.nc'
    process_flux_vars(flux_files, tmp_flux_file)

    print("Remapping flux file.")
    processed_file_flux = 'processed_flux.nc'
    check_call(["ncremap",
                "-i", tmp_flux_file,
                "-o", processed_file_flux,
                "-m", mapping_file,
                "-P", "mpas"])

    print("Writing processed and remapped flux fields to ISMIP7 file format.")
    generate_output_2d_flux_vars(processed_file_flux,
                                 ismip7_grid_file,
                                 output_path, metadata)

    os.remove(tmp_flux_file)
    os.remove(processed_file_flux)
