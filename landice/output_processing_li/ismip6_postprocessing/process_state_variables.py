"""
# this script has functions that are needed to post-process and write state
# output variables from ISMIP6 simulations.
# The input files (i.e., MALI output files) need to have been
# concatenated to have yearly data, which can be done using 'ncrcat' command
# before using this script.
"""

from netCDF4 import Dataset
import xarray as xr
import numpy as np
import shutil
import os


def process_state_vars(inputfile_state, tmp_file, inputfile_temperature=None):
    """
    inputfile_state: output file copy from MALI simulations
    tmp_file: temporary file name
    inputfile_temperature: output temperature file from MALI simulations
    """

    if inputfile_temperature is None:
        print("temperature file is not provided. Reading the output all file")
        input_temp = xr.open_dataset(inputfile_state, engine="netcdf4")
    else:
        input_temp = xr.open_dataset(inputfile_temperature,
                                     engine="netcdf4")  # this might not be needed if basal and surface temps are already recorded in the output_all.nc file

    inputfile_state_vars = xr.open_dataset(inputfile_state, engine="netcdf4")

    # get the mesh description data
    nCells = inputfile_state_vars.dims['nCells']
    nTime = inputfile_state_vars.dims['Time']
    nLayer = inputfile_state_vars.dims['nVertLevels']
    nInterface = nLayer + 1  # inputfile_state_vars.dims['nVertInterfaces']
    cellMask = inputfile_state_vars['cellMask'][:, :]
    basalTemperature = input_temp['basalTemperature'][:, :]
    # surfaceTemperature = input_temp['surfaceTemperature'][:,:] #this doesnt seem necessary
    # basalSpeed = inputfile_state_vars['basalSpeed'][:,:] # check: this if this field will be outputted in the final version.
    betaSolve = inputfile_state_vars['betaSolve'][:, :]  # check: this field only exists in the initial file

    inputfile_state_vars['litempbotfl'] = basalTemperature * (cellMask[:, :] & 4) / 4
    inputfile_state_vars['litempbotgr'] = basalTemperature * (1 - (cellMask[:, :] & 4) / 4)

    uxbase = inputfile_state_vars['uReconstructX'][:, :, nInterface - 1]
    uybase = inputfile_state_vars['uReconstructY'][:, :, nInterface - 1]
    uxsurf = inputfile_state_vars['uReconstructX'][:, :, 0]
    uysurf = inputfile_state_vars['uReconstructY'][:, :, 0]
    xvelmean = inputfile_state_vars['xvelmean'][:, :]
    yvelmean = inputfile_state_vars['yvelmean'][:, :]
    interfaceLayer = np.zeros(nInterface)
    thkLayer = inputfile_state_vars['layerThicknessFractions'][:]

    interfaceLayer[0] = 0.5 * thkLayer[0]
    for i in range(1, nLayer - 1):
        interfaceLayer[i] = 0.5 * (thkLayer[i] + thkLayer[i - 1])
    interfaceLayer[nLayer - 1] = 0.5 * thkLayer[nLayer - 1]

    tmp_ux_mean = inputfile_state_vars['thickness'].copy()
    tmp_uy_mean = inputfile_state_vars['thickness'].copy()
    for i in range(nTime):
        tmp_ux_mean[i, :] = np.sum(inputfile_state_vars['uReconstructX'][i, :, :] * interfaceLayer, axis=1)
        tmp_uy_mean[i, :] = np.sum(inputfile_state_vars['uReconstructY'][i, :, :] * interfaceLayer, axis=1)

    inputfile_state_vars['uReconstructX_mean'] = xvelmean
    inputfile_state_vars['uReconstructY_mean'] = yvelmean
    inputfile_state_vars['uReconstructX_sfc'] = uxsurf
    inputfile_state_vars['uReconstructY_sfc'] = uysurf
    inputfile_state_vars['uReconstructX_base'] = uxbase
    inputfile_state_vars['uReconstructY_base'] = uybase

    inputfile_state_vars['sftflf'] = (cellMask[:, :] & 4) / 4 * (cellMask[:, :] & 2) / 2
    inputfile_state_vars['sftgrf'] = ((cellMask[:, :] * 0 + 1) - (cellMask[:, :] & 4) / 4) * (cellMask[:, :] & 2) / 2
    inputfile_state_vars['sftgif'] = (cellMask[:, :] & 2) / 2
    # inputfile_state_vars['strbasemag'] = (betaSolve[:, :] * (basalSpeed) * 365.0 * 24.0 * 3600.0)  # HH: comment out for now until we know what it is. * (cellMask[:,:]*0+1-(cellMask[:,:]&4)/4) * (cellMask[:,:]&2)/2
    inputfile_state_vars['strbasemag'] = (betaSolve[:, :] * ((uxbase[:, :]) ** 2 + (uybase[:, :]) ** 2 + 1.0e-10) ** (
        0.5) * 365.0 * 24.0 * 3600.0) * (cellMask[:, :] * 0 + 1 - (cellMask[:, :] & 4) / 4) * (cellMask[:, :] & 2) / 2

    inputfile_state_vars.to_netcdf(tmp_file)
    inputfile_state_vars.close()
    input_temp.close()


def write_netcdf_2d_state_vars(mali_var_name, ismip6_var_name, var_std_name,
                               var_units, var_varname,
                               remapped_mali_outputfile, ismip6_grid_file,
                               exp, output_path):
    """
    mali_var_name: variable name on MALI side
    ismip6_var_name: variable name required by ISMIP6
    var_std_name: standard variable name
    var_units: variable units
    var_varname: variable variable name
    remapped_mali_outputfile: mali state file remapped on the ISMIP6 grid
    ismip6_grid_file: original ISMIP6 file
    exp: experiment name
    output_path: output path to which the output files will be saved
    """

    spy = 365.0
    timeSpan = 1  # HH:check - is something that's needed?
    data_ismip6 = Dataset(ismip6_grid_file, 'r')
    var_x = data_ismip6.variables['x'][:]
    var_y = data_ismip6.variables['y'][:]

    data = Dataset(remapped_mali_outputfile, 'r')
    data.set_auto_mask(False)
    var_sftgif = data.variables['sftgif'][:, :, :]
    var_mali = data.variables[mali_var_name][:,:,:]
    var_mali[np.where(abs(var_mali + 1e34) < 1e33)] = np.NAN
    timeSteps, latN, lonN = np.shape(var_mali)

    dataOut = Dataset(f'{output_path}/{ismip6_var_name}_AIS_DOE_MALI_{exp}.nc',
                      'w', format='NETCDF4_CLASSIC')
    dataOut.createDimension('time', timeSteps)
    dataOut.createDimension('x', lonN)
    dataOut.createDimension('y', latN)
    dataValues = dataOut.createVariable(ismip6_var_name, 'f4',
                                       ('time', 'y', 'x'), fill_value=np.NAN)
    xValues = dataOut.createVariable('x', 'f4', ('x'))
    yValues = dataOut.createVariable('y', 'f4', ('y'))
    timeValues = dataOut.createVariable('time', 'f4', ('time'))
    AUTHOR_STR = 'Matthew Hoffman, Trevor Hillebrand, Holly Kyeore Han'
    DATE_STR = '23-Oct-2022'

    for i in range(timeSteps):

        timeValues[i] = i * timeSpan * spy

        if not ismip6_var_name == 'sftgif':
            mask = var_sftgif[i, :, :]
            tmp = var_mali[i, :, :]
            tmp[mask == 0] = np.NAN
            dataValues[i, :, :] = tmp
        else:
            dataValues[i, :, :] = var_mali[i, :, :]

    for i in range(latN):
        xValues[i] = var_x[i]

    for i in range(lonN):
        yValues[i] = var_y[i]

    dataValues.standard_name = var_std_name
    dataValues.units = var_units
    timeValues.units = 'days since 2015-01-01'
    timeValues.calendar = '365_day'
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
    dataOut.GROUP = 'Los Alamos National Laboratory'
    dataOut.VARIABLE = var_varname
    dataOut.DATE = DATE_STR
    dataOut.close()


def generate_output_2d_state_vars(file_remapped_mali_state, ismip6_grid_file,
                                  exp, output_path):
    """
    file_remapped_mali_state: output files on mali mesh remapped
    on the ismip6 grid
    ismip6_grid_file: ismip6 original file
    exp: ISMIP6 experiment name
    output_path: path to which the final output files are saved
    """


    # ----------- lithk ------------------
    write_netcdf_2d_state_vars('thickness','lithk', 'land_ice_thickness',
                               'm', 'Ice thickness',
                               file_remapped_mali_state, ismip6_grid_file,
                               exp, output_path)

    # ----------- orog ------------------
    write_netcdf_2d_state_vars('upperSurface','orog', 'surface_altitude', 'm',
                               'Surface elevation',
                               file_remapped_mali_state, ismip6_grid_file,
                               exp, output_path)

    # ----------- base ------------------
    write_netcdf_2d_state_vars('lowerSurface','base', 'base_altitude', 'm',
                               'Base elevation',
                               file_remapped_mali_state, ismip6_grid_file,
                               exp, output_path)

    # ----------- topg ------------------
    write_netcdf_2d_state_vars('bedTopography','topg', 'bedrock_altitude', 'm',
                               'Bedrock elevation',
                               file_remapped_mali_state, ismip6_grid_file,
                               exp, output_path)

    # ----------- hfgeoubed------------------
    # Note: even though this is a flux variable, we are taking a snapshot of it
    # as it does not change with time #Uncomment the function all once basalHeatFlux is outputted in the output stream
    # write_netcdf_2d_state_vars('basalHeatFlux', 'hfgeoubed',
    #                            'upward_geothermal_heat_flux_in_land_ice',
    #                            'W m^-2', 'Geothermal heat flux',
    #                            file_remapped_mali_state, ismip6_grid_file,
    #                            exp, output_path)

    # ----------- xvelsurf ------------------
    write_netcdf_2d_state_vars('uReconstructX_sfc', 'xvelsurf',
                               'land_ice_surface_x_velocity',
                               'm s^-1', 'Surface velocity in x',
                               file_remapped_mali_state, ismip6_grid_file,
                               exp, output_path)

    # -----------yxvelsurf ------------------
    write_netcdf_2d_state_vars('uReconstructY_sfc', 'yvelsurf',
                               'land_ice_surface_y_velocity',
                               'm s^-1', 'Surface velocity in x',
                               file_remapped_mali_state, ismip6_grid_file,
                               exp, output_path)

    # ----------- xvelbase ------------------
    write_netcdf_2d_state_vars('uReconstructX_base', 'xvelbase',
                               'land_ice_basal_x_velocity',
                               'm s^-1', 'Basal velocity in x',
                               file_remapped_mali_state, ismip6_grid_file,
                               exp, output_path)

    # ----------- yvelbase ------------------
    write_netcdf_2d_state_vars('uReconstructY_base', 'yvelbase',
                               'land_ice_basal_y_velocity',
                               'm s^-1', 'Basal velocity in y',
                               file_remapped_mali_state, ismip6_grid_file,
                               exp, output_path)

    # ----------- zvelsurf & zvelbase ------------------
        # ISMIP6 requires these variables, but MALI does not output them.
        # So, we are not processing/writing these variables out.

    # ----------- xvelmean ------------------
    write_netcdf_2d_state_vars('xvelmean', 'xvelmean',
                               'land_ice_vertical_mean_x_velocity',
                               'm s^-1', 'Mean velocity in x',
                               file_remapped_mali_state, ismip6_grid_file,
                               exp, output_path)

    # ----------- yvelmean ------------------
    write_netcdf_2d_state_vars('yvelmean', 'yvelmean',
                               'land_ice_vertical_mean_y_velocity',
                               'm s^-1', 'Mean velocity in y',
                               file_remapped_mali_state, ismip6_grid_file,
                               exp, output_path)

    # ----------- litemptop ------------------
    write_netcdf_2d_state_vars('surfaceTemperature', 'litemptop',
                               'temperature_at_top_of_ice_sheet_model', 'K',
                               'Surface temperature',
                               file_remapped_mali_state, ismip6_grid_file,
                               exp, output_path)

    # ----------- litempbotgr ------------------
    write_netcdf_2d_state_vars('litempbotgr', 'litempbotgr',
                               'temperature_at_base_of_ice_sheet_model', 'K',
                               'Basal temperature beneath grounded ice sheet',
                               file_remapped_mali_state, ismip6_grid_file,
                               exp, output_path)

    # ----------- litempbotfl ------------------
    write_netcdf_2d_state_vars('litempbotfl', 'litempbotfl',
                               'temperature_at_base_of_ice_sheet_model', 'K',
                               'Basal temperature beneath floating ice shelf',
                               file_remapped_mali_state, ismip6_grid_file,
                               exp, output_path)

    # ----------- strbasemag ------------------
    write_netcdf_2d_state_vars('strbasemag', 'strbasemag',
                               'land_ice_basal_drag ', 'Pa',
                               'Basal drag',
                               file_remapped_mali_state, ismip6_grid_file,
                               exp, output_path)

    # ----------- sftgif ------------------
    write_netcdf_2d_state_vars('sftgif','sftgif',
                               'land_ice_area_fraction', '1',
                               'Land ice area fraction',
                               file_remapped_mali_state, ismip6_grid_file,
                               exp, output_path)

    # ----------- sftgrf ------------------
    write_netcdf_2d_state_vars('sftgrf', 'sftgrf',
                               'grounded_ice_sheet_area_fraction', '1',
                               'Grounded ice sheet area fraction',
                               file_remapped_mali_state, ismip6_grid_file,
                               exp, output_path)

    # ----------- sftflf ------------------
    write_netcdf_2d_state_vars('sftflf','sftflf',
                               'floating_ice_shelf_area_fraction', '1',
                               'Floating ice shelf area fraction',
                               file_remapped_mali_state, ismip6_grid_file,
                               exp, output_path)


def generate_output_1d_vars(global_stats_file, exp, output_path=None):
    """
    This code processes both state and flux 1D variables
    global_stats_file: MALI globalStats.nc output file
    exp: ISMIP6 experiment number
    output_path:
    """

    if not os.path.exists(output_path):
        output_path = os.getcwd()

    AUTHOR_STR = 'Matthew Hoffman, Trevor Hillebrand, Holly Kyeore Han'
    DATE_STR = '23-Oct-2022'

    # make a copy of the original globalStats file
    shutil.copy(global_stats_file, "copy_globalStats.nc")

    data = Dataset(global_stats_file, 'r+')
    xtime = data.variables['xtime'][:, :]
    dt = data.variables['deltat'][:]

    timeSteps = len(xtime)  # total length of data
    timeSteps_out = len(np.arange(2015, 2301))  # yearly timestep

    # read in state variables
    vol = data.variables['totalIceVolume'][:]
    vaf = data.variables['volumeAboveFloatation'][:]
    gia = data.variables['groundedIceArea'][:]
    fia = data.variables['floatingIceArea'][:]

    # read in flux variables over which yearly average will be taken
    smb = data.variables['totalSfcMassBal'][:]
    bmb = data.variables['totalBasalMassBal'][:]
    cfx = data.variables['totalCalvingFlux'][:]
    gfx = data.variables['groundingLineFlux'][:]

    # initialize empty arrays for with possible maximum dimensions
    # below variables (X_yearly) will store all available data values (column)
    # for each year (row)
    vol_yearly = np.zeros((timeSteps, timeSteps_out))
    vaf_yearly = np.zeros((timeSteps, timeSteps_out))
    gia_yearly = np.zeros((timeSteps, timeSteps_out))
    fia_yearly = np.zeros((timeSteps, timeSteps_out))
    smb_yearly = np.zeros((timeSteps, timeSteps_out))
    bmb_yearly = np.zeros((timeSteps, timeSteps_out))
    cfx_yearly = np.zeros((timeSteps, timeSteps_out))
    gfx_yearly = np.zeros((timeSteps, timeSteps_out))
    dt_yearly = np.zeros((timeSteps, timeSteps_out))

    # initialize 1D variables that will store data value on the
    # January 1st of each year
    vol_snapshot = np.zeros(timeSteps_out)
    vaf_snapshot = np.zeros(timeSteps_out)
    gia_snapshot = np.zeros(timeSteps_out)
    fia_snapshot = np.zeros(timeSteps_out)
    smb_avg = np.zeros(timeSteps_out)
    bmb_avg = np.zeros(timeSteps_out)
    cfx_avg = np.zeros(timeSteps_out)
    gfx_avg = np.zeros(timeSteps_out)

    # initialize counter variables
    yearIndex = 0
    yearOld = 0
    for i in range(timeSteps):

        year = int(xtime[i][0:4].tobytes())
        if year > yearOld:
            dataIndex = 0
            yearIndex = yearIndex + 1

        vol_yearly[dataIndex, yearIndex - 1] = vol[i]
        vaf_yearly[dataIndex, yearIndex - 1] = vaf[i]
        gia_yearly[dataIndex, yearIndex - 1] = gia[i]
        fia_yearly[dataIndex, yearIndex - 1] = fia[i]
        smb_yearly[dataIndex, yearIndex - 1] = smb[i]
        bmb_yearly[dataIndex, yearIndex - 1] = bmb[i]
        cfx_yearly[dataIndex, yearIndex - 1] = cfx[i]
        gfx_yearly[dataIndex, yearIndex - 1] = gfx[i]

        dt_yearly[dataIndex, yearIndex - 1] = dt[i]

        dataIndex = dataIndex + 1
        yearOld = year

    for i in range(timeSteps_out):
        vol_snapshot[i] = vol_yearly[0, i]
        vaf_snapshot[i] = vaf_yearly[0, i]
        gia_snapshot[i] = gia_yearly[0, i]
        fia_snapshot[i] = fia_yearly[0, i]
        smbi = smb_yearly[:, i]
        bmbi = bmb_yearly[:, i]
        cfxi = cfx_yearly[:, i]
        gfxi = gfx_yearly[:, i]

        # take the average of the flux variables
        smb_avg[i] = np.mean(smbi[np.where(smbi != 0)])
        bmb_avg[i] = np.mean(bmbi[np.where(bmbi != 0)])
        cfx_avg[i] = np.mean(cfxi[np.where(cfxi != 0)])
        gfx_avg[i] = np.mean(gfxi[np.where(gfxi != 0)])

    # -------------- lim ------------------
    data_scalar = Dataset(f'{output_path}/lim_AIS_DOE_MALI_{exp}.nc', 'w', format='NETCDF4_CLASSIC')
    data_scalar.createDimension('time', timeSteps_out)
    limValues = data_scalar.createVariable('lim', 'f4', ('time'))
    timeValues = data_scalar.createVariable('time', 'f4', ('time'))
    for i in range(timeSteps_out):
        limValues[i] = vol_snapshot[i] * 910
        timeValues[i] = i * 365.0
    timeValues.units = 'days since 2015-01-01'
    timeValues.calendar = '365_day'
    timeValues.standard_name = 'time'
    timeValues.long_name = 'time'
    limValues.standard_name = 'land_ice_mass'
    limValues.units = 'kg'
    data_scalar.AUTHORS = AUTHOR_STR
    data_scalar.MODEL = 'MALI (MPAS-Albany Land Ice model)'
    data_scalar.GROUP= 'Los Alamos National Laboratory'
    data_scalar.VARIABLE = 'Total ice mass'
    data_scalar.DATE = DATE_STR
    data_scalar.close()

    # -------------- limnsw ------------------
    data_scalar = Dataset(f'{output_path}/limnsw_AIS_DOE_MALI_{exp}.nc', 'w', format='NETCDF4_CLASSIC')
    data_scalar.createDimension('time', timeSteps_out)
    limnswValues = data_scalar.createVariable('limnsw', 'f4', ('time'))
    timeValues = data_scalar.createVariable('time', 'f4', ('time'))
    for i in range(timeSteps_out):
        limnswValues[i] = vaf_snapshot[i] * 910
        timeValues[i] = i * 365.0
    timeValues.units = 'days since 2015-01-01'
    timeValues.calendar = '365_day'
    timeValues.standard_name = 'time'
    timeValues.long_name = 'time'
    limnswValues.standard_name = 'land_ice_mass_not_displacing_sea_water'
    limnswValues.units = 'kg'
    data_scalar.AUTHORS = AUTHOR_STR
    data_scalar.MODEL = 'MALI (MPAS-Albany Land Ice model)'
    data_scalar.GROUP = 'Los Alamos National Laboratory'
    data_scalar.VARIABLE = 'Mass above floatation'
    data_scalar.DATE = DATE_STR
    data_scalar.close()

    # -------------- iareagr ------------------
    data_scalar = Dataset(f'{output_path}/iareagr_AIS_DOE_MALI_{exp}.nc', 'w', format='NETCDF4_CLASSIC')
    data_scalar.createDimension('time', timeSteps_out)
    iareagrValues = data_scalar.createVariable('iareagr', 'f4', ('time'))
    timeValues = data_scalar.createVariable('time', 'f4', ('time'))
    for i in range(timeSteps_out):
        iareagrValues[i] = gia_snapshot[i]
        timeValues[i] = i * 365.0
    timeValues.units = 'days since 2015-01-01'
    timeValues.calendar = '365_day'
    timeValues.standard_name = 'time'
    timeValues.long_name = 'time'
    iareagrValues.standard_name = 'grounded_ice_sheet_area'
    iareagrValues.units = 'm^2'
    data_scalar.AUTHORS= AUTHOR_STR
    data_scalar.MODEL= 'MALI (MPAS-Albany Land Ice model)'
    data_scalar.GROUP = 'Los Alamos National Laboratory'
    data_scalar.VARIABLE = 'Grounded ice area'
    data_scalar.DATE = DATE_STR
    data_scalar.close()

    # -------------- iareafl ------------------
    data_scalar = Dataset(f'{output_path}/iareafl_AIS_DOE_MALI_{exp}.nc', 'w', format='NETCDF4_CLASSIC')
    data_scalar.createDimension('time', timeSteps_out)
    iareaflValues = data_scalar.createVariable('iareafl', 'f4', ('time'))
    timeValues = data_scalar.createVariable('time', 'f4', ('time'))
    for i in range(timeSteps_out):
        iareaflValues[i] = fia_snapshot[i]
        timeValues[i] = i * 365.0
    timeValues.units = 'days since 2015-01-01'
    timeValues.calendar = '365_day'
    timeValues.standard_name = 'time'
    timeValues.long_name = 'time'
    iareaflValues.standard_name = 'floating_ice_shelf_area'
    iareaflValues.units = 'm^2'
    data_scalar.AUTHORS= AUTHOR_STR
    data_scalar.MODEL= 'MALI (MPAS-Albany Land Ice model)'
    data_scalar.GROUP = 'Los Alamos National Laboratory'
    data_scalar.VARIABLE = 'Floating ice area'
    data_scalar.DATE = DATE_STR
    data_scalar.close()

    # -------------- tendacabf: this is a flux var
    data_scalar = Dataset(f'{output_path}/tendacabf_AIS_DOE_MALI_{exp}.nc', 'w', format='NETCDF4_CLASSIC')
    data_scalar.createDimension('time', timeSteps_out)
    tendacabfValues = data_scalar.createVariable('tendacabf', 'f4', ('time'))
    timeValues = data_scalar.createVariable('time', 'f4', ('time'))
    data_scalar.createDimension('bnds', 2)
    timebndsValues = data_scalar.createVariable('time_bnds', 'f4', ('time', 'bnds'))
    for i in range(timeSteps_out):
        tendacabfValues[i] = smb_avg[i] / 31536000
        timeValues[i] = (i + 0.5) * 365.0
        timebndsValues[i, 0] = i * 365.0
        timebndsValues[i, 1] = (i + 1) * 365.0
    timeValues.units = 'days since 2015-01-01'
    timeValues.calendar = '365_day'
    timeValues.standard_name = 'time'
    timeValues.long_name = 'time'
    tendacabfValues.standard_name = 'tendency_of_land_ice_mass_due_to_surface_mass_balance'
    tendacabfValues.units = 'kg s^-1'
    data_scalar.AUTHORS= AUTHOR_STR
    data_scalar.MODEL= 'MALI (MPAS-Albany Land Ice model)'
    data_scalar.GROUP = 'Los Alamos National Laboratory'
    data_scalar.VARIABLE = 'Total SMB flux'
    data_scalar.DATE = DATE_STR
    data_scalar.close()

    # -------------- tendlibmassbf: this is a flux var
    data_scalar = Dataset(f'{output_path}/tendlibmassbf_AIS_DOE_MALI_{exp}.nc', 'w', format='NETCDF4_CLASSIC')
    data_scalar.createDimension('time', timeSteps_out)
    tendlibmassbfValues = data_scalar.createVariable('tendlibmassbf', 'f4', ('time'))
    timeValues = data_scalar.createVariable('time', 'f4', ('time'))
    data_scalar.createDimension('bnds', 2)
    timebndsValues = data_scalar.createVariable('time_bnds', 'f4', ('time', 'bnds'))
    for i in range(timeSteps_out):
        tendlibmassbfValues[i] = bmb_avg[i] / 31536000
        timeValues[i] = (i + 0.5) * 365.0
        timebndsValues[i, 0] = i * 365.0
        timebndsValues[i, 1] = (i + 1) * 365.0
    timeValues.units = 'days since 2015-01-01'
    timeValues.calendar = '365_day'
    timeValues.standard_name = 'time'
    timeValues.long_name = 'time'
    tendlibmassbfValues.standard_name = 'tendency_of_land_ice_mass_due_to_basal_mass_balance '
    tendlibmassbfValues.units = 'kg s^-1'
    data_scalar.AUTHORS= AUTHOR_STR
    data_scalar.MODEL= 'MALI (MPAS-Albany Land Ice model)'
    data_scalar.GROUP = 'Los Alamos National Laboratory'
    data_scalar.VARIABLE = 'Total BMB flux'
    data_scalar.DATE = DATE_STR
    data_scalar.close()

    # -------------- tendlibmassbffl: this is a flux var
    data_scalar = Dataset(f'{output_path}/tendlibmassbffl_AIS_DOE_MALI_{exp}.nc', 'w',
                          format='NETCDF4_CLASSIC')
    data_scalar.createDimension('time', timeSteps_out)
    tendlibmassbfflValues = data_scalar.createVariable('tendlibmassbffl', 'f4', ('time'))
    timeValues = data_scalar.createVariable('time', 'f4', ('time'))
    data_scalar.createDimension('bnds', 2)
    timebndsValues = data_scalar.createVariable('time_bnds', 'f4', ('time', 'bnds'))
    for i in range(timeSteps_out):
        tendlibmassbfflValues[i] = bmb_avg[i] / 31536000
        timeValues[i] = (i + 0.5) * 365.0
        timebndsValues[i, 0] = i * 365.0
        timebndsValues[i, 1] = (i + 1) * 365.0
    timeValues.units = 'days since 2015-01-01'
    timeValues.calendar = '365_day'
    timeValues.standard_name = 'time'
    timeValues.long_name = 'time'
    tendlibmassbfflValues.standard_name = 'tendency_of_land_ice_mass_due_to_basal_mass_balance'
    tendlibmassbfflValues.units = 'kg s^-1'
    data_scalar.AUTHORS= AUTHOR_STR
    data_scalar.MODEL= 'MALI (MPAS-Albany Land Ice model)'
    data_scalar.GROUP = 'Los Alamos National Laboratory'
    data_scalar.VARIABLE = 'Total BMB flux beneath floating ice'
    data_scalar.DATE = DATE_STR
    data_scalar.close()

    # -------------- tendlicalvf: this is a flux var
    data_scalar = Dataset(f'{output_path}/tendlicalvf_AIS_DOE_MALI_{exp}.nc', 'w', format='NETCDF4_CLASSIC')
    data_scalar.createDimension('time', timeSteps_out)
    tendlicalvfValues = data_scalar.createVariable('tendlicalvf', 'f4', ('time'))
    timeValues = data_scalar.createVariable('time', 'f4', ('time'))
    data_scalar.createDimension('bnds', 2)
    timebndsValues = data_scalar.createVariable('time_bnds', 'f4', ('time', 'bnds'))
    for i in range(timeSteps_out):
        tendlicalvfValues[i] = -cfx_avg[i] / 31536000
        timeValues[i] = (i + 0.5) * 365.0
        timebndsValues[i, 0] = i * 365.0
        timebndsValues[i, 1] = (i + 1) * 365.0
    timeValues.units = 'days since 2015-01-01'
    timeValues.calendar = '365_day'
    timeValues.standard_name = 'time'
    timeValues.long_name = 'time'
    tendlicalvfValues.standard_name = 'tendency_of_land_ice_mass_due_to_calving'
    tendlicalvfValues.units = 'kg s^-1'
    data_scalar.AUTHORS= AUTHOR_STR
    data_scalar.MODEL= 'MALI (MPAS-Albany Land Ice model)'
    data_scalar.GROUP = 'Los Alamos National Laboratory'
    data_scalar.VARIABLE = 'Total calving flux'
    data_scalar.DATE = DATE_STR
    data_scalar.close()

    # -------------- tendlifmassbf: this is a flux var
        # this variable (Total calving and ice front melting front is not calculated in MALI) ??

    # -------------- tendligroundf: this is a flux var
    data_scalar = Dataset(f'{output_path}/tendligroundf_AIS_DOE_MALI_{exp}.nc', 'w', format='NETCDF4_CLASSIC')
    data_scalar.createDimension('time', timeSteps_out)
    tendligroundfValues = data_scalar.createVariable('tendligroundf', 'f4', ('time'))
    timeValues = data_scalar.createVariable('time', 'f4', ('time'))
    data_scalar.createDimension('bnds', 2)
    timebndsValues = data_scalar.createVariable('time_bnds', 'f4', ('time', 'bnds'))
    for i in range(timeSteps_out):
        tendligroundfValues[i] = gfx_avg[i] / 31536000
        timeValues[i] = (i + 0.5) * 365.0
        timebndsValues[i, 0] = i * 365.0
        timebndsValues[i, 1] = (i + 1) * 365.0
    timeValues.units = 'days since 2015-01-01'
    timeValues.calendar = '365_day'
    timeValues.standard_name = 'time'
    timeValues.long_name = 'time'
    tendligroundfValues.standard_name = 'tendency_of_grounded_ice_mass'
    tendligroundfValues.units = 'kg s^-1'
    data_scalar.AUTHORS= AUTHOR_STR
    data_scalar.MODEL= 'MALI (MPAS-Albany Land Ice model)'
    data_scalar.GROUP = 'Los Alamos National Laboratory'
    data_scalar.VARIABLE = 'Total grounding line flux'
    data_scalar.DATE = DATE_STR
    data_scalar.close()

    data.close()
