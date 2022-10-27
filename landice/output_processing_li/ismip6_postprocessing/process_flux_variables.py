"""
This script has functions that are needed to post-process and write flux
output variables from ISMIP6 simulations.
"""

from netCDF4 import Dataset
import xarray as xr
import numpy as np
from datetime import date
from subprocess import check_call
import os, sys


def do_time_avg_flux_vars(input_file, output_file):
    """
    input_file: MALI simulation flux file that has the all time levels
    output_file: file with time-averaged fluxes
    """
    print("Starting time averaging of flux variables")
    input_file_tmp = 'flux_input_tmp.nc'
    dataIn = xr.open_dataset(input_file, chunks={'Time': 5}, decode_cf=False) # need decode_cf=False to prevent xarray from reading daysSinceStart as a timedelta type.
    #print(dataIn.info)
    #print(dataIn.data_vars)
    time = dataIn.dims['Time']
    nCells = dataIn.dims['nCells']
    xtimeIn = dataIn['xtime'][:].values
    #print(xtimeIn)
    xtime = []
    for i in range(time):
        xtime.append(xtimeIn[i].tostring().decode('utf-8').strip().strip('\x00'))
    #print(xtime)
    deltat = dataIn['deltat'][:]
    del dataIn.daysSinceStart.attrs['units'] # need this line to prevent xarray from reading daysSinceStart as a timedelta type.
    daysSinceStart = dataIn['daysSinceStart'][:]
    cellMask = dataIn['cellMask'][:,:]
    sfcMassBal = dataIn['sfcMassBalApplied'][:, :]
    basalMassBal = dataIn['basalMassBalApplied'][:, :]
    dHdt = dataIn['dHdt'][:,:] / (3600.0 * 24.0 * 365.0) # convert units to m/s
    glFlux = dataIn['fluxAcrossGroundingLineOnCells'][:, :]

    # calculate averaged BMB flux beneath grounded and floating ice
    libmassbffl = basalMassBal * (cellMask[:, :] & 4) / 4
    libmassbfgr = basalMassBal * (1 - (cellMask[:, :] & 4) / 4)
    iceMask = (cellMask[:, :] & 2) / 2  # grounded: dynamic ice
    dataIn['libmassbffl'] = libmassbffl
    dataIn['libmassbfgr'] = libmassbfgr
    dataIn['iceMask'] = iceMask
    #print("    saving modified input file")
    #dataIn.to_netcdf(input_file_tmp, mode='w')
    #print("    finished saving modified input file")


    # Figure out some timekeeping stuff - using netCDF4 b/c xarray is a nightmare
    fin = Dataset(input_file, 'r')
    simulationStartTime = fin.variables['simulationStartTime'][:].tostring().decode('utf-8').strip().strip('\x00')
    fin.close()
    simulationStartDate = simulationStartTime.split("_")[0]
    if simulationStartDate[5:10] != '01-01':
        sys.exit("Error: simulationStartTime for flux file is not on Jan. 1.")
    refYear = int(simulationStartDate[0:4])
    startYr = refYear + np.floor(daysSinceStart[0] / 365.0) # using floor here because we might not have output at jan 1, but we'll definitely have at least one time level per year
    finalYr = refYear + daysSinceStart[-1] / 365.0
    if (daysSinceStart[-1] / 365.0 != daysSinceStart[-1] // 365):
        sys.exit(f"Error: final time of flux output file is not on Jan. 1.: daysSinceStart={daysSinceStart[-1]}, xtime={xtime[-1]}" )
    print(f"simulationStartTime={simulationStartTime}; simulationStartDate={simulationStartDate}; refYear={refYear}")
    print(f"start year={startYr}; final year={finalYr}")

    # get an array of years that are not duplicative
    decYears = refYear + daysSinceStart/365.0
    #years = np.floor(decYears - 1.0e-10)  # this is the "owning" year; Jan 1 belongs to the previous year, so offset decYears by small amount
    #years[0] = int(xtime[0].decode("utf-8")[0:4])
    ##years = np.trim_zeros(years)
    years = np.arange(startYr, finalYr) # we don't want the final year in the time array as a year to process - it's actually the end point of the previous year

    timeBndsMin = np.ones((len(years),)) * 1.0e36
    timeBndsMax = np.ones((len(years),)) * -1.0e36

    avgSmb = np.zeros((len(years), nCells)) * np.nan
    avgBmbfl = np.zeros((len(years), nCells)) * np.nan
    avgBmbgr = np.zeros((len(years), nCells)) * np.nan
    avgDHdt = np.zeros((len(years), nCells)) * np.nan
    avgGF = np.zeros((len(years), nCells)) * np.nan
    maxIceMask = np.zeros((len(years), nCells), dtype=np.int) * np.nan

    print("    begin looping over years")
    for j in range(len(years)):
        # we want time bounds to span the full year
        timeBndsMin[j] = (years[j] - refYear) * 365.0
        timeBndsMax[j] = (years[j]+1.0 - refYear) * 365.0
        print(f"     year index: {j}, year={years[j]}; timeBindsMin={timeBndsMin[j]}, timeBndsMax={timeBndsMax[j]}")
        sumYearSmb = 0
        sumYearBmb = 0
        sumYearDHdt = 0
        # sumYearCF = 0
        sumYearGF = 0
        sumYearBHF = 0
        sumYearTime = 0
        sumYearBmbfl = 0
        sumYearBmbgr = 0
        sumIceMask = 0

        timeBndMin = 1.0e36
        timeBndMax = -1.0e36
        for i in range(time):

            if decYears[i] > years[j] and decYears[i] <= years[j]+1.0:
                sumYearSmb = sumYearSmb + sfcMassBal[i, :] * deltat[i]
                sumYearBmbfl = sumYearBmbfl + libmassbffl[i, :] * deltat[i]
                sumYearBmbgr = sumYearBmbgr + libmassbfgr[i, :] * deltat[i]
                sumYearDHdt = sumYearDHdt + dHdt[i, :] * deltat[i]
                # sumYearCF = sumYearCF + calvingFlux[i,:]*deltat[i]
                sumYearGF = sumYearGF + glFlux[i, :] * deltat[i]
                sumYearTime = sumYearTime + deltat[i]

                sumIceMask = sumIceMask + iceMask[i,:]

                print(f"         year={years[j]}, decYears={decYears[i]}, daysSinceStart={daysSinceStart[j]}, xtime={xtime[i]}")


        avgSmb[j,:] = sumYearSmb / sumYearTime
        avgBmbfl[j,:] = sumYearBmbfl / sumYearTime
        avgBmbgr[j,:] = sumYearBmbgr / sumYearTime
        avgDHdt[j,:] = sumYearDHdt / sumYearTime
        # avgYearCF = sumYearCF/sumYearTime
        avgGF[j,:] = sumYearGF / sumYearTime
        maxIceMask[j,:] = (sumIceMask>0) # Get mask for anywhere that had ice during this year


    print("    write time averaged values")

    print(f"avg shape={avgSmb.shape}, time shape={timeBndsMin.shape}")
    out_data_vars = {
                     'sfcMassBalApplied': (['Time', 'nCells'], avgSmb),
                     'libmassbffl':       (['Time', 'nCells'], avgBmbfl),
                     'libmassbfgr':       (['Time', 'nCells'], avgBmbgr),
                     'dHdt':              (['Time', 'nCells'], avgDHdt),
                     'fluxAcrossGroundingLineOnCells': (['Time', 'nCells'], avgGF),
                     'iceMask':        (['Time', 'nCells'], maxIceMask),
                     'timeBndsMin': (['Time'], timeBndsMin),
                     'timeBndsMax': (['Time'], timeBndsMax),
                     'simulationStartTime': dataIn['simulationStartTime']
                     }
    out_coords = {
                  'Time':   (['Time'], (timeBndsMin+timeBndsMax)/2.0)
                 }


    dataOut = xr.Dataset(data_vars=out_data_vars, coords=out_coords)
    dataOut.to_netcdf(output_file, mode='w')
    dataIn.close()

    os.remove(input_file_tmp)


def translate_GL_and_calving_flux_edge2cell(file_flux_time_avged,
                                            file_flux_on_cell):
    """
    file_flux_time_avged: time-averaged flux variables in MALI mesh
    (i.e., output file of the function do_time_avg_flux_vars)
    file_flux_on_cell: files with flux variables on the cell centers
    file_state: MALI state output file
    """
    print("Starting translation of GL and calving fluxes")

    data = xr.open_dataset(file_flux_time_avged, engine="netcdf4")
    nCells = data.dims['nCells']
    time = data.dims['Time']
    nEdgesOnCell = data['nEdgesOnCell'][:].values
    edgesOnCell = data['edgesOnCell'][:].values
    cellsOnEdge = data['cellsOnEdge'][:].values
    dvEdge = data['dvEdge'][:].values
    areaCell = data['areaCell'][:].values
    deltat = data['deltat'][:].values
    # edgeMask = data['edgeMask'][:, :].values # this needs to be outputted as well. Once commented out, comment out L195 as well, and indent L197
    # dHdt = data['dHdt'][:, :].values # Uncomment this line once the dHdt variable is outputted in the stream and delete the line below
    dHdt = data['calvingThickness'][:, :].values  # delete this line once the dHdt variable is outputted in the stream. This line is used just for now to check the code
    fluxGLEdge = data['fluxAcrossGroundingLine'].load()
    fluxGLCell_array = np.zeros((time, nCells))

    print("=== starting the grounding line flux processing ===")
    for i in range(nCells):
        edgeId = edgesOnCell[i, :nEdgesOnCell[i]] - 1
        dvCell = np.sum(fluxGLEdge[:, edgeId], axis=1) * 910.0
        fluxGLCell_array[:, i] = dvCell

    data['fluxAcrossGroundingLineCell'] = (('Time', 'nCells'), fluxGLCell_array)
    print("=== done grounding line processing ===")

    doCalving = False
    # Calving processing will eventually be moved online to MALI
    # Disable for now because we don't want to waste time on this for
    # fixed calving front runs anyway.
    if doCalving:
        print("===starting the calving flux processing===")
        rho_i = 910.0
        h0 = 11.1

        # get calving thickness values
        calvingThickness = data['calvingThickness'][:, :].values
        assert time == len(deltat)

        # create and initialize a new data array for calvingFluxArray
        calvingFluxArray = data['calvingThickness'].copy()
        calvingFluxArray = calvingFluxArray * 0.0

        # This part seems to have been used to get `edgeMask`
        # thickness = data_state['thickness'][:, :].values
        # thicknessRecovered = np.copy(dHdt)
        # for i in range(time):
        #     thickness[:] = thickness[:] + dHdt[i, :] * deltat[i]
        #     thickness[np.where(thickness < 0)] = 0.0
        #     thicknessRecovered[i, :] = thickness[:]

        for t in range(time):

            index_cf = np.where(calvingThickness[t, :] > 0)[0]
            for i in index_cf:

                ne = nEdgesOnCell[i]
                edgeId = edgesOnCell[i, :nEdgesOnCell[i]] - 1

                dvEdgeSum = 0.0
                for j in range(ne):
                    neighborCellId = cellsOnEdge[edgeId[j], :] - 1

                    # if ((edgeMask[0, edgeId[j]] & 16) / 16 and (edgeMask[0, edgeId[j]] & 4) / 4):  # uncomment this line and indent L197 once edgeMask is acquired.
                    # 16: dynamic margin & 4: floating
                    dvEdgeSum = dvEdgeSum + dvEdge[edgeId[j]]

                    # below is Tong's script that uses thickness recovered to calculate dvEdgeSum. Tong writes "I haven't tested the below part yet"
                    # else:
                    #    boolArray1 =  ((thicknessRecovered[:,neighborCellId[0]] > h0) | (thicknessRecovered[:,neighborCellId[1]] > h0)) \
                    #            & ~((thicknessRecovered[:,neighborCellId[0]] > h0) & (thicknessRecovered[:,neighborCellId[1]] > h0))
                    #    boolArray2 = (smbApplied[:,neighborCellId[0]] > 0) & (smbApplied[:,neighborCellId[1]] > 0)
                    #    boolArray = boolArray1 & boolArray2
                    #    dvEdgeSum = dvEdgeSum + dvEdge[edgeId[j]]*boolArray

                if dvEdgeSum > 0:
                    calvingVelDivideThickness = calvingThickness[t, i] * areaCell[i] / (deltat[t] * (dvEdgeSum + 1e-10))
                    # put 1e-10 here in case of dvEdgeSum = 0, and change the value back to zero in below
                else:
                    calvingVelDivideThickness = 0.0
                    # unit: m^2/s

                calvingFlux = rho_i * calvingVelDivideThickness
                # unit: kg/m/s
                calvingFluxArray[t, i] = calvingFlux

        data['calvingFlux'] = calvingFluxArray

        print("===done calving flux processing!===")

    data.to_netcdf(
    file_flux_on_cell)  # writing out to a netCDF file seems to be needed to save the newly added variable `fluxAcrossGroundingLineCell`.
    data.close()

def write_netcdf_2d_flux_vars(mali_var_name, ismip6_var_name, var_std_name,
                              var_units, var_varname,
                              remapped_mali_flux_file,
                              ismip6_grid_file,
                              exp, output_path):

    """
    mali_var_name: variable name on MALI side
    ismip6_var_name: variable name required by ISMIP6
    var_std_name: standard variable name
    var_units: variable units
    var_varname: variable variable name
    remapped_mali_flux_file: mali flux file remapped on the ISMIP6 grid
    ismip6_grid_file: original ISMIP6 file
    exp: experiment name
    output_path: output path to which the output files will be saved
    """


    data_ismip6 = Dataset(ismip6_grid_file, 'r')
    var_x = data_ismip6.variables['x'][:]
    var_y = data_ismip6.variables['y'][:]

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

    dataOut = Dataset(f'{output_path}/{ismip6_var_name}_AIS_DOE_MALI_{exp}.nc',
                      'w', format='NETCDF4_CLASSIC')
    dataOut.createDimension('time', timeSteps)
    dataOut.createDimension('bnds', 2)
    timebndsValues = dataOut.createVariable('time_bnds', 'f4', ('time', 'bnds'))
    dataOut.createDimension('x', lonN)
    dataOut.createDimension('y', latN)
    dataValues = dataOut.createVariable(ismip6_var_name, 'f4',
                                       ('time', 'y', 'x'), fill_value=np.NAN)
    xValues = dataOut.createVariable('x', 'f4', ('x'))
    yValues = dataOut.createVariable('y', 'f4', ('y'))
    timeValues = dataOut.createVariable('time', 'f4', ('time'))

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
                                 ismip6_grid_file,
                                 exp, output_path):
    """
    file_remapped_mali_flux: flux output file on mali mesh remapped
    onto the ismip6 grid
    ismip6 grid
    ismip6_grid_file: ismip6 original file
    exp: ISMIP6 experiment name
    output_path: path to which the final output files are saved
    """

    print("Writing 2d flux variables")
    # ----------- acabf ------------------
    write_netcdf_2d_flux_vars('sfcMassBalApplied', 'acabf',
                              'land_ice_surface_specific_mass_balance_flux',
                              'kg m-2 s-1', 'Surface mass balance flux',
                              file_remapped_mali_flux,
                              ismip6_grid_file,
                              exp, output_path)

    # ----------- libmassbffl ------------------
    write_netcdf_2d_flux_vars('libmassbffl', 'libmassbffl',
                              'land_ice_basal_specific_mass_balance_flux',
                              'kg m-2 s-1',
                              'Basal mass balance flux beneath floating ice',
                              file_remapped_mali_flux,
                              ismip6_grid_file,
                              exp, output_path)

    # ----------- libmassbfgr ------------------
    write_netcdf_2d_flux_vars('libmassbfgr', 'libmassbfgr',
                              'land_ice_basal_specific_mass_balance_flux',
                              'kg m-2 s-1',
                              'Basal mass balance flux beneath grounded ice',
                              file_remapped_mali_flux,
                              ismip6_grid_file,
                              exp, output_path)

    # ----------- dlithkdt ------------------
    write_netcdf_2d_flux_vars('dHdt', 'dlithkdt',
                              'tendency_of_land_ice_thickness',
                              'm s-1',
                              'Ice thickness imbalance',
                              file_remapped_mali_flux,
                              ismip6_grid_file,
                              exp, output_path)

    # ----------- licalvf ------------------
    write_netcdf_2d_flux_vars('calvingFlux', 'licalvf',
                              'land_ice_specific_mass_flux_due_to_calving',
                              'kg m-2 s-1',
                              'Calving flux',
                              file_remapped_mali_flux,
                              ismip6_grid_file,
                              exp, output_path)

    # ----------- lifmassbf ------------------
        # Need to calculated it: ice front melt and calving flux


    # ----------- ligroundf ------------------
    write_netcdf_2d_flux_vars('fluxAcrossGroundingLineOnCells', 'ligroundf',
                              'land_ice_specific_mass_flux_at_grounding_line',
                              'kg m-2 s-1',
                              'Grounding line flux',
                              file_remapped_mali_flux,
                              ismip6_grid_file,
                              exp, output_path)
