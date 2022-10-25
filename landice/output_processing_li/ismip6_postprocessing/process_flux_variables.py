"""
This script has functions that are needed to post-process and write flux
output variables from ISMIP6 simulations.
"""

from netCDF4 import Dataset
import xarray as xr
import numpy as np
from datetime import date
from subprocess import check_call
import os


def do_time_avg_flux_vars(input_file, output_file):
    """
    input_file: MALI simulation flux file that has the all time levels
    output_file: file with time-averaged fluxes
    """
    dataIn = xr.open_dataset(input_file)
    time = dataIn.dims['Time']
    xtime = dataIn['xtime'][:].values
    deltat = dataIn['deltat'][:]
    cellMask = dataIn['cellMask'][:,:]
    sfcMassBal = dataIn['sfcMassBalApplied'][:, :]
    basalMassBal = dataIn['basalMassBalApplied'][:, :]
    # basalHeatFlux = dataIn['basalHeatFlux'][:,:] HH: this needs to be outputted in the flux stream. Uncomment once the variable is in the output file. But there's also chance that this will be outputted in the state file
    # dHdt = dataIn['dHdt'][:,:] HH: this needs to be outputted in the flux stream. Uncomment once the variable is in the output file
    glFlux = dataIn['fluxAcrossGroundingLine'][:, :]

    # calculate averaged BMB flux beneath grounded and floating ice
    libmassbffl = basalMassBal * (cellMask[:, :] & 4) / 4
    libmassbfgr = basalMassBal * (1 - (cellMask[:, :] & 4) / 4)
    iceMask = (cellMask[:, :] & 2) / 2  # grounded: dynamic ice
    dataIn['libmassbffl'] = libmassbffl
    dataIn['libmassbfgr'] = libmassbfgr
    dataIn['iceMask'] = iceMask
    dataIn.to_netcdf(input_file, mode='a')

    # get an array of years that are not duplicative
    years = np.zeros(len(xtime))
    years[0] = int(xtime[0].decode("utf-8")[0:4])
    m = 1
    for n in np.arange(1, len(xtime)):
        yr_current = int(xtime[n].decode("utf-8")[0:4])
        if yr_current > years[m - 1]:
            years[m] = int(xtime[n].decode("utf-8")[0:4])
            m = m + 1
    years = np.trim_zeros(years)

    # create a temporary file to which processed data will be saved with the
    # right dimension
    ind = len(years) - 1
    command = ["ncks", "-O",
               "-d", f"Time,0,{str(ind)}",
               input_file,
               "tmp.nc"]

    check_call(command)

    # reduce the dim to yearTotal/yearInterval, and create the ouput file
    dataOut = xr.open_dataset('tmp.nc')
    for j in range(len(years)):

        sumYearSmb = 0
        sumYearBmb = 0
        # sumYearDHdt = 0
        # sumYearCF = 0
        sumYearGF = 0
        sumYearBHF = 0
        sumYearTime = 0
        sumYearBmbfl = 0
        sumYearBmbgr = 0
        sumIceMask = 0

        for i in range(time):

            yr_current = int(xtime[i].decode("utf-8")[0:4])
            if yr_current == years[j]:
                sumYearSmb = sumYearSmb + sfcMassBal[i, :] * deltat[i]
                sumYearBmbfl = sumYearBmbfl + libmassbffl[i, :] * deltat[i]
                sumYearBmbgr = sumYearBmbgr + libmassbfgr[i, :] * deltat[i]
                # sumYearDHdt = sumYearDHdt + dHdt[i, :] * deltat[i]
                # sumYearCF = sumYearCF + calvingFlux[i,:]*deltat[i]
                sumYearGF = sumYearGF + glFlux[i, :] * deltat[i]
                # sumYearBHF = sumYearBHF + basalHeatFlux[i,:]*deltat[i] # this might not be used because our BHF does not evolve with time
                sumYearTime = sumYearTime + deltat[i]

                sumIceMask = sumIceMask + iceMask[i,;]

        avgYearSmb = sumYearSmb / sumYearTime
        avgYearBmbfl = sumYearBmbfl / sumYearTime
        avgYearBmbgr = sumYearBmbgr / sumYearTime
        # avgYearDHdt = sumYearDHdt / sumYearTime
        # avgYearCF = sumYearCF/sumYearTime
        avgYearGF = sumYearGF / sumYearTime
        # avgYearBHF= sumYearBHF / sumYearTime
        maxIceMask = (sumIceMask>0) # Get mask for anywhere that had ice during this year

        dataOut['sfcMassBalApplied'][j, :] = avgYearSmb
        dataOut['libmassbffl'][j, :] = avgYearBmbfl
        dataOut['libmassbfgr'][j, :] = avgYearBmbgr
        # dataOut['dHdt'][j, :] = avgYearDHdt
        dataOut['fluxAcrossGroundingLine'][j, :] = avgYearGF
        # dataOut['basalHeatFlux'][j, :] = avgYearBHF

        # Generate a mask for ice extent
        dataOut['iceMask'][j,:] = maxIceMask

    dataOut.to_netcdf(output_file, mode='w')
    dataIn.close()
    dataOut.close()
    os.remove('tmp.nc')


def translate_GL_and_calving_flux_edge2cell(file_flux_time_avged,
                                            file_flux_on_cell,
                                            file_state):
    """
    file_flux_time_avged: time-averaged flux variables in MALI mesh
    (i.e., output file of the function do_time_avg_flux_vars)
    file_flux_on_cell: files with flux variables on the cell centers
    file_state: MALI state output file
    """

    data = xr.open_dataset(file_flux_time_avged, engine="netcdf4")
    data_state = xr.open_dataset(file_state, engine="netcdf4")
    nCells = data.dims['nCells']
    time = data.dims['Time']
    nEdgesOnCell = data['nEdgesOnCell'][:].values
    edgesOnCell = data['edgesOnCell'][:]
    cellsOnEdge = data['cellsOnEdge'][:].values
    dvEdge = data['dvEdge'][:].values
    areaCell = data['areaCell'][:].values
    deltat = data['deltat'][:].values
    # edgeMask = data['edgeMask'][:, :].values # this needs to be outputted as well. Once commented out, comment out L195 as well, and indent L197
    # dHdt = data['dHdt'][:, :].values # Uncomment this line once the dHdt variable is outputted in the stream and delete the line below
    dHdt = data['calvingThickness'][:, :].values  # delete this line once the dHdt variable is outputted in the stream. This line is used just for now to check the code
    fluxGLEdge = data['fluxAcrossGroundingLine'][:, :]
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
        data.to_netcdf(
            file_flux_on_cell)  # writing out to a netCDF file seems to be needed to save the newly added variable `fluxAcrossGroundingLineCell`.
        data.close()

        print("===done calving flux processing!===")


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

    spy = 365.0
    timeSpan = 1  # HH:check - is something that's needed?
    data_ismip6 = Dataset(ismip6_grid_file, 'r')
    var_x = data_ismip6.variables['x'][:]
    var_y = data_ismip6.variables['y'][:]

    data = Dataset(remapped_mali_flux_file, 'r')
    data.set_auto_mask(False)
    ice_mask = data.variables['sftgif'][:, :, :]
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
    DATE_STR = daTE.TOday().strftime("%d-%b-%Y")

    for i in range(timeSteps):
        mask = ice_mask[i, :, :]
        tmp = var_mali[i, :, :]
        tmp[mask == 0] = np.NAN
        dataValues[i, :, :] = tmp
        timeValues[i] = (i + 0.5) * timeSpan * spy
        timebndsValues[i, 0] = i * timeSpan * spy
        timebndsValues[i, 1] = (i + 1) * timeSpan * spy

    for i in range(latN):
        xValues[i] = var_x[i]

    for i in range(lonN):
        yValues[i] = var_y[i]

    dataValues.standard_name = var_std_name
    dataValues.units = var_units
    timeValues.bounds = 'time_bnds'
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
    dataOut.GROUP = 'Los Alamos National Laboratory, Department of Energy'
    dataOut.VARIABLE = var_varname
    dataOut.DATE = DATE_STR
    dataOut.close()


def generate_output_2d_flux_vars(file_remapped_mali_flux,
                                 file_remapped_mali_state,
                                 ismip6_grid_file,
                                 exp, output_path):
    """
    file_remapped_mali_flux: flux output file on mali mesh remapped
    onto the ismip6 grid
    file_remapped_mali_state: flux output file on mali mesh remapped onto the
    ismip6 grid
    ismip6_grid_file: ismip6 original file
    exp: ISMIP6 experiment name
    output_path: path to which the final output files are saved
    """

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

    # ----------- dlithkdt ------------------ # Uncomment the section once dHdt is outputted in the stream
    # write_netcdf_2d_flux_vars('dHdt', 'dlithkdt',
    #                           'tendency_of_land_ice_thickness',
    #                           'm s-1',
    #                           'Ice thickness imbalance',
    #                           file_remapped_mali_flux,
    #                           ismip6_grid_file,
    #                           exp, output_path)

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
    write_netcdf_2d_flux_vars('fluxAcrossGroundingLineCell', 'ligroundf',
                              'land_ice_specific_mass_flux_at_grounding_line',
                              'kg m-2 s-1',
                              'Grounding line flux',
                              file_remapped_mali_flux,
                              ismip6_grid_file,
                              exp, output_path)
