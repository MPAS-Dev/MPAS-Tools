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
import warnings


def do_time_avg_flux_vars(input_file, output_file):
    """
    input_file: MALI simulation flux file that has the all time levels
    output_file: file with time-averaged fluxes
    """
    print("Starting time averaging of flux variables")
    dataIn = xr.open_dataset(input_file, decode_cf=False) # need decode_cf=False to prevent xarray from reading daysSinceStart as a timedelta type.
    if 'units' in dataIn.daysSinceStart.attrs:  # make have been removed in a previous step, so check if it exists
        del dataIn.daysSinceStart.attrs['units'] # need this line to prevent xarray from reading daysSinceStart as a timedelta type.

    time = dataIn.dims['Time']
    nCells = dataIn.dims['nCells']
    xtimeIn = dataIn['xtime'][:].values
    #print(xtimeIn)
    xtime = []
    for i in range(time):
        xtime.append(xtimeIn[i].tostring().decode('utf-8').strip().strip('\x00'))
    #print(xtime)
    deltat = dataIn['deltat'][:]
    daysSinceStart = dataIn['daysSinceStart'][:]
    cellMask = dataIn['cellMask'][:,:]
    sfcMassBal = dataIn['sfcMassBalApplied'][:, :]
    floatingBasalMassBalApplied = dataIn['floatingBasalMassBalApplied'][:, :]
    groundedBasalMassBalApplied = dataIn['groundedBasalMassBalApplied'][:, :]
    dHdt = dataIn['dHdt'][:,:] / (3600.0 * 24.0 * 365.0) # convert units to m/s
    glFlux = dataIn['fluxAcrossGroundingLineOnCells'][:, :]
    calvingFlux = dataIn['calvingFlux'][:, :]
    faceMeltAndCalvingFlux = dataIn['faceMeltAndCalvingFlux'][:, :]

    iceMask = (cellMask[:, :] & 2) / 2  # grounded: dynamic ice

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
    avgCF = np.zeros((len(years), nCells)) * np.nan
    avgCFandFM = np.zeros((len(years), nCells)) * np.nan
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
        sumYearCF = 0
        sumYearCFandFM = 0
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
                sumYearBmbfl = sumYearBmbfl + floatingBasalMassBalApplied[i, :] * deltat[i]
                sumYearBmbgr = sumYearBmbgr + groundedBasalMassBalApplied[i, :] * deltat[i]
                sumYearDHdt = sumYearDHdt + dHdt[i, :] * deltat[i]
                sumYearCF = sumYearCF + calvingFlux[i,:] * deltat[i]
                sumYearCFandFM = sumYearCFandFM + faceMeltAndCalvingFlux[i,:] * deltat[i]
                sumYearGF = sumYearGF + glFlux[i, :] * deltat[i]
                sumYearTime = sumYearTime + deltat[i]

                sumIceMask = sumIceMask + iceMask[i,:]

                print(f"         year={years[j]}, decYears={decYears[i]}, daysSinceStart={daysSinceStart[j]}, xtime={xtime[i]}")


        avgSmb[j,:] = sumYearSmb / sumYearTime
        avgBmbfl[j,:] = sumYearBmbfl / sumYearTime
        avgBmbgr[j,:] = sumYearBmbgr / sumYearTime
        avgDHdt[j,:] = sumYearDHdt / sumYearTime
        avgCF[j,:] = sumYearCF / sumYearTime
        avgCFandFM[j,:] = sumYearCFandFM / sumYearTime
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
                     'calvingFlux': (['Time', 'nCells'], avgCF),
                     'faceMeltAndCalvingFlux': (['Time', 'nCells'], avgCFandFM),
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


def clean_flux_fields_before_time_averaging(file_input, file_mesh,
                                            file_output):
    """
    Convert the MALI output field calvingThickness to the ISMIP6 variable
    licalvf and apply bounds checking on BMB, where some crazy values occasionally occur.
    """

    debug_face_melt_flux = False
    data = xr.open_dataset(file_input, decode_cf=False) # need decode_cf=False to prevent xarray from reading daysSinceStart as a timedelta type.
    if 'units' in data.daysSinceStart.attrs:
        del data.daysSinceStart.attrs['units'] # need this line to prevent xarray from reading daysSinceStart as a timedelta type.
    time = data.dims['Time']
    nCells = data.dims['nCells']
    nEdgesOnCell = data['nEdgesOnCell'][:].values
    edgesOnCell = data['edgesOnCell'][:].values
    cellsOnCell = data['cellsOnCell'][:].values
    dvEdge = data['dvEdge'][:].values
    areaCell = data['areaCell'][:].values
    xCell = data['xCell'][:].values
    yCell = data['yCell'][:].values
    deltat = data['deltat'][:].values
    thickness = data['thickness'][:].values
    surfaceSpeed = data['surfaceSpeed'][:].values
    if 'bedTopography' in data:
        bedTopography = data['bedTopography'][:].values
        print('bedTopography field found; using bedTopography at all time levels.')
    else:
        data_mesh = xr.open_dataset(file_mesh)
        bedTopography = data_mesh['bedTopography'][:].values
        print('No bedTopography field found; using bedTopography from mesh file.')
    if 'calvingThicknessFromThreshold' in data:
        calvingThicknessFromThreshold = data['calvingThicknessFromThreshold'][:, :].values
    else:
        print('WARNING: No calvingThicknessFromThreshold field found; creating a field populated with zeros.')
        calvingThicknessFromThreshold = thickness.copy() * 0.0

    rho_i = 910.0

    print("===starting cleaning floatingBasalMassBalApplied===")
    # We've encountered a few enormous BMB values.  Until we solve where that is coming from,
    # set them to something reasonable.
    floatingBasalMassBalApplied = data['floatingBasalMassBalApplied'][:, :].values
    for t in range(1, time):
        if t%20 == 0:
            print(f"    Time: {t+1} / {time}")
        # Set large negative BMB values (1 m/s) to the equivalent of the thickness from the previous time step
        # (Commented line here picked up too many places)
        #ind = np.nonzero((-floatingBasalMassBalApplied[t,:]/rho_i*deltat[t] > thickness[t-1,:]) * (thickness[t-1,:]>0.0))[0]
        ind = np.nonzero(floatingBasalMassBalApplied[t,:]/rho_i < -1.0)[0]
        if len(ind) > 0:
            print(f"Fixing {len(ind)} cells with large negative floating BMB values.", ind)
            floatingBasalMassBalApplied[t, ind] = -thickness[t-1, ind] / deltat[t] * rho_i

        ind = np.nonzero(floatingBasalMassBalApplied[t,:]/rho_i > 1.0)[0]
        if len(ind) > 0:
            print(f"Fixing {len(ind)} cells with large positive floating BMB values.", ind)
            ind2 = np.nonzero(floatingBasalMassBalApplied[t,:]/rho_i <= 1.0)[0]
            maxGoodBMB = floatingBasalMassBalApplied[t, ind2].max()
            floatingBasalMassBalApplied[t, ind] = maxGoodBMB

    print("===done cleaning floatingBasalMassBalApplied===")

    assert time == len(deltat)

    calvingVelocity = data['calvingVelocity'][:, :].values

    # create and initialize a new data array for calvingFluxArray
    calvingFluxArray = data['calvingVelocity'].copy() * 0.0
    thresholdFlux = data['calvingVelocity'].copy() * 0.0
    calvingThickness = data['calvingThickness'][:, :].values
    print("===starting facemelt flux processing===")

    # create and initialize a new data array for faceMeltFluxArray
    # (copied from calving code below)
    # Some runs won't have this output field, so assume if field is not present
    # that facemelting was not enabled
    faceMeltFluxArray = data['calvingVelocity'].copy() * 0.0
    if 'faceMeltSpeed' in data:
        faceMeltSpeed = data['faceMeltSpeed'][:, :].values
        # faceMeltSpeed is defined below the water line, but face-melting is
        # applied to the full ice thickness, so the effective speed is
        # averaged over the full thickness from the previous time step.
        # Note that this calculation assumes that bedTopography is constant in time,
        # that config_sea_level = 0, and that faceMeltSpeed is only valid for
        # grounded cells, i.e., that bedTopography and lowerSurface are equivalent
        # (which is currently the case).

        faceMeltingThickness = data['faceMeltingThickness'][:, :].values
        faceMeltSpeedVertAvg = faceMeltingThickness.copy() * 0.0
        # Fields for validation and debugging
        if debug_face_melt_flux:
            deltat_array = np.tile(deltat,  (np.shape(faceMeltSpeed)[1],1)).transpose()
            # Cleaned field for debugging and validation
            faceMeltingThicknessCleaned = faceMeltingThickness.copy()
        for t in range(time):
            if t%20 == 0:
                print(f"    Time: {t+1} / {time}")

            if 'bedTopography' in data:
                bed = bedTopography[t,:] # have value per time level
            else:
                bed = bedTopography[0,:] # just have a single value

            prev_t = max(t-1, 0)  # ensure that index_cf never uses thickness from last (-1) time step
            index_cf = np.where((faceMeltingThickness[t, :] > 0.0) * (bed[:] < 0.0) *
                                (faceMeltingThickness[t, :] != thickness[prev_t, :]) * 
                                (thickness[prev_t, :] > 0.))[0]
            for i in index_cf:
                # faceMeltSpeed is calculated for ice below water line, but needs to be aplied
                # to full ice thickness, so we need a vertically averaged speed. Also ensure that
                # the vertically averaged speed is never > faceMeltSpeed due to small ice thickness.
                # This may be slightly inaccurate on the very first time step.
                faceMeltSpeedVertAvg[t,i] = faceMeltSpeed[t, i] * np.abs(bed[i] / thickness[prev_t, i])
                faceMeltSpeedVertAvg[t,i] = min(faceMeltSpeedVertAvg[t,i], faceMeltSpeed[t, i])
                # Use this cell if it has nonzero faceMeltingThickness because faceMeltSpeed
                # is defined everywhere, but only applied on grounded ice
                if faceMeltingThickness[t,i] > 0.0:
                    faceMeltFluxArray[t,i] = faceMeltSpeedVertAvg[t,i] * rho_i # convert to proper units
            # Push mass removed from stranded non-dynamic cells into calving
            index_stranded_cell_cleanup = np.where(faceMeltingThickness[t, :] == thickness[prev_t, :])[0]
            calvingThicknessFromThreshold[t, index_stranded_cell_cleanup] += faceMeltingThickness[t, index_stranded_cell_cleanup]
            if debug_face_melt_flux:
                faceMeltingThicknessCleaned[t, index_stranded_cell_cleanup] -= faceMeltingThickness[t, index_stranded_cell_cleanup]
            # This is just for debugging and validation
    print("===done facemelt flux processing!===")

    print("===starting the calving flux processing===")

    for t in range(time):
        if t%20 == 0:
            print(f"    Time: {t+1} / {time}")

        if 'bedTopography' in data:
            bed = bedTopography[t,:] # have value per time level
        else:
            bed = bedTopography[0,:] # just have a single value

        index_cf = np.where((calvingVelocity[t, :] > 0.0) * (bed[:] < 0.0))[0]
        for i in index_cf:
            ne = nEdgesOnCell[i]
            for j in range(ne):
                neighborCellId = cellsOnCell[i, j] - 1
                # Use this cell if it has a neighbor with zero calvingVelocity that is below sea level
                if calvingVelocity[t,neighborCellId] == 0.0 and bed[neighborCellId] < 0.0:
                    calvingFluxArray[t,i] = calvingVelocity[t,i] * rho_i # convert to proper units
                    continue # no need to keep searching the neighbors of this cell


        # we may need to add on threshold calving too
        if 'calvingThicknessFromThreshold' in data:
            index_cf = np.where(calvingThicknessFromThreshold[t, :] > 0.0)[0]
        else:
            index_cf = []
        if len(index_cf) > 0:
            thresholdBoundary = np.zeros((nCells,), 'i')
            thresholdBoundaryAssignedVolume = np.zeros((nCells,))
            thresholdBoundarySummedThickness = np.zeros((nCells,))
            thresholdBoundaryContributors = np.zeros((nCells,))
            thresholdBoundaryLength = np.zeros((nCells,))
            thresholdSpeed = np.zeros((nCells,))
            # First make list of boundary cells calved
            for i in index_cf:
                ne = nEdgesOnCell[i]
                for j in range(ne):
                    neighborCellId = cellsOnCell[i, j] - 1
                    if thickness[t,neighborCellId] > 0.0 and bed[neighborCellId] < 0.0 and calvingThicknessFromThreshold[t,neighborCellId] == 0.0:
                        thresholdBoundary[i] = 1
                        thresholdBoundaryLength[i] += dvEdge[edgesOnCell[i,j]-1]
            bdyIndices = np.where(thresholdBoundary == 1)[0]
            print(f"Found {len(index_cf)} cells with threshold calving at time {t}; {len(bdyIndices)} are boundary cells.")
            if len(bdyIndices) == 0:
                print(f"0 boundary cells were found; skipping to next time step")
                continue
            # Now loop over all threshold cells and assign their volume to the nearest boundary cell
            for i in index_cf:
                if thresholdBoundary[i] == 1:
                    ownerIdx = i # often the cell is its own owner, so check before doing the more expensive search
                else:
                    ownerIdx = bdyIndices[np.argmin((xCell[i]-xCell[bdyIndices])**2 + (yCell[i]-yCell[bdyIndices])**2)]
                thresholdBoundaryAssignedVolume[ownerIdx] += calvingThicknessFromThreshold[t,i] * areaCell[i]
                thresholdBoundarySummedThickness[ownerIdx] += calvingThicknessFromThreshold[t,i]
                thresholdBoundaryContributors[ownerIdx] += 1
            #print(thresholdBoundaryAssignedVolume.sum(), (calvingThicknessFromThreshold[t,:]*areaCell[:]).sum())
            diff = np.absolute(thresholdBoundaryAssignedVolume.sum() - (calvingThicknessFromThreshold[t,:]*areaCell[:]).sum())
            if diff < 1.0:
                warnings.warn(f"Difference between assigned `thresholdBoundaryAssignedVolume` value and "
                              f"`calvingThicknessFromThreshold` threshold value is less than 1: {diff} < 1.0")
            #for i in bdyIndices:
                #print(f"length={thresholdBoundaryLength[i]}, vol={thresholdBoundaryAssignedVolume[i]}, sumthk={thresholdBoundarySummedThickness[i]}, num={thresholdBoundaryContributors[i]}, meanthk={thresholdBoundarySummedThickness[i]/thresholdBoundaryContributors[i]}")
            # Finally calculate licalvf for each boundary cell and add to whatever was already there
            thresholdSpeed[bdyIndices] = thresholdBoundaryAssignedVolume[bdyIndices] / \
                    (thresholdBoundarySummedThickness[bdyIndices] / thresholdBoundaryContributors[bdyIndices] * \
                     thresholdBoundaryLength[bdyIndices]) / \
                    deltat[t]    # units of m/s
            # Our estimated threshold speed is really a retreat speed. So to get calving speed, add on the advective speed
            thresholdFlux[t, bdyIndices] += (thresholdSpeed[bdyIndices] + surfaceSpeed[t,bdyIndices]) * rho_i
            calvingFluxArray[t,bdyIndices] += thresholdFlux[t,bdyIndices]

    data['calvingFlux'] = calvingFluxArray  # Note: thresholdFlux was already added in above
    data['thresholdFlux'] = thresholdFlux  # this is just written for diagnostic purposes.  It's not actually sent to ISMIP6.
    data['faceMeltAndCalvingFlux'] = faceMeltFluxArray + calvingFluxArray  # ismip6 only wants the combined fields for face-melt
    print("===done calving flux processing!===")
    if debug_face_melt_flux:
        print('debug_face_melt_flux is True, so I assume you want a breakpoint' +
              ' to check fluxes. Just type continue when you want to proceed.')
        breakpoint()
    data.to_netcdf(file_output) # copy of the input file with new vars added
    data.close()

def write_netcdf_2d_flux_vars(mali_var_name, ismip6_var_name, var_std_name,
                              var_units, var_varname, remapped_mali_flux_file,
                              res_ismip6_grid, ismip6_grid_file,
                              exp, output_path):

    """
    mali_var_name: variable name on MALI side
    ismip6_var_name: variable name required by ISMIP6
    var_std_name: standard variable name
    var_units: variable units
    var_varname: variable variable name
    remapped_mali_flux_file: mali flux file remapped on the ISMIP6 grid
    res_ismip6_grid: resolution of the ISMIP6 grid in kilometers
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
    timebndsValues = dataOut.createVariable('time_bnds', 'd', ('time', 'bnds'))
    dataOut.createDimension('x', lonN)
    dataOut.createDimension('y', latN)
    dataValues = dataOut.createVariable(ismip6_var_name, 'd',
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


def generate_output_2d_flux_vars(file_remapped_mali_flux, res_ismip6_grid,
                                 ismip6_grid_file, exp, output_path):
    """
    file_remapped_mali_flux: flux output file on mali mesh remapped
    onto the ismip6 grid
    ismip6 grid
    res_ismip6_grid: ismip6 grid resolution in kilometers
    ismip6_grid_file: ismip6 original file
    exp: ISMIP6 experiment name
    output_path: path to which the final output files are saved
    """

    print("Writing 2d flux variables")
    # ----------- acabf ------------------
    write_netcdf_2d_flux_vars('sfcMassBalApplied', 'acabf',
                              'land_ice_surface_specific_mass_balance_flux',
                              'kg m-2 s-1', 'Surface mass balance flux',
                              file_remapped_mali_flux, res_ismip6_grid,
                              ismip6_grid_file, exp, output_path)

    # ----------- libmassbffl ------------------
    write_netcdf_2d_flux_vars('libmassbffl', 'libmassbffl',
                              'land_ice_basal_specific_mass_balance_flux',
                              'kg m-2 s-1',
                              'Basal mass balance flux beneath floating ice',
                              file_remapped_mali_flux, res_ismip6_grid,
                              ismip6_grid_file, exp, output_path)

    # ----------- libmassbfgr ------------------
    write_netcdf_2d_flux_vars('libmassbfgr', 'libmassbfgr',
                              'land_ice_basal_specific_mass_balance_flux',
                              'kg m-2 s-1',
                              'Basal mass balance flux beneath grounded ice',
                              file_remapped_mali_flux, res_ismip6_grid,
                              ismip6_grid_file, exp, output_path)

    # ----------- dlithkdt ------------------
    write_netcdf_2d_flux_vars('dHdt', 'dlithkdt',
                              'tendency_of_land_ice_thickness',
                              'm s-1',
                              'Ice thickness imbalance',
                              file_remapped_mali_flux, res_ismip6_grid,
                              ismip6_grid_file, exp, output_path)

    # ----------- licalvf ------------------
    write_netcdf_2d_flux_vars('calvingFlux', 'licalvf',
                              'land_ice_specific_mass_flux_due_to_calving',
                              'kg m-2 s-1',
                              'Calving flux',
                              file_remapped_mali_flux, res_ismip6_grid,
                              ismip6_grid_file, exp, output_path)

    # ----------- lifmassbf ------------------
    # Note: facemelting and calving flux are combined above
    write_netcdf_2d_flux_vars('faceMeltAndCalvingFlux', 'lifmassbf',
                              'land_ice_specific_mass_flux_due_to_calving_and_ice_front_melting',
                              'kg m-2 s-1',
                              'Ice front melt and calving flux',
                              file_remapped_mali_flux, res_ismip6_grid,
                              ismip6_grid_file, exp, output_path)

    # ----------- ligroundf ------------------
    write_netcdf_2d_flux_vars('fluxAcrossGroundingLineOnCells', 'ligroundf',
                              'land_ice_specific_mass_flux_at_grounding_line',
                              'kg m-2 s-1',
                              'Grounding line flux',
                              file_remapped_mali_flux, res_ismip6_grid,
                              ismip6_grid_file, exp, output_path)
