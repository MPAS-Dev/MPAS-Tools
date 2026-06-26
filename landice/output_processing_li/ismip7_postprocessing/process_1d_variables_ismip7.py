"""
Functions for processing and writing 1D (scalar time-series) output variables
for ISMIP7 submissions from MALI globalStats output.
"""

from netCDF4 import Dataset
from datetime import date
import numpy as np
import os
import sys
import glob
import xarray as xr


EXPECTED_VARIABLES = [
    'daysSinceStart', 'deltat', 'simulationStartTime',
    'totalIceVolume', 'volumeAboveFloatation',
    'groundedIceArea', 'floatingIceArea',
    'totalSfcMassBal', 'totalGroundedBasalMassBal',
    'totalFloatingBasalMassBal', 'totalCalvingFlux',
    'totalFaceMeltingFlux', 'groundingLineFlux',
]


def check_global_stats_files(files):
    """
    Validate a list of globalStats files before processing.

    Checks that:
    - The list is not empty
    - Each file exists
    - Each file contains the expected variables
    - simulationStartTime is consistent across all files
    - No time overlaps exist between consecutive files
    - No unexpectedly large time gaps (> 366 days) exist between consecutive files

    Parameters
    ----------
    files : list of str
        Sorted list of globalStats file paths.

    Raises
    ------
    ValueError
        If any validation check fails.
    FileNotFoundError
        If any file does not exist.
    """
    if len(files) == 0:
        raise ValueError(
            "No globalStats files matched the provided glob pattern.")

    for f in files:
        if not os.path.exists(f):
            raise FileNotFoundError(f"globalStats file not found: {f}")

    # Check required variables in each file
    for f in files:
        with xr.open_dataset(f, decode_cf=False) as ds:
            missing = [v for v in EXPECTED_VARIABLES if v not in ds]
        if missing:
            raise ValueError(
                f"File '{f}' is missing expected variables: {missing}")

    # Check simulationStartTime consistency across all files
    start_times = []
    for f in files:
        with xr.open_dataset(f, decode_cf=False) as ds:
            start_times.append(
                ds['simulationStartTime'].values.tobytes()
                .decode('utf-8').strip().strip('\x00'))
    if len(set(start_times)) > 1:
        raise ValueError(
            f"Inconsistent simulationStartTime across globalStats files: "
            f"{set(start_times)}")

    # Check for time overlaps or large gaps between consecutive files
    for i in range(len(files) - 1):
        with xr.open_dataset(files[i], decode_cf=False) as ds_a:
            end_a = float(ds_a['daysSinceStart'].values[-1])
        with xr.open_dataset(files[i + 1], decode_cf=False) as ds_b:
            start_b = float(ds_b['daysSinceStart'].values[0])
        if start_b <= end_a:
            raise ValueError(
                f"Time overlap detected between files:\n"
                f"  {files[i]} (ends at day {end_a})\n"
                f"  {files[i + 1]} (starts at day {start_b})")
        gap_days = start_b - end_a
        if gap_days > 366:
            print(f"WARNING: Gap of {gap_days:.1f} days between files:\n"
                  f"  {files[i]}\n  {files[i + 1]}")

    print(f"Validated {len(files)} globalStats file(s).")


def generate_output_1d_vars(files, exp, output_path=None):
    """
    Process and write 1D (scalar time-series) state and flux variables.

    Parameters
    ----------
    files : list of str
        Sorted list of globalStats.nc file paths to process.
    exp : str
        ISMIP7 experiment name (e.g. 'C001').
    output_path : str, optional
        Directory for output files. Defaults to current working directory.
    """

    if output_path is None or not os.path.exists(output_path):
        output_path = os.getcwd()

    AUTHOR_STR = 'Matthew Hoffman, Trevor Hillebrand, Holly Kyeore Han'
    DATE_STR = date.today().strftime("%d-%b-%Y")

    ds = xr.open_mfdataset(files, combine='nested', concat_dim='Time',
                           decode_cf=False, data_vars='minimal',
                           coords='minimal', compat='override')
    with xr.open_dataset(files[0], decode_cf=False) as ds_first:
        simulationStartTime = (ds_first['simulationStartTime'].values
                               .tobytes().decode('utf-8').strip().strip('\x00'))
    daysSinceStart = ds['daysSinceStart'].values
    dt = ds['deltat'].values
    simulationStartDate = simulationStartTime.split("_")[0]
    if simulationStartDate[5:10] != '01-01':
        ds.close()
        sys.exit("Error: simulationStartTime for globalStats file is not on Jan. 1.")
    refYear = int(simulationStartDate[0:4])
    decYears = refYear + daysSinceStart / 365.0
    endYr = decYears[-1]
    if endYr != np.round(endYr):
        ds.close()
        sys.exit("Error: end year not an even year in globalStats file.")

    # Determine processed time levels for state and flux fields
    # The historical state fields should include the initial time (Jan. 1).
    # Projection state fields should not include the initial time (Jan. 1)
    # of the projection because it's a restart from the historical.
    # Flux fields should never use the Jan. 1 time level at the start of the
    # year as part of the averaging.
    # For year conventions here, for state fields, the year is the snapshot at
    # the start of the year, e.g., state year 2000 means the snapshot at Jan. 1, 2000.
    # For flux fields, the years is the calendar year being averaged over,
    # e.g., flux year 2000 is the average between Jan. 1, 2000, and Jan. 1, 2001.
    # Note this year convention differs from the first column in table in A2.3.2 at
    # https://www.climate-cryosphere.org/wiki/index.php?title=ISMIP7-Projections2300-Antarctica#A2.3.3_Table_A1:_Variable_request_for_ISMIP6
    # but that year indexing convention ultimately doesn't matter because the
    # time coordinates in these files uses units of days since a reference date,
    # and it does not use a year indexing convention at all.
    if decYears[0] == np.round(decYears[0]):
        # The initial time level will only be on an even year (Jan. 1)
        # for the hist run.  In that case, we want to include that initial
        # even year in the state processing.  We also want the state snapshot
        # at the final (even) year in the output.
        # The flux processing should start with the first year, which covers a
        # full 12 months.  We exclude the final year, which is just a Jan. 1 posting.
        years_state = np.arange(decYears[0], endYr + 1)
        years_flux = np.arange(decYears[0], endYr)
    else:
        # For projection runs, the first state snapshot we want is the first Jan. 1,
        # which we be the first even year after the initial time in the file.
        # For flux files, the first full year we want to process is the year of the
        # first time level in the file.  As with hist, we exclude the final year,
        # which is just a Jan. 1 posting.
        years_state = np.arange(np.ceil(decYears[0]), endYr + 1)
        years_flux = np.arange(np.floor(decYears[0]), endYr)
    nt_state = len(years_state)
    nt_flux = len(years_flux)
    print(f'For state processing, using start year={years_state[0]} and end year={years_state[-1]}.')
    print(f'For flux  processing, using start year={years_flux[0]} and end year={years_flux[-1]}.')

    # read in state variables
    vol = ds['totalIceVolume'].values
    vaf = ds['volumeAboveFloatation'].values
    gia = ds['groundedIceArea'].values
    fia = ds['floatingIceArea'].values

    # read in flux variables over which yearly average will be taken
    smb = ds['totalSfcMassBal'].values
    bmbGr = ds['totalGroundedBasalMassBal'].values.copy()
    # clean out some garbage values we can't account for
    ind = np.nonzero(bmbGr > 1.0e18)[0]
    if len(ind) > 0:
        print(f"WARNING: Found {len(ind)} values of totalGroundedBasalMassBal>1.0e18")
        bmbGr[ind] = np.nan
    ind = np.nonzero(bmbGr < -1.0e18)[0]
    if len(ind) > 0:
        print(f"WARNING: Found {len(ind)} values of totalGroundedBasalMassBal<-1.0e18")
        bmbGr[ind] = np.nan
    bmbFlt = ds['totalFloatingBasalMassBal'].values.copy()
    # clean out some garbage values we can't account for
    ind = np.nonzero(bmbFlt>1.0e18)[0]
    if len(ind) > 0:
        print(f"WARNING: Found {len(ind)} values of totalFloatingBasalMassBal>1.0e18")
        bmbFlt[ind] = np.nan
    ind = np.nonzero(bmbFlt<-1.0e18)[0]
    if len(ind) > 0:
        print(f"WARNING: Found {len(ind)} values of totalFloatingBasalMassBal<-1.0e18")
        bmbFlt[ind] = np.nan
    cfx = ds['totalCalvingFlux'].values
    fmfx = ds['totalFaceMeltingFlux'].values
    gfx = ds['groundingLineFlux'].values

    # initialize 1D variables that will store data value on the
    # January 1st of each year
    vol_snapshot = np.zeros(nt_state) * np.nan
    vaf_snapshot = np.zeros(nt_state) * np.nan
    gia_snapshot = np.zeros(nt_state) * np.nan
    fia_snapshot = np.zeros(nt_state) * np.nan
    days_snapshot = np.zeros(nt_state) * np.nan
    smb_avg = np.zeros(nt_flux) * np.nan
    bmbGr_avg = np.zeros(nt_flux) * np.nan
    bmbFlt_avg = np.zeros(nt_flux) * np.nan
    cfx_avg = np.zeros(nt_flux) * np.nan
    fmfx_avg = np.zeros(nt_flux) * np.nan
    gfx_avg = np.zeros(nt_flux) * np.nan
    days_min = np.zeros(nt_flux) * np.nan
    days_max = np.zeros(nt_flux) * np.nan

    # this is for the state variables
    for i in range(nt_state):
        ind_snap = np.where(decYears==years_state[i])[0]

        vol_snapshot[i] = vol[ind_snap]
        vaf_snapshot[i] = vaf[ind_snap]
        gia_snapshot[i] = gia[ind_snap]
        fia_snapshot[i] = fia[ind_snap]
        days_snapshot[i] = daysSinceStart[ind_snap]

        if decYears[ind_snap] == endYr:
            break

    # this is for the flux variables
    for i in range(nt_flux):
        ind_avg = np.where(np.logical_and(decYears > years_flux[i],
                                          decYears <= (years_flux[i] + 1.0)))[0]
        smbi = smb[ind_avg]
        bmbGri = bmbGr[ind_avg]
        bmbFlti = bmbFlt[ind_avg]
        cfxi = cfx[ind_avg]
        fmfxi = fmfx[ind_avg]
        gfxi = gfx[ind_avg]
        dti  = dt[ind_avg]

        # take the average of the flux variables
        smb_avg[i] = np.nansum(smbi * dti) / np.nansum(dti)
        bmbGr_avg[i] = np.nansum(bmbGri * dti) / np.nansum(dti)
        bmbFlt_avg[i] = np.nansum(bmbFlti * dti) / np.nansum(dti)
        cfx_avg[i] = np.nansum(cfxi * dti) / np.nansum(dti)
        fmfx_avg[i] = np.nansum(fmfxi * dti) / np.nansum(dti)
        gfx_avg[i] = np.nansum(gfxi * dti) / np.nansum(dti)
        days_min[i] = (years_flux[i] - refYear) * 365.0
        days_max[i] = (years_flux[i] + 1.0 - refYear) * 365.0

        if decYears[ind_avg][-1] == endYr:
            break

    # -------------- lim ------------------
    data_scalar = Dataset(f'{output_path}/lim_AIS_DOE_MALI_{exp}.nc', 'w', format='NETCDF4_CLASSIC')
    data_scalar.createDimension('time', nt_state)
    limValues = data_scalar.createVariable('lim', 'd', ('time'))
    timeValues = data_scalar.createVariable('time', 'd', ('time'))
    for i in range(nt_state):
        limValues[i] = vol_snapshot[i] * 910
        timeValues[i] = days_snapshot[i]
    timeValues.units = f'days since {simulationStartDate}'
    timeValues.calendar = 'noleap'
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
    data_scalar.createDimension('time', nt_state)
    limnswValues = data_scalar.createVariable('limnsw', 'd', ('time'))
    timeValues = data_scalar.createVariable('time', 'd', ('time'))
    for i in range(nt_state):
        limnswValues[i] = vaf_snapshot[i] * 910
        timeValues[i] = days_snapshot[i]
    timeValues.units = f'days since {simulationStartDate}'
    timeValues.calendar = 'noleap'
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
    data_scalar.createDimension('time', nt_state)
    iareagrValues = data_scalar.createVariable('iareagr', 'd', ('time'))
    timeValues = data_scalar.createVariable('time', 'd', ('time'))
    for i in range(nt_state):
        iareagrValues[i] = gia_snapshot[i]
        timeValues[i] = days_snapshot[i]
    timeValues.units = f'days since {simulationStartDate}'
    timeValues.calendar = 'noleap'
    timeValues.standard_name = 'time'
    timeValues.long_name = 'time'
    iareagrValues.standard_name = 'grounded_ice_sheet_area'
    iareagrValues.units = 'm2'
    data_scalar.AUTHORS= AUTHOR_STR
    data_scalar.MODEL= 'MALI (MPAS-Albany Land Ice model)'
    data_scalar.GROUP = 'Los Alamos National Laboratory'
    data_scalar.VARIABLE = 'Grounded ice area'
    data_scalar.DATE = DATE_STR
    data_scalar.close()

    # -------------- iareafl ------------------
    data_scalar = Dataset(f'{output_path}/iareafl_AIS_DOE_MALI_{exp}.nc', 'w', format='NETCDF4_CLASSIC')
    data_scalar.createDimension('time', nt_state)
    iareaflValues = data_scalar.createVariable('iareafl', 'd', ('time'))
    timeValues = data_scalar.createVariable('time', 'd', ('time'))
    for i in range(nt_state):
        iareaflValues[i] = fia_snapshot[i]
        timeValues[i] = days_snapshot[i]
    timeValues.units = f'days since {simulationStartDate}'
    timeValues.calendar = 'noleap'
    timeValues.standard_name = 'time'
    timeValues.long_name = 'time'
    iareaflValues.standard_name = 'floating_ice_shelf_area'
    iareaflValues.units = 'm2'
    data_scalar.AUTHORS= AUTHOR_STR
    data_scalar.MODEL= 'MALI (MPAS-Albany Land Ice model)'
    data_scalar.GROUP = 'Los Alamos National Laboratory'
    data_scalar.VARIABLE = 'Floating ice area'
    data_scalar.DATE = DATE_STR
    data_scalar.close()

    # -------------- tendacabf: this is a flux var
    data_scalar = Dataset(f'{output_path}/tendacabf_AIS_DOE_MALI_{exp}.nc', 'w', format='NETCDF4_CLASSIC')
    data_scalar.createDimension('time', nt_flux)
    tendacabfValues = data_scalar.createVariable('tendacabf', 'd', ('time'))
    timeValues = data_scalar.createVariable('time', 'd', ('time'))
    data_scalar.createDimension('bnds', 2)
    timebndsValues = data_scalar.createVariable('time_bnds', 'd', ('time', 'bnds'))
    for i in range(nt_flux):
        tendacabfValues[i] = smb_avg[i] / 31536000.0
        timeValues[i] = (days_min[i] + days_max[i]) / 2.0
        timebndsValues[i, 0] = days_min[i]
        timebndsValues[i, 1] = days_max[i]
    timeValues.units = f'days since {simulationStartDate}'
    timeValues.calendar = 'noleap'
    timeValues.standard_name = 'time'
    timeValues.long_name = 'time'
    tendacabfValues.standard_name = 'tendency_of_land_ice_mass_due_to_surface_mass_balance'
    tendacabfValues.units = 'kg s-1'
    data_scalar.AUTHORS= AUTHOR_STR
    data_scalar.MODEL= 'MALI (MPAS-Albany Land Ice model)'
    data_scalar.GROUP = 'Los Alamos National Laboratory'
    data_scalar.VARIABLE = 'Total SMB flux'
    data_scalar.DATE = DATE_STR
    data_scalar.close()

    # -------------- tendlibmassbfgr: this is a flux var
    data_scalar = Dataset(f'{output_path}/tendlibmassbfgr_AIS_DOE_MALI_{exp}.nc', 'w', format='NETCDF4_CLASSIC')
    data_scalar.createDimension('time', nt_flux)
    tendlibmassbfgrValues = data_scalar.createVariable('tendlibmassbfgr', 'd', ('time'))
    timeValues = data_scalar.createVariable('time', 'd', ('time'))
    data_scalar.createDimension('bnds', 2)
    timebndsValues = data_scalar.createVariable('time_bnds', 'd', ('time', 'bnds'))
    for i in range(nt_flux):
        tendlibmassbfgrValues[i] = bmbGr_avg[i] / 31536000.0
        timeValues[i] = (days_min[i] + days_max[i]) / 2.0
        timebndsValues[i, 0] = days_min[i]
        timebndsValues[i, 1] = days_max[i]
    timeValues.units = f'days since {simulationStartDate}'
    timeValues.calendar = 'noleap'
    timeValues.standard_name = 'time'
    timeValues.long_name = 'time'
    tendlibmassbfgrValues.standard_name = 'tendency_of_land_ice_mass_due_to_basal_mass_balance'
    tendlibmassbfgrValues.units = 'kg s-1'
    data_scalar.AUTHORS= AUTHOR_STR
    data_scalar.MODEL= 'MALI (MPAS-Albany Land Ice model)'
    data_scalar.GROUP = 'Los Alamos National Laboratory'
    data_scalar.VARIABLE = 'Total BMB flux beneath grounded ice'
    data_scalar.DATE = DATE_STR
    data_scalar.close()

    # -------------- tendlibmassbffl: this is a flux var
    data_scalar = Dataset(f'{output_path}/tendlibmassbffl_AIS_DOE_MALI_{exp}.nc', 'w',
                          format='NETCDF4_CLASSIC')
    data_scalar.createDimension('time', nt_flux)
    tendlibmassbfflValues = data_scalar.createVariable('tendlibmassbffl', 'd', ('time'))
    timeValues = data_scalar.createVariable('time', 'd', ('time'))
    data_scalar.createDimension('bnds', 2)
    timebndsValues = data_scalar.createVariable('time_bnds', 'd', ('time', 'bnds'))
    for i in range(nt_flux):
        tendlibmassbfflValues[i] = bmbFlt_avg[i] / 31536000
        timeValues[i] = (days_min[i] + days_max[i]) / 2.0
        timebndsValues[i, 0] = days_min[i]
        timebndsValues[i, 1] = days_max[i]
    timeValues.units = f'days since {simulationStartDate}'
    timeValues.calendar = 'noleap'
    timeValues.standard_name = 'time'
    timeValues.long_name = 'time'
    tendlibmassbfflValues.standard_name = 'tendency_of_land_ice_mass_due_to_basal_mass_balance'
    tendlibmassbfflValues.units = 'kg s-1'
    data_scalar.AUTHORS= AUTHOR_STR
    data_scalar.MODEL= 'MALI (MPAS-Albany Land Ice model)'
    data_scalar.GROUP = 'Los Alamos National Laboratory'
    data_scalar.VARIABLE = 'Total BMB flux beneath floating ice'
    data_scalar.DATE = DATE_STR
    data_scalar.close()

    # -------------- tendlicalvf: this is a flux var
    data_scalar = Dataset(f'{output_path}/tendlicalvf_AIS_DOE_MALI_{exp}.nc', 'w', format='NETCDF4_CLASSIC')
    data_scalar.createDimension('time', nt_flux)
    tendlicalvfValues = data_scalar.createVariable('tendlicalvf', 'd', ('time'))
    timeValues = data_scalar.createVariable('time', 'd', ('time'))
    data_scalar.createDimension('bnds', 2)
    timebndsValues = data_scalar.createVariable('time_bnds', 'd', ('time', 'bnds'))
    for i in range(nt_flux):
        tendlicalvfValues[i] = -cfx_avg[i] / 31536000
        timeValues[i] = (days_min[i] + days_max[i]) / 2.0
        timebndsValues[i, 0] = days_min[i]
        timebndsValues[i, 1] = days_max[i]
    timeValues.units = f'days since {simulationStartDate}'
    timeValues.calendar = 'noleap'
    timeValues.standard_name = 'time'
    timeValues.long_name = 'time'
    tendlicalvfValues.standard_name = 'tendency_of_land_ice_mass_due_to_calving'
    tendlicalvfValues.units = 'kg s-1'
    data_scalar.AUTHORS= AUTHOR_STR
    data_scalar.MODEL= 'MALI (MPAS-Albany Land Ice model)'
    data_scalar.GROUP = 'Los Alamos National Laboratory'
    data_scalar.VARIABLE = 'Total calving flux'
    data_scalar.DATE = DATE_STR
    data_scalar.close()

    # -------------- tendlifmassbf: this is a flux var
    # In ISMIP6, this variable used to be 'Total calving and ice front melting flux'
    # In ISMIP7, it represents 'Total ice front melting flux' only, without calving flux
    data_scalar = Dataset(f'{output_path}/tendlifmassbf_AIS_DOE_MALI_{exp}.nc', 'w', format='NETCDF4_CLASSIC')
    data_scalar.createDimension('time', nt_flux)
    tendlifmassbfValues = data_scalar.createVariable('tendlifmassbf', 'd', ('time'))
    timeValues = data_scalar.createVariable('time', 'd', ('time'))
    data_scalar.createDimension('bnds', 2)
    timebndsValues = data_scalar.createVariable('time_bnds', 'd', ('time', 'bnds'))
    for i in range(nt_flux):
        tendlifmassbfValues[i] = -fmfx_avg[i] / 31536000
        timeValues[i] = (days_min[i] + days_max[i]) / 2.0
        timebndsValues[i, 0] = days_min[i]
        timebndsValues[i, 1] = days_max[i]
    timeValues.units = f'days since {simulationStartDate}'
    timeValues.calendar = 'noleap'
    timeValues.standard_name = 'time'
    timeValues.long_name = 'time'
    tendlifmassbfValues.standard_name = 'tendency_of_land_ice_mass_due_to_ice_front_melting'
    tendlifmassbfValues.units = 'kg s-1'
    data_scalar.AUTHORS= AUTHOR_STR
    data_scalar.MODEL= 'MALI (MPAS-Albany Land Ice model)'
    data_scalar.GROUP = 'Los Alamos National Laboratory'
    data_scalar.VARIABLE = 'Total ice front melting flux'
    data_scalar.DATE = DATE_STR
    data_scalar.close()

    # -------------- tendligroundf: this is a flux var
    data_scalar = Dataset(f'{output_path}/tendligroundf_AIS_DOE_MALI_{exp}.nc', 'w', format='NETCDF4_CLASSIC')
    data_scalar.createDimension('time', nt_flux)
    tendligroundfValues = data_scalar.createVariable('tendligroundf', 'd', ('time'))
    timeValues = data_scalar.createVariable('time', 'd', ('time'))
    data_scalar.createDimension('bnds', 2)
    timebndsValues = data_scalar.createVariable('time_bnds', 'd', ('time', 'bnds'))
    for i in range(nt_flux):
        tendligroundfValues[i] = gfx_avg[i] / 31536000
        timeValues[i] = (days_min[i] + days_max[i]) / 2.0
        timebndsValues[i, 0] = days_min[i]
        timebndsValues[i, 1] = days_max[i]
    timeValues.units = f'days since {simulationStartDate}'
    timeValues.calendar = 'noleap'
    timeValues.standard_name = 'time'
    timeValues.long_name = 'time'
    tendligroundfValues.standard_name = 'tendency_of_grounded_ice_mass'
    tendligroundfValues.units = 'kg s-1'
    data_scalar.AUTHORS= AUTHOR_STR
    data_scalar.MODEL= 'MALI (MPAS-Albany Land Ice model)'
    data_scalar.GROUP = 'Los Alamos National Laboratory'
    data_scalar.VARIABLE = 'Total grounding line flux'
    data_scalar.DATE = DATE_STR
    data_scalar.close()

    ds.close()
