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


def _write_state_var(varname, data_values, time_days, standard_name, units,
                     variable_desc, output_path, simulationStartDate, metadata):
    """Write a single 1D state (snapshot) variable to a NETCDF4_CLASSIC file."""
    nt = len(data_values)
    ds_out = Dataset(f'{output_path}/{varname}_{metadata["icesheet"]}_{metadata["group_nickname"]}_MALI_{metadata["exp"]}.nc',
                     'w', format='NETCDF4_CLASSIC')
    ds_out.createDimension('time', nt)
    var_out = ds_out.createVariable(varname, 'd', ('time',))
    time_out = ds_out.createVariable('time', 'd', ('time',))
    var_out[:] = data_values
    time_out[:] = time_days
    time_out.units = f'days since {simulationStartDate}'
    time_out.calendar = 'noleap'
    time_out.standard_name = 'time'
    time_out.long_name = 'time'
    var_out.standard_name = standard_name
    var_out.units = units
    ds_out.AUTHORS = metadata['authors']
    ds_out.MODEL = metadata['model']
    ds_out.GROUP = metadata['group']
    ds_out.VARIABLE = variable_desc
    ds_out.DATE = metadata['date']
    ds_out.close()


def _write_flux_var(varname, data_values, days_min, days_max, standard_name,
                    units, variable_desc, output_path, simulationStartDate,
                    metadata):
    """Write a single 1D flux (time-averaged) variable to a NETCDF4_CLASSIC file."""
    nt = len(data_values)
    ds_out = Dataset(f'{output_path}/{varname}_{metadata["icesheet"]}_{metadata["group_nickname"]}_MALI_{metadata["exp"]}.nc',
                     'w', format='NETCDF4_CLASSIC')
    ds_out.createDimension('time', nt)
    ds_out.createDimension('bnds', 2)
    var_out = ds_out.createVariable(varname, 'd', ('time',))
    time_out = ds_out.createVariable('time', 'd', ('time',))
    bnds_out = ds_out.createVariable('time_bnds', 'd', ('time', 'bnds'))
    var_out[:] = data_values
    time_out[:] = (days_min + days_max) / 2.0
    bnds_out[:, 0] = days_min
    bnds_out[:, 1] = days_max
    time_out.units = f'days since {simulationStartDate}'
    time_out.calendar = 'noleap'
    time_out.standard_name = 'time'
    time_out.long_name = 'time'
    var_out.standard_name = standard_name
    var_out.units = units
    ds_out.AUTHORS = metadata['authors']
    ds_out.MODEL = metadata['model']
    ds_out.GROUP = metadata['group']
    ds_out.VARIABLE = variable_desc
    ds_out.DATE = metadata['date']
    ds_out.close()


def generate_output_1d_vars(files, output_path, metadata):
    """
    Process and write 1D (scalar time-series) state and flux variables.

    Parameters
    ----------
    files : list of str
        Sorted list of globalStats.nc file paths to process.
    output_path : str
        Directory for output files.
    metadata : dict
        Submission metadata with keys: exp, icesheet, authors, group, model, date.
    """

    if output_path is None or not os.path.exists(output_path):
        output_path = os.getcwd()

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
        # Use isclose to avoid floating-point equality issues.
        ind_snap = np.where(np.isclose(decYears, years_state[i]))[0]
        if len(ind_snap) == 0:
            raise ValueError(f"No state snapshot found for year {years_state[i]}.")
        if len(ind_snap) > 1:
            print(f"WARNING: Found {len(ind_snap)} snapshots for year "
                  f"{years_state[i]}; using the first one.")
        idx_snap = ind_snap[0]

        vol_snapshot[i] = vol[idx_snap]
        vaf_snapshot[i] = vaf[idx_snap]
        gia_snapshot[i] = gia[idx_snap]
        fia_snapshot[i] = fia[idx_snap]
        days_snapshot[i] = daysSinceStart[idx_snap]

        if decYears[idx_snap] == endYr:
            break

    # this is for the flux variables
    for i in range(nt_flux):
        ind_avg = np.where(np.logical_and(decYears > years_flux[i],
                                          decYears <= (years_flux[i] + 1.0)))[0]
        if len(ind_avg) == 0:
            raise ValueError(f"No flux averaging samples found for year "
                             f"{years_flux[i]}.")
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

    # shared keyword arguments for both writer helpers
    common = dict(
        output_path=output_path,
        simulationStartDate=simulationStartDate,
        metadata=metadata,
    )

    # --- state (snapshot) variables ---
    _write_state_var('lim', vol_snapshot * 910, days_snapshot,
                     'land_ice_mass', 'kg', 'Total ice mass', **common)
    _write_state_var('limnsw', vaf_snapshot * 910, days_snapshot,
                     'land_ice_mass_not_displacing_sea_water', 'kg',
                     'Mass above floatation', **common)
    _write_state_var('iareagr', gia_snapshot, days_snapshot,
                     'grounded_ice_sheet_area', 'm2', 'Grounded ice area', **common)
    _write_state_var('iareafl', fia_snapshot, days_snapshot,
                     'floating_ice_shelf_area', 'm2', 'Floating ice area', **common)

    # --- flux (time-averaged) variables ---
    _write_flux_var('tendacabf', smb_avg / 31536000.0, days_min, days_max,
                    'tendency_of_land_ice_mass_due_to_surface_mass_balance',
                    'kg s-1', 'Total SMB flux', **common)
    _write_flux_var('tendlibmassbfgr', bmbGr_avg / 31536000.0, days_min, days_max,
                    'tendency_of_land_ice_mass_due_to_basal_mass_balance',
                    'kg s-1', 'Total BMB flux beneath grounded ice', **common)
    _write_flux_var('tendlibmassbffl', bmbFlt_avg / 31536000.0, days_min, days_max,
                    'tendency_of_land_ice_mass_due_to_basal_mass_balance',
                    'kg s-1', 'Total BMB flux beneath floating ice', **common)
    # tendlicalvf: sign convention — calving removes mass, so negate
    _write_flux_var('tendlicalvf', -cfx_avg / 31536000.0, days_min, days_max,
                    'tendency_of_land_ice_mass_due_to_calving',
                    'kg s-1', 'Total calving flux', **common)
    # tendlifmassbf: in ISMIP7 this is ice-front melting only (not calving)
    _write_flux_var('tendlifmassbf', -fmfx_avg / 31536000.0, days_min, days_max,
                    'tendency_of_land_ice_mass_due_to_ice_front_melting',
                    'kg s-1', 'Total ice front melting flux', **common)
    _write_flux_var('tendligroundf', gfx_avg / 31536000.0, days_min, days_max,
                    'tendency_of_grounded_ice_mass',
                    'kg s-1', 'Total grounding line flux', **common)

    ds.close()
