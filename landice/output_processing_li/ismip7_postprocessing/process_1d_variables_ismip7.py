"""
Functions for processing and writing 1D (scalar time-series) output variables
for ISMIP7 submissions from MALI globalStats output.
"""

from netCDF4 import Dataset, num2date, date2num, default_fillvals
import numpy as np
import os
import sys
import xarray as xr
from validate import validate_mali_files
from output_naming import build_output_filename


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

    See validate.validate_mali_files for full details of checks performed.
    """
    validate_mali_files(files, EXPECTED_VARIABLES, label='globalStats')


def _write_state_var(
        varname,
        data_values,
        time_days,
        standard_name,
        units,
        variable_desc,
        output_path,
        metadata):
    """Write one 1D state (snapshot) variable to NETCDF4_CLASSIC."""
    nt = len(data_values)
    filename = build_output_filename(output_path, varname, metadata)
    ds_out = Dataset(filename, 'w', format='NETCDF4_CLASSIC')
    ds_out.createDimension('time', None)  # unlimited dimension (record)
    var_out = ds_out.createVariable(
        varname, 'f4', ('time',), fill_value=default_fillvals['f4'])
    time_out = ds_out.createVariable(
        'time', 'f4', ('time',), fill_value=default_fillvals['f4'])
    var_out[:] = data_values
    time_out[:] = time_days
    time_out.units = 'days since 1850-01-01'
    time_out.calendar = 'standard'
    time_out.standard_name = 'time'
    time_out.long_name = 'time'
    var_out.standard_name = standard_name
    var_out.units = units
    # Write all metadata as global attributes
    for key, value in metadata.items():
        if key != 'time_range':  # time_range is added by caller
            setattr(ds_out, key, value)
    ds_out.close()


def _write_flux_var(varname, data_values, days_min, days_max, standard_name,
                    units, variable_desc, output_path, metadata):
    """Write one 1D flux (time-averaged) variable to NETCDF4_CLASSIC."""
    nt = len(data_values)
    filename = build_output_filename(output_path, varname, metadata)
    ds_out = Dataset(filename, 'w', format='NETCDF4_CLASSIC')
    ds_out.createDimension('time', None)  # unlimited dimension (record)
    ds_out.createDimension('bnds', 2)
    var_out = ds_out.createVariable(
        varname, 'f4', ('time',), fill_value=default_fillvals['f4'])
    time_out = ds_out.createVariable(
        'time', 'f4', ('time',), fill_value=default_fillvals['f4'])
    bnds_out = ds_out.createVariable(
        'time_bounds', 'f4', ('time', 'bnds'),
        fill_value=default_fillvals['f4'])
    var_out[:] = data_values
    time_out[:] = (days_min + days_max) / 2.0
    bnds_out[:, 0] = days_min
    bnds_out[:, 1] = days_max
    time_out.units = 'days since 1850-01-01'
    time_out.calendar = 'standard'
    time_out.bounds = 'time_bounds'
    time_out.standard_name = 'time'
    time_out.long_name = 'time'
    var_out.standard_name = standard_name
    var_out.units = units
    # Write all metadata as global attributes
    for key, value in metadata.items():
        if key != 'time_range':  # time_range is added by caller
            setattr(ds_out, key, value)
    ds_out.VARIABLE = variable_desc
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
        Submission metadata with keys:
        exp, icesheet, authors, group, model, date.
        time_range will be used as already computed by the main driver.
    """
    metadata = metadata.copy()

    if output_path is None or not os.path.exists(output_path):
        output_path = os.getcwd()

    ds = xr.open_mfdataset(files, combine='nested', concat_dim='Time',
                           decode_cf=False, data_vars='minimal',
                           coords='minimal', compat='override')
    with xr.open_dataset(files[0], decode_cf=False) as ds_first:
        simulationStartTime = (
            ds_first['simulationStartTime']
            .values.tobytes().decode('utf-8').strip().strip('\x00')
        )
    daysSinceStart = ds['daysSinceStart'].values
    dt = ds['deltat'].values
    simulationStartDate = simulationStartTime.split("_")[0]
    if simulationStartDate[5:10] != '01-01':
        ds.close()
        sys.exit(
            "Error: simulationStartTime for globalStats file "
            "is not on Jan. 1."
        )
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
    # the start of the year, e.g., state year 2000 means the snapshot at
    # Jan. 1, 2000.
    # For flux fields, the years is the calendar year being averaged over,
    # e.g., flux year 2000 is the average between Jan. 1, 2000,
    # and Jan. 1, 2001.
    # Note this year convention differs from the first column in table in
    # A2.3.2 at
    # https://www.climate-cryosphere.org/wiki/index.php?title=ISMIP7-Projections2300-Antarctica#A2.3.3_Table_A1:_Variable_request_for_ISMIP6
    # but that year indexing convention ultimately doesn't matter because the
    # time coordinates in these files uses units of days since a
    # reference date,
    # and it does not use a year indexing convention at all.
    if decYears[0] == np.round(decYears[0]):
        # The initial time level will only be on an even year (Jan. 1)
        # for the hist run.  In that case, we want to include that initial
        # even year in the state processing.  We also want the state snapshot
        # at the final (even) year in the output.
        # The flux processing should start with the first year, which covers a
        # full 12 months.  We exclude the final year, which is just a Jan. 1
        # posting.
        years_state = np.arange(decYears[0], endYr + 1)
        years_flux = np.arange(decYears[0], endYr)
    else:
        # For projection runs, the first state snapshot we want is the
        # first Jan. 1,
        # which we be the first even year after the initial time in the file.
        # For flux files, the first full year we want to process is the
        # year of the
        # first time level in the file. As with hist, we exclude the
        # final year,
        # which is just a Jan. 1 posting.
        years_state = np.arange(np.ceil(decYears[0]), endYr + 1)
        years_flux = np.arange(np.floor(decYears[0]), endYr)
    nt_state = len(years_state)
    nt_flux = len(years_flux)

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
        print(
            f"WARNING: Found {
                len(ind)} values of totalGroundedBasalMassBal>1.0e18")
        bmbGr[ind] = np.nan
    ind = np.nonzero(bmbGr < -1.0e18)[0]
    if len(ind) > 0:
        print(
            f"WARNING: Found {
                len(ind)} values of totalGroundedBasalMassBal<-1.0e18")
        bmbGr[ind] = np.nan
    bmbFlt = ds['totalFloatingBasalMassBal'].values.copy()
    # clean out some garbage values we can't account for
    ind = np.nonzero(bmbFlt > 1.0e18)[0]
    if len(ind) > 0:
        print(
            f"WARNING: Found {
                len(ind)} values of totalFloatingBasalMassBal>1.0e18")
        bmbFlt[ind] = np.nan
    ind = np.nonzero(bmbFlt < -1.0e18)[0]
    if len(ind) > 0:
        print(
            f"WARNING: Found {
                len(ind)} values of totalFloatingBasalMassBal<-1.0e18")
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
    # Skip any snapshot at the simulation start time (daysSinceStart ≈ 0)
    snapshot_idx = 0
    for year_idx, year in enumerate(years_state):
        # Use isclose to avoid floating-point equality issues.
        # isclose relative tolerance not meaningful here
        ind_snap = np.where(np.isclose(decYears, year, rtol=0.0))[0]
        if len(ind_snap) == 0:
            raise ValueError(
                f"No state snapshot found for year {year}.")
        if len(ind_snap) > 1:
            print(f"WARNING: Found {len(ind_snap)} snapshots for year "
                  f"{year}; using the first one. "
                  f"Snapshot values: {decYears[ind_snap]}")
        idx_snap = ind_snap[0]

        # Skip the initial snapshot at simulation start time
        if np.isclose(daysSinceStart[idx_snap], 0.0, rtol=0.0):
            print(f"Skipping state snapshot at simulation start time (year {year}).")
            continue

        vol_snapshot[snapshot_idx] = vol[idx_snap]
        vaf_snapshot[snapshot_idx] = vaf[idx_snap]
        gia_snapshot[snapshot_idx] = gia[idx_snap]
        fia_snapshot[snapshot_idx] = fia[idx_snap]
        days_snapshot[snapshot_idx] = daysSinceStart[idx_snap]
        snapshot_idx += 1

        if decYears[idx_snap] == endYr:
            break
    
    # Trim arrays to actual number of snapshots written
    vol_snapshot = vol_snapshot[:snapshot_idx]
    vaf_snapshot = vaf_snapshot[:snapshot_idx]
    gia_snapshot = gia_snapshot[:snapshot_idx]
    fia_snapshot = fia_snapshot[:snapshot_idx]
    days_snapshot = days_snapshot[:snapshot_idx]

    # this is for the flux variables
    for i in range(nt_flux):
        ind_avg = np.where(
            np.logical_and(
                decYears > years_flux[i],
                decYears <= (
                    years_flux[i] +
                    1.0)))[0]
        if len(ind_avg) == 0:
            raise ValueError(f"No flux averaging samples found for year "
                             f"{years_flux[i]}.")
        smbi = smb[ind_avg]
        bmbGri = bmbGr[ind_avg]
        bmbFlti = bmbFlt[ind_avg]
        cfxi = cfx[ind_avg]
        fmfxi = fmfx[ind_avg]
        gfxi = gfx[ind_avg]
        dti = dt[ind_avg]

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

    # Convert from MALI noleap time axis to Gregorian days since 1850-01-01
    input_units = f'days since {simulationStartDate}'
    output_units = 'days since 1850-01-01'
    days_snapshot = date2num(
        num2date(days_snapshot, units=input_units, calendar='noleap'),
        units=output_units,
        calendar='standard',
    )
    days_min = date2num(
        num2date(days_min, units=input_units, calendar='noleap'),
        units=output_units,
        calendar='standard',
    )
    days_max = date2num(
        num2date(days_max, units=input_units, calendar='noleap'),
        units=output_units,
        calendar='standard',
    )

    # shared keyword arguments for both writer helpers
    common = dict(
        output_path=output_path,
        metadata=metadata,
    )

    # --- state (snapshot) variables ---
    _write_state_var('lim', vol_snapshot * 910, days_snapshot,
                     'land_ice_mass', 'kg', 'Total ice mass', **common)
    _write_state_var('limnsw', vaf_snapshot * 910, days_snapshot,
                     'land_ice_mass_not_displacing_sea_water', 'kg',
                     'Mass above floatation', **common)
    _write_state_var(
        'iareagr',
        gia_snapshot,
        days_snapshot,
        'grounded_ice_sheet_area',
        'm^2',
        'Grounded ice area',
        **common)
    _write_state_var(
        'iareafl',
        fia_snapshot,
        days_snapshot,
        'floating_ice_shelf_area',
        'm^2',
        'Floating ice area',
        **common)

    # --- flux (time-averaged) variables ---
    _write_flux_var('tendacabf', smb_avg / 31536000.0, days_min, days_max,
                    'tendency_of_land_ice_mass_due_to_surface_mass_balance',
                    'kg s-1', 'Total SMB flux', **common)
    _write_flux_var(
        'tendlibmassbfgr',
        bmbGr_avg / 31536000.0,
        days_min,
        days_max,
        'tendency_of_land_ice_mass_due_to_basal_mass_balance',
        'kg s-1',
        'Total BMB flux beneath grounded ice',
        **common)
    _write_flux_var(
        'tendlibmassbffl',
        bmbFlt_avg / 31536000.0,
        days_min,
        days_max,
        'tendency_of_land_ice_mass_due_to_basal_mass_balance',
        'kg s-1',
        'Total BMB flux beneath floating ice',
        **common)
    # tendlicalvf: sign convention — calving removes mass, so negate
    _write_flux_var('tendlicalvf', -cfx_avg / 31536000.0, days_min, days_max,
                    'tendency_of_land_ice_mass_due_to_calving',
                    'kg s-1', 'Total calving flux', **common)
    # tendlifmassbf: in ISMIP7 this is ice-front melting only (not calving)
    _write_flux_var('tendlifmassbf', -fmfx_avg / 31536000.0, days_min,
                    days_max,
                    'tendency_of_land_ice_mass_due_to_ice_front_melting',
                    'kg s-1', 'Total ice front melting flux', **common)
    _write_flux_var('tendligroundf', gfx_avg / 31536000.0, days_min, days_max,
                    'tendency_of_grounded_ice_mass',
                    'kg s-1', 'Total grounding line flux', **common)

    ds.close()
