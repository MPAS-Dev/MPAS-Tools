"""
This script has functions that are needed to post-process and write flux
output variables from ISMIP7 simulations.
"""

from netCDF4 import Dataset, num2date, date2num, default_fillvals
import xarray as xr
import numpy as np
from subprocess import check_call
import os
from validate import validate_mali_files
from gl_flux_from_restart_ismip7 import gl_flux_from_restarts
from output_naming import build_output_filename
from netcdf_utils import to_masked, zero_where_invalid
import datetime
from time_utils import (get_calendar, decimal_years, days_at_jan1,
                        days_to_ismip7_time)


EXPECTED_FLUX_VARIABLES = [
    'daysSinceStart', 'simulationStartTime', 'cellMask',
    'avgSMBFlux', 'avgFloatingBMBFlux', 'avgGroundedBMBFlux',
    'avgDhdt', 'avgCalvingFlux', 'avgFaceMeltFlux',
    'avgGroundingLineFlux',
]

REQUIRED_FLUX_REMAPPING_VARIABLES = [
    'daysSinceStart', 'simulationStartTime',
    'avgSMBFlux', 'avgFloatingBMBFlux', 'avgGroundedBMBFlux',
    'avgDhdt', 'avgCalvingFlux', 'avgFaceMeltFlux',
    'avgGroundingLineFlux',
    # sftgrf/sftflf are used only to set the basal mass balance fluxes to 0
    # where there is no grounded/floating ice
    'sftgrf', 'sftflf',
    # masks not currently used - see note below in process_flux_vars
    #'ice_mask', 'floating_mask', 'dynamic_mask', 'grounded_mask',
    #'gl_mask'
]


def check_flux_files(files):
    """
    Validate a list of MALI flux output files before processing.

    See validate.validate_mali_files for full details of checks performed.
    """
    validate_mali_files(files, EXPECTED_FLUX_VARIABLES, label='flux')


def process_flux_vars(files, tmp_file, input_file_mesh=None,
                      gl_restart_files=None):
    """
    Prepare flux files into a temporary file for remapping.
    This is very simple, because flux variables are already
    time-averaged by MALI.

    Parameters
    ----------
    files : list of str
        Sorted list of MALI flux output file paths.
    tmp_file : str
        Temporary output file name.
    input_file_mesh : str, optional
        MALI file with the time-invariant basalHeatFlux, used for hfgeoubed.
    gl_restart_files : list of str, optional
        If given, avgGroundingLineFlux is replaced by an approximation
        computed from restart files (for runs whose
        native avgGroundingLineFlux is unreliable).
    """
    ds_flux = xr.open_mfdataset(files, combine='nested',
                                concat_dim='Time',
                                engine='netcdf4',
                                decode_cf=False,
                                data_vars='minimal',
                                coords='minimal',
                                compat='override')
    if 'units' in ds_flux['daysSinceStart'].attrs:
        del ds_flux['daysSinceStart'].attrs['units']

    # variables that need to be modified prior to remapping
    # these masks are not currently used, but would be needed if the
    # flux vars are meant to be masked.
    # However there is a potential issue that fluxes could have occurred
    # in locations that fall outside of the final mask, because
    # the fluxes are time-averaged over a year, and the mask is only
    # valid at the end of the year.
    #ds_flux['ice_mask'] = (ds_flux['cellMask'][:, :] & 2) / 2  # grounded: dynamic ice
    #ds_flux['floating_mask'] = (ds_flux['cellMask'][:, :] & 4) / 4
    #ds_flux['dynamic_mask'] = (ds_flux['cellMask'][:, :] & 2) / 2
    #ds_flux['grounded_mask'] = (ds_flux['cellMask'][:, :] * 0 + 1) - ds_flux['floating_mask']
    #ds_flux['gl_mask'] = (ds_flux['cellMask'][:, :] & 256) / 256

    # same definitions as in the state processing
    cell_mask = ds_flux['cellMask'][:, :]
    floating = (cell_mask & 4) / 4
    dynamic = (cell_mask & 2) / 2
    ds_flux['sftflf'] = floating * dynamic
    ds_flux['sftgrf'] = (1.0 - floating) * dynamic

    ds_flux['avgCalvingFlux'] = -1.0 * ds_flux['avgCalvingFlux']  # flip sign of calving flux to match ISMIP7 convention
    ds_flux['avgFaceMeltFlux'] = -1.0 * ds_flux['avgFaceMeltFlux']  # flip sign of face melt flux to match ISMIP7 convention
    if gl_restart_files is not None:
        if input_file_mesh is None:
            raise ValueError("--input_mesh is required to approximate "
                             "ligroundf from restart files.")
        print("Replacing avgGroundingLineFlux with an approximation from "
              "restart snapshots.")
        gl = gl_flux_from_restarts(gl_restart_files, input_file_mesh,
                                   ds_flux['daysSinceStart'].values)
        print("Approximate grounding-line flux min, max (kg m-2 s-1): ",
              np.nanmin(gl), np.nanmax(gl))
        ds_flux['avgGroundingLineFlux'] = xr.DataArray(
            gl, dims=('Time', 'nCells'))

    # MALI reports the grounding-line flux as positive from grounded to
    # floating ice; ISMIP7 uses a mass-change convention (like calving), so
    # loss from the grounded ice is negative.  Negative MALI values (flow
    # onto ice rumples) become small positive values and are kept.
    ds_flux['avgGroundingLineFlux'] = -1.0 * ds_flux['avgGroundingLineFlux']
    ds_flux['avgDhdt'] = ds_flux['avgDhdt'] / (3600.0 * 24.0 * 365.0)  # convert from m/yr to m/s

    # hfgeoubed is constant in time, so no averaging is needed; broadcast it
    # to the flux times so it follows the FL time convention like the fluxes.
    have_hfgeo = False
    if input_file_mesh is not None:
        with xr.open_dataset(input_file_mesh, decode_cf=False) as ds_mesh:
            if 'basalHeatFlux' in ds_mesh:
                hfgeo = ds_mesh['basalHeatFlux'].isel(Time=0).values
                ntime = ds_flux.sizes['Time']
                ds_flux['hfgeoubed'] = xr.DataArray(
                    np.broadcast_to(hfgeo, (ntime,) + hfgeo.shape).copy(),
                    dims=('Time', 'nCells'), attrs={'units': 'W m-2'})
                have_hfgeo = True
            else:
                print(f"basalHeatFlux not found in {input_file_mesh}; "
                      "hfgeoubed will not be written.")

    missing = [
        var for var in REQUIRED_FLUX_REMAPPING_VARIABLES
        if var not in ds_flux
    ]
    if missing:
        raise ValueError(
            "Processed flux dataset is missing required output variables: "
            f"{missing}"
        )

    # subset to only required variables to keep file small for remapping
    keep = list(REQUIRED_FLUX_REMAPPING_VARIABLES)
    if have_hfgeo:
        keep.append('hfgeoubed')
    ds_flux_out = ds_flux[keep]
    # remove the first time step (which is always 0)
    time_mask = ~np.isclose(ds_flux_out['daysSinceStart'].values, 0.0, rtol=0.0)
    if not np.any(time_mask):
        raise ValueError("No flux time records remain after dropping daysSinceStart==0.")
    ds_flux_out = ds_flux_out.isel(Time=time_mask)

    ds_flux_out.to_netcdf(tmp_file)

    ds_flux_out.close()
    ds_flux.close()


def write_netcdf_2d_flux_vars(mali_var_name, ismip7_var_name, var_std_name,
                              var_units, var_varname, remapped_mali_flux_file,
                              ismip7_grid_file, output_path, metadata,
                              calendar='gregorian'):
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
    simulationStartTime = data.variables['simulationStartTime'][:].tobytes(
    ).decode('utf-8').strip().strip('\x00')
    simulationStartDate = simulationStartTime.split("_")[0]
    daysSinceStart = data.variables['daysSinceStart'][:]
    decYears = decimal_years(daysSinceStart, simulationStartDate, calendar)
    years_flux = np.floor(decYears + 1.0e-6)

    # Flux outputs are annual means over each calendar year.
    # Bounds are Jan 1 of the previous year through Jan 1 of the year indexed.
    timeBndsMin = days_at_jan1(years_flux - 1, simulationStartDate, calendar)
    timeBndsMax = days_at_jan1(years_flux, simulationStartDate, calendar)
    # Convert from the MALI model-calendar time axis to Gregorian days since
    # 1850-01-01
    timeBndsMin = days_to_ismip7_time(timeBndsMin, simulationStartDate,
                                      calendar)
    timeBndsMax = days_to_ismip7_time(timeBndsMax, simulationStartDate,
                                      calendar)
    if mali_var_name not in data.variables:
        print(f"WARNING: {mali_var_name} not present.  Skipping.")
        data.close()
        return
    var_mali = data.variables[mali_var_name][:, :, :]
    var_mali[np.where(abs(var_mali + 1e34) < 1e33)] = np.nan
    # Basal mass balance fluxes must be exactly 0 where the corresponding area
    # fraction (as written to the sftgrf/sftflf files, in float32) is 0.
    bmb_mask_name = {'libmassbfgr': 'sftgrf',
                     'libmassbffl': 'sftflf'}.get(ismip7_var_name)
    if bmb_mask_name is not None:
        bmb_mask = np.asarray(data.variables[bmb_mask_name][:, :, :],
                              dtype=np.float32)
        bmb_mask = np.where(np.abs(bmb_mask) > 1e30, np.float32(0.0),
                            bmb_mask)
    timeSteps, latN, lonN = np.shape(var_mali)

    output_filename = build_output_filename(
        output_path,
        ismip7_var_name,
        metadata,
    )
    dataOut = Dataset(output_filename, 'w', format='NETCDF4_CLASSIC')
    dataOut.createDimension('time', None)  # unlimited dimension (record)
    dataOut.createDimension('bnds', 2)
    timebndsValues = dataOut.createVariable(
        'time_bounds', 'f4', ('time', 'bnds'),
        fill_value=default_fillvals['f4'])
    dataOut.createDimension('x', lonN)
    dataOut.createDimension('y', latN)
    dataValues = dataOut.createVariable(ismip7_var_name, 'f4',
                                        ('time', 'y', 'x'),
                                        fill_value=default_fillvals['f4'])
    xValues = dataOut.createVariable('x', 'f4', ('x',),
                                     fill_value=default_fillvals['f4'])
    yValues = dataOut.createVariable('y', 'f4', ('y',),
                                     fill_value=default_fillvals['f4'])
    timeValues = dataOut.createVariable('time', 'f4', ('time',),
                                        fill_value=default_fillvals['f4'])

    for i in range(timeSteps):
        # mask not currently used - see note above in process_flux_vars
        #mask = data.variables['ice_mask'][i, :, :]
        tmp = var_mali[i, :, :]
        #tmp[mask == 0] = np.nan
        if bmb_mask_name is not None:
            tmp = np.where(bmb_mask[i, :, :] > 0, tmp, 0.0)
        # all ISMIP7 2D flux variables must be 0 (not missing) where there is
        # no ice
        if ismip7_var_name == 'hfgeoubed':
            # defined throughout the domain, missing only outside it
            dataValues[i, :, :] = to_masked(tmp)
        else:
            dataValues[i, :, :] = zero_where_invalid(tmp)
        # ISMIP7 convention for flux (FL) variables: the time coordinate is
        # Jul 1 of the nominal (averaging) year, which is the year before the
        # Jan 1 record that holds the annual mean.  The averaging interval is
        # given by time_bounds.
        nominal_year = int(years_flux[i]) - 1
        timeValues[i] = date2num(datetime.datetime(nominal_year, 7, 1),
                                 units='days since 1850-01-01',
                                 calendar='standard')
        timebndsValues[i, 0] = timeBndsMin[i]
        timebndsValues[i, 1] = timeBndsMax[i]

    xValues[:] = var_x
    yValues[:] = var_y

    dataValues.standard_name = var_std_name
    dataValues.units = var_units
    timeValues.bounds = 'time_bounds'
    timeValues.units = 'days since 1850-01-01'
    timeValues.calendar = 'standard'
    timeValues.standard_name = 'time'
    timeValues.long_name = 'time'
    timebndsValues.units = 'days since 1850-01-01'
    timebndsValues.calendar = 'standard'
    xValues.units = 'm'
    xValues.standard_name = 'x'
    xValues.long_name = 'x'
    yValues.units = 'm'
    yValues.standard_name = 'y'
    yValues.long_name = 'y'
    # Write all metadata as global attributes
    for key, value in metadata.items():
        if key != 'time_range':  # time_range is added by caller
            setattr(dataOut, key, value)
    dataOut.VARIABLE = var_varname
    dataOut.close()
    data.close()


def generate_output_2d_flux_vars(file_remapped_mali_flux,
                                 ismip7_grid_file, output_path, metadata,
                                 calendar='gregorian'):
    """
    file_remapped_mali_flux: flux output file on mali mesh remapped
    onto the ismip7 grid
    ismip7 grid
    ismip7_grid_file: ismip7 original file
    exp: ISMIP7 experiment name
    output_path: path to which the final output files are saved
    """

    # ----------- hfgeoubed ------------------
    with Dataset(file_remapped_mali_flux, 'r') as ds_check:
        have_hfgeo = 'hfgeoubed' in ds_check.variables
    if have_hfgeo:
        print("* Writing hfgeoubed")
        write_netcdf_2d_flux_vars('hfgeoubed', 'hfgeoubed',
                                  'upward_geothermal_heat_flux_in_land_ice',
                                  'W m-2', 'Geothermal heat flux',
                                  file_remapped_mali_flux,
                                  ismip7_grid_file, output_path, metadata,
                                  calendar=calendar)
    else:
        print("* Skipping hfgeoubed (no basalHeatFlux from input mesh file)")

    # ----------- acabf ------------------
    print("* Writing acabf")
    write_netcdf_2d_flux_vars('avgSMBFlux', 'acabf',
                              'land_ice_surface_specific_mass_balance_flux',
                              'kg m-2 s-1', 'Surface mass balance flux',
                              file_remapped_mali_flux,
                              ismip7_grid_file, output_path, metadata,
                              calendar=calendar)

    # ----------- libmassbffl ------------------
    print("* Writing libmassbffl")
    write_netcdf_2d_flux_vars('avgFloatingBMBFlux', 'libmassbffl',
                              'land_ice_basal_specific_mass_balance_flux',
                              'kg m-2 s-1',
                              'Basal mass balance flux beneath floating ice',
                              file_remapped_mali_flux,
                              ismip7_grid_file, output_path, metadata,
                              calendar=calendar)

    # ----------- libmassbfgr ------------------
    print("* Writing libmassbfgr")
    write_netcdf_2d_flux_vars('avgGroundedBMBFlux', 'libmassbfgr',
                              'land_ice_basal_specific_mass_balance_flux',
                              'kg m-2 s-1',
                              'Basal mass balance flux beneath grounded ice',
                              file_remapped_mali_flux,
                              ismip7_grid_file, output_path, metadata,
                              calendar=calendar)

    # ----------- dlithkdt ------------------
    write_netcdf_2d_flux_vars('avgDhdt', 'dlithkdt',
                              'tendency_of_land_ice_thickness',
                              'm s-1',
                              'Ice thickness imbalance',
                              file_remapped_mali_flux,
                              ismip7_grid_file, output_path, metadata,
                              calendar=calendar)

    # ----------- licalvf ------------------
    print("* Writing licalvf")
    write_netcdf_2d_flux_vars('avgCalvingFlux', 'licalvf',
                              'land_ice_specific_mass_flux_due_to_calving',
                              'kg m-2 s-1',
                              'Calving flux',
                              file_remapped_mali_flux,
                              ismip7_grid_file, output_path, metadata,
                              calendar=calendar)

    # ----------- lifmassbf ------------------
    print("* Writing lifmassbf")
    # Note: facemelting and calving flux are combined above
    write_netcdf_2d_flux_vars('avgFaceMeltFlux', 'lifmassbf',
                              'TBD by ISMIP7',
                              'kg m-2 s-1',
                              'Ice front melt flux',
                              file_remapped_mali_flux,
                              ismip7_grid_file, output_path, metadata,
                              calendar=calendar)

    # ----------- ligroundf ------------------
    print("* Writing ligroundf")
    write_netcdf_2d_flux_vars('avgGroundingLineFlux', 'ligroundf',
                              'TBD by ISMIP7',
                              'kg m-2 s-1',
                              'Grounding line flux',
                              file_remapped_mali_flux,
                              ismip7_grid_file, output_path, metadata,
                              calendar=calendar)


def process_flux_pipeline(flux_files, mapping_file, ismip7_grid_file,
                          output_path, metadata, input_file_mesh=None,
                          gl_restart_files=None):
    """
    Full flux-variable processing pipeline:
    validate, concatenate, remap, write.

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
        Submission metadata dict with time_range computed from globalStats.
    input_file_mesh : str, optional
        MALI file with basalHeatFlux, used for hfgeoubed.
    gl_restart_files : list of str, optional
        Restart files; if given, ligroundf is approximated from them
        instead of using the native avgGroundingLineFlux.
    """
    metadata = metadata.copy()

    check_flux_files(flux_files)

    with xr.open_dataset(flux_files[0], decode_cf=False) as ds_first:
        calendar = get_calendar(ds_first.attrs)
    print(f"Model calendar: {calendar}")

    print("Preparing concatenated flux file for remapping.")
    tmp_flux_file = 'tmp_flux.nc'
    process_flux_vars(flux_files, tmp_flux_file, input_file_mesh,
                      gl_restart_files)

    print("Remapping flux file.")
    remapped_file_flux = 'remapped_flux.nc'
    check_call(["ncremap",
                "-i", tmp_flux_file,
                "-o", remapped_file_flux,
                "-m", mapping_file,
                "-P", "mpas"])

    print("Writing processed and remapped flux fields to ISMIP7 file format.")
    generate_output_2d_flux_vars(remapped_file_flux,
                                 ismip7_grid_file,
                                 output_path, metadata, calendar=calendar)

    #os.remove(tmp_flux_file)
    #os.remove(remapped_file_flux)
