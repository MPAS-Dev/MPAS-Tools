#!/usr/bin/env python
"""
Interpolate dh/dt and mean surface elevation from an ITS_LIVE ice
elevation change dataset to a MALI mesh over a user-specified time window.

The ITS_LIVE dataset contains cumulative surface elevation change (dh)
and absolute surface elevation (h) relative to a reference epoch. This
script selects the time window between start_date and end_date, computes
dh/dt (m/s) from the endpoints, computes the time-mean h over the window,
and remaps both fields to the destination MALI mesh using conservative
interpolation.

The interpolated fields appended to the MALI mesh are:
  - observedThicknessTendency: dh/dt in m/s
  - observedThicknessTendencyUncertainty: propagated RMS uncertainty in m/s
  - observedSurfaceElevation: time-mean surface elevation in m
  - observedSurfaceElevationUncertainty: mean RMS uncertainty in m

Note that this should be run on a compute node unless the weights file already exists.
"""

import argparse
import datetime
import os
import subprocess
import sys

import numpy as np
import pandas as pd
import xarray

from mpas_tools.io import write_netcdf
from mpas_tools.scrip.from_mpas import scrip_from_mpas


def mask_source_scrip(source_scrip, masked_scrip, data_file):
    """
    Create a masked copy of the source SCRIP file with grid_imask set to 0
    for cells that contain invalid data (NaN/fill) in any field.

    Parameters
    ----------
    source_scrip : str
        Input source SCRIP file.
    masked_scrip : str
        Output masked SCRIP file path.
    data_file : str
        Path to the gridded data file (itslive_dhdt.nc) used to determine
        which cells have valid data.
    """
    # Determine data-validity mask from the source data file
    with xarray.open_dataset(data_file) as ds_data:
        # A cell is valid only if ALL fields have finite values
        valid_data = np.ones(ds_data.dims['y'] * ds_data.dims['x'], dtype=bool)
        for var in ds_data.data_vars:
            vals = ds_data[var].values
            if vals.ndim == 3:
                vals = vals[0]  # remove Time dimension
            valid_data &= np.isfinite(vals).ravel()

    with xarray.open_dataset(source_scrip) as ds_src:
        n_total = valid_data.size
        n_active = int(valid_data.sum())
        n_invalid = n_total - n_active
        print(f'  Source SCRIP masking: {n_active}/{n_total} cells active '
              f'({n_invalid} masked for invalid data)')

        ds_out = ds_src.copy()
        ds_out['grid_imask'] = xarray.DataArray(
            valid_data.astype(np.int32),
            dims=('grid_size',),
            attrs={'long_name': '0/1 mask for active source cells'})
        ds_out.to_netcdf(masked_scrip, mode='w')


def run_command(args):
    """Run a subprocess command and raise on failure."""
    print(f"  Running: {' '.join(args)}")
    subprocess.check_call(args)


def interpolate_itslive_dhdt(itslive_file, mali_file, start_date, end_date,
                             proj='gis-gimp', method='conserve',
                             mali_scrip_file=None,
                             source_scrip_file=None, weights_file=None):
    """
    Compute dh/dt and mean surface elevation from ITS_LIVE and interpolate
    to a MALI mesh.

    Parameters
    ----------
    itslive_file : str
        Path to the ITS_LIVE ice elevation change netCDF file.

    mali_file : str
        Path to the MALI mesh netCDF file. Interpolated fields will be
        appended to this file.

    start_date : str
        Start date for the time window (format: YYYY-MM-DD).

    end_date : str
        End date for the time window (format: YYYY-MM-DD).

    proj : str, optional
        Projection of the source dataset (default: 'gis-gimp').

    method : str, optional
        Interpolation method for ESMF_RegridWeightGen. Options include
        'conserve' and 'bilinear' (default: 'conserve').

    mali_scrip_file : str or None, optional
        Path to an existing SCRIP file for the MALI mesh. If None, one
        will be generated.

    source_scrip_file : str or None, optional
        Path to an existing SCRIP file for the ITS_LIVE grid. If None,
        one will be generated.

    weights_file : str or None, optional
        Path to an existing ESMF remapping weights file. If None, one
        will be generated.
    """

    print(f'Computing dh/dt and mean h from ITS_LIVE dataset between '
          f'{start_date} and {end_date}')

    ds = xarray.open_dataset(itslive_file)

    # Convert user dates to the dataset's time coordinate
    t_start = pd.Timestamp(start_date)
    t_end = pd.Timestamp(end_date)

    # Select nearest time steps for dh/dt endpoints
    h_start = ds['h'].sel(time=t_start, method='nearest')
    h_end = ds['h'].sel(time=t_end, method='nearest')

    # Get actual selected times for accurate dt calculation
    t0 = pd.Timestamp(h_start.time.values)
    t1 = pd.Timestamp(h_end.time.values)
    dt_seconds = (t1 - t0).total_seconds()
    dt_years = dt_seconds / (365.25 * 24 * 3600)

    print(f'  Nearest start time: {t0.strftime("%Y-%m-%d")}')
    print(f'  Nearest end time:   {t1.strftime("%Y-%m-%d")}')
    print(f'  dt = {dt_years:.4f} years')

    if dt_years <= 0:
        raise ValueError(
            f'end_date ({end_date}) must be after start_date ({start_date}). '
            f'Selected time steps: {t0} and {t1}')

    # Compute dh/dt in m/s
    dhdt = (h_end.values - h_start.values) / dt_seconds

    # Compute dh/dt uncertainty by propagating RMS errors from endpoints
    rms_start = ds['rms'].sel(time=t_start, method='nearest').values
    rms_end = ds['rms'].sel(time=t_end, method='nearest').values
    dhdt_err = np.sqrt(rms_start**2 + rms_end**2) / dt_seconds

    # Compute time-mean surface elevation over the specified window
    h_window = ds['h'].sel(time=slice(t_start, t_end))
    h_mean = h_window.mean(dim='time').values
    n_times = h_window.sizes['time']
    print(f'  Averaged h over {n_times} time steps')

    # Compute surface elevation uncertainty as mean RMS over the window
    rms_window = ds['rms'].sel(time=slice(t_start, t_end))
    h_err = rms_window.mean(dim='time').values

    # Write all fields to a temporary gridded netCDF file
    dhdt_filename = 'itslive_dhdt.nc'
    ds_out = xarray.Dataset(
        {'observedThicknessTendency':
         (['Time', 'y', 'x'], dhdt[np.newaxis].astype(np.float64)),
         'observedThicknessTendencyUncertainty':
         (['Time', 'y', 'x'], dhdt_err[np.newaxis].astype(np.float64)),
         'observedSurfaceElevation':
         (['Time', 'y', 'x'], h_mean[np.newaxis].astype(np.float64)),
         'observedSurfaceElevationUncertainty':
         (['Time', 'y', 'x'], h_err[np.newaxis].astype(np.float64))},
        coords={'x': ds['x'], 'y': ds['y'], 'Time': [0.0]})
    ds_out['observedThicknessTendency'].attrs = {
        'units': 'm/s',
        'long_name': 'Observed thickness tendency (dh/dt)',
        'comment': f'dh/dt computed from ITS_LIVE between '
                   f'{t0.strftime("%Y-%m-%d")} and '
                   f'{t1.strftime("%Y-%m-%d")}'}
    ds_out['observedThicknessTendencyUncertainty'].attrs = {
        'units': 'm/s',
        'long_name': 'Uncertainty in observed thickness tendency',
        'comment': 'Propagated from ITS_LIVE RMS errors at endpoints: '
                   'sqrt(rms_start^2 + rms_end^2) / dt'}
    ds_out['observedSurfaceElevation'].attrs = {
        'units': 'm',
        'long_name': 'Mean observed surface elevation',
        'comment': f'Time-mean surface elevation from ITS_LIVE between '
                   f'{t0.strftime("%Y-%m-%d")} and '
                   f'{t1.strftime("%Y-%m-%d")} '
                   f'({n_times} time steps)'}
    ds_out['observedSurfaceElevationUncertainty'].attrs = {
        'units': 'm',
        'long_name': 'Uncertainty in observed surface elevation',
        'comment': f'Mean RMS error from ITS_LIVE over '
                   f'{t0.strftime("%Y-%m-%d")} to '
                   f'{t1.strftime("%Y-%m-%d")} '
                   f'({n_times} time steps)'}
    fill = 1.0e36
    encoding = {v: {'_FillValue': fill} for v in ds_out.data_vars}
    ds_out.to_netcdf(dhdt_filename, encoding=encoding)
    ds.close()

    print(f'  Wrote dh/dt and mean h fields to {dhdt_filename}')

    # Generate or reuse remapping weights
    # SCRIP files are only needed if we must generate the weights
    if weights_file is not None:
        weights_filename = weights_file
        print(f'Using existing weights file: {weights_filename}')
        # Check that the weights file method is consistent with requested method.
        # ESMF stores e.g. 'Conservative remapping' or 'Bilinear remapping'.
        _esmf_method_keywords = {
            'conserve': 'conservative',
            'bilinear': 'bilinear',
        }
        with xarray.open_dataset(weights_filename) as ds_wgt:
            file_method = ds_wgt.attrs.get('map_method', None)
            if file_method is not None:
                keyword = _esmf_method_keywords.get(method.lower(), method.lower())
                if keyword not in file_method.lower():
                    raise ValueError(
                        f'Weights file was generated with method '
                        f'"{file_method}" but --method={method} was '
                        f'requested. Regenerate weights or change --method.')
    else:
        # Create or reuse destination SCRIP file
        if mali_scrip_file is not None:
            mali_scrip = mali_scrip_file
            print(f'Using existing MALI SCRIP file: {mali_scrip}')
        else:
            mesh_base = os.path.splitext(mali_file)[0]
            mali_scrip = f'{mesh_base}_scrip.nc'
            print('Creating SCRIP file for destination mesh')
            scrip_from_mpas(mali_file, mali_scrip)

        # Create or reuse source SCRIP file
        if source_scrip_file is not None:
            source_scrip = source_scrip_file
            print(f'Using existing ITS_LIVE SCRIP file: {source_scrip}')
        else:
            source_scrip = 'itslive_dhdt.scrip.nc'
            print('Creating SCRIP file for ITS_LIVE dataset')
            args = ['create_scrip_file_from_planar_rectangular_grid',
                    '-i', dhdt_filename,
                    '-s', source_scrip,
                    '-p', proj,
                    '-r', '2']
            run_command(args)

        # Mask source SCRIP: exclude cells with invalid (NaN/fill) data
        masked_source_scrip = 'itslive_dhdt.scrip_masked.nc'
        print('Masking source SCRIP for data validity')
        mask_source_scrip(source_scrip, masked_source_scrip, dhdt_filename)

        weights_filename = f'itslive_to_MPAS_{method}_weights.nc'
        print(f'Generating ITS_LIVE -> MPAS remapping weights '
              f'(method={method})')
        args = ['ESMF_RegridWeightGen',
                '--source', masked_source_scrip,
                '--destination', mali_scrip,
                '--weight', weights_filename,
                '--method', method,
                '--netcdf4',
                '--dst_regional',
                '--src_regional',
                '--ignore_unmapped']
        run_command(args)

    # Apply weights using ncremap
    remapped_filename = 'itslive_dhdt_remapped.nc'
    print('Applying weights with ncremap')
    args = ['ncremap',
            '-m', weights_filename,
            '-i', dhdt_filename,
            '-o', remapped_filename]
    run_command(args)

    # Rename ncremap's default 'ncol' dimension to match MALI's 'nCells'
    args = ['ncrename', '-d', 'ncol,nCells', remapped_filename]
    run_command(args)

    # Clean up unphysical values and merge remapped data into MALI file.
    # Read remapped variables directly (avoids ncks importing extra ncremap
    # variables like area, lat, lon, lat_vertices, lon_vertices).
    itslive_var_list = ['observedThicknessTendency',
                        'observedThicknessTendencyUncertainty',
                        'observedSurfaceElevation',
                        'observedSurfaceElevationUncertainty']
    print('Reading remapped ITS_LIVE fields and merging into MALI file')

    # Read remapped data (xarray decodes fill values to NaN)
    ds_remap = xarray.open_dataset(remapped_filename)

    # Open MALI file
    ds_mali = xarray.open_dataset(mali_file)

    # Assign remapped variables with proper (Time, nCells) shape
    for vname in itslive_var_list:
        data = ds_remap[vname].values
        if data.ndim == 1:
            data = data[np.newaxis, :]  # add Time dimension
        ds_mali[vname] = (('Time', 'nCells'), data)

    ds_remap.close()

    # Extract arrays for cleanup
    bed = ds_mali['bedTopography'].values
    thickness = ds_mali['thickness'].values
    dhdt = ds_mali['observedThicknessTendency'].values.copy()
    dhdt_err = ds_mali['observedThicknessTendencyUncertainty'].values.copy()
    h_obs = ds_mali['observedSurfaceElevation'].values.copy()
    h_err = ds_mali['observedSurfaceElevationUncertainty'].values.copy()

    # Unphysical mask: NaN or fill-value sentinel
    dhdt_bad = np.isnan(dhdt) | (np.abs(dhdt) > 1.0e10)
    h_bad = np.isnan(h_obs) | (np.abs(h_obs) > 1.0e10)
    no_ice = (thickness == 0.0)

    print(f'  {int(dhdt_bad.sum())} unphysical observedThicknessTendency '
          f'cells set to 0')
    print(f'  {int(h_bad.sum())} unphysical observedSurfaceElevation '
          f'cells set to bedTopography')
    print(f'  {int(no_ice.sum())} cells with thickness==0 assigned max '
          f'observedThicknessTendencyUncertainty')

    dhdt[dhdt_bad] = 0.0
    dhdt_err = np.fmin(dhdt_err, 1.0)
    dhdt_err[dhdt_bad] = 1.0
    dhdt_err[no_ice] = 1.0
    h_obs[h_bad] = bed[h_bad]

    ds_mali['observedThicknessTendency'].values[:] = dhdt
    ds_mali['observedThicknessTendencyUncertainty'].values[:] = dhdt_err
    ds_mali['observedSurfaceElevation'].values[:] = h_obs
    ds_mali['observedSurfaceElevationUncertainty'].values[:] = h_err

    # Append command to netCDF history attribute
    timestamp = datetime.datetime.now(
        datetime.timezone.utc).strftime('%Y-%m-%d %H:%M:%S UTC')
    cmd_str = ' '.join(sys.argv)
    new_history = f'{timestamp}: {cmd_str}'
    old_history = ds_mali.attrs.get('history', '')
    if old_history:
        ds_mali.attrs['history'] = f'{new_history}\n{old_history}'
    else:
        ds_mali.attrs['history'] = new_history

    write_netcdf(ds_mali, mali_file)
    ds_mali.close()

    print('Done.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i', '--itslive', required=True,
                        help='Path to ITS_LIVE ice elevation change '
                             'netCDF file')
    parser.add_argument('-m', '--mali', required=True,
                        help='Path to MALI mesh netCDF file (fields will '
                             'be appended)')
    parser.add_argument('-s', '--start-date', required=True,
                        help='Start date for dh/dt calculation '
                             '(format: YYYY-MM-DD)')
    parser.add_argument('-e', '--end-date', required=True,
                        help='End date for dh/dt calculation '
                             '(format: YYYY-MM-DD)')
    parser.add_argument('-p', '--proj', default='gis-gimp',
                        help='Projection of source dataset '
                             '(default: gis-gimp)')
    parser.add_argument('--method', default='conserve',
                        choices=['conserve', 'bilinear'],
                        help='Interpolation method: conserve or bilinear '
                             '(default: conserve)')
    parser.add_argument('--mali-scrip', default=None,
                        help='Path to existing SCRIP file for MALI mesh '
                             '(skip generation if provided)')
    parser.add_argument('--source-scrip', default=None,
                        help='Path to existing SCRIP file for ITS_LIVE grid '
                             '(skip generation if provided)')
    parser.add_argument('--weights', default=None,
                        help='Path to existing ESMF remapping weights file '
                             '(skip generation if provided)')

    args = parser.parse_args()

    interpolate_itslive_dhdt(
        itslive_file=args.itslive,
        mali_file=args.mali,
        start_date=args.start_date,
        end_date=args.end_date,
        proj=args.proj,
        method=args.method,
        mali_scrip_file=args.mali_scrip,
        source_scrip_file=args.source_scrip,
        weights_file=args.weights)
