"""
This script has functions that are needed to post-process and write output
variables that are derived from MALI restart files (currently litempavg, the
depth-average ice temperature, which is not in the 2d output streams).
"""

import os
from subprocess import check_call

import cftime
import numpy as np
import xarray as xr
from netCDF4 import date2num

from process_state_variables_ismip7 import write_netcdf_2d_state_vars
from time_utils import get_calendar, _cftime_calendar

EXPECTED_RESTART_VARIABLES = [
    'xtime', 'simulationStartTime', 'cellMask', 'temperature',
]


def check_restart_files(files):
    """Check that restart files exist and contain the needed variables.
    (validate_mali_files is not used because restart files have no
    daysSinceStart.)"""
    if len(files) == 0:
        raise ValueError("No restart files provided.")
    for f in files:
        if not os.path.exists(f):
            raise FileNotFoundError(f"Restart file not found: {f}")
        with xr.open_dataset(f, decode_cf=False) as ds:
            missing = [v for v in EXPECTED_RESTART_VARIABLES if v not in ds]
        if missing:
            raise ValueError(
                f"Restart file '{f}' is missing variables: {missing}")
    print(f"Validated {len(files)} restart file(s).")


def _decode_str(values):
    return values.tobytes().decode('utf-8').strip().strip('\x00')


def process_restart_vars(files, input_file_mesh, tmp_file):
    """
    Depth-average the 3D temperature of Jan. 1 restart files and write a
    file on the MALI mesh, with state-style time and mask variables, for
    remapping.

    Parameters
    ----------
    files : list of str
        Restart files, each with one time level.
    input_file_mesh : str
        MALI file with the time-invariant layerThicknessFractions.
    tmp_file : str
        Temporary output file name.
    """
    with xr.open_dataset(input_file_mesh, decode_cf=False) as ds_mesh:
        if 'layerThicknessFractions' not in ds_mesh:
            raise ValueError(
                f"layerThicknessFractions not found in {input_file_mesh}")
        fractions = ds_mesh['layerThicknessFractions'].values.astype(float)
    if not np.isclose(fractions.sum(), 1.0):
        raise ValueError(
            f"layerThicknessFractions sum to {fractions.sum()}, not 1.")

    with xr.open_dataset(files[0], decode_cf=False) as ds_first:
        calendar = get_calendar(ds_first.attrs)
        start_str = ds_first['simulationStartTime'].values
        start_time = _decode_str(start_str[0] if start_str.ndim > 1
                                 else start_str)
        start_date = start_time.split('_')[0]
        start_chars = ds_first['simulationStartTime'].values.copy()
    if start_chars.ndim > 1:
        start_chars = start_chars[0]
    cal = _cftime_calendar(calendar)

    days = []
    litempavg = []
    dynamic = []
    floating = []
    for f in files:
        with xr.open_dataset(f, decode_cf=False) as ds:
            xtime = _decode_str(ds['xtime'].values[0])
            date_part, time_part = xtime.split('_')
            year, month, day = [int(v) for v in date_part.split('-')]
            if (month, day) != (1, 1):
                print(f"Skipping restart file not at Jan. 1: {f}")
                continue
            d = date2num(cftime.datetime(year, month, day, calendar=cal),
                         units=f'days since {start_date}', calendar=cal)
            if np.isclose(d, 0.0):
                print(f"Skipping restart data at simulation start: {f}")
                continue
            temp = ds['temperature'].values[0]  # nCells, nVertLevels
            mask = ds['cellMask'].values[0]
            dyn = ((mask & 2) / 2).astype(float)
            flt = ((mask & 4) / 4).astype(float)
            avg = np.tensordot(temp, fractions, axes=([1], [0]))
            avg = np.where(dyn > 0, avg, 0.0)
            days.append(float(d))
            litempavg.append(avg)
            dynamic.append(dyn)
            floating.append(flt)

    if not days:
        raise ValueError("No Jan. 1 restart records available.")
    order = np.argsort(days)
    days = np.asarray(days)[order]
    litempavg = np.stack(litempavg)[order]
    dynamic = np.stack(dynamic)[order]
    floating = np.stack(floating)[order]

    ds_out = xr.Dataset(
        {
            'daysSinceStart': (('Time',), days),
            'simulationStartTime': (('StrLen',), start_chars),
            'litempavg': (('Time', 'nCells'), litempavg * dynamic),
            'sftgif': (('Time', 'nCells'), dynamic),
            'sftflf': (('Time', 'nCells'), floating * dynamic),
            'sftgrf': (('Time', 'nCells'), (1.0 - floating) * dynamic),
        })
    ds_out.to_netcdf(tmp_file)
    ds_out.close()
    return calendar


def process_restart_pipeline(restart_files, input_file_mesh, mapping_file,
                             ismip7_grid_file, output_path, metadata):
    """
    Full restart-derived variable pipeline: validate, depth-average, remap,
    write.  Writes litempavg.
    """
    metadata = metadata.copy()
    check_restart_files(restart_files)

    tmp_file = 'tmp_restart.nc'
    remapped_file = 'processed_and_remapped_restart.nc'
    calendar = process_restart_vars(restart_files, input_file_mesh, tmp_file)

    print("Remapping restart-derived file.")
    check_call(["ncremap", "-i", tmp_file, "-o", remapped_file,
                "-m", mapping_file, "-P", "mpas"])

    print("* Writing litempavg")
    write_netcdf_2d_state_vars('litempavg', 'litempavg', '', 'K',
                               'Depth average temperature',
                               remapped_file, ismip7_grid_file,
                               output_path, metadata, calendar=calendar)

    os.remove(tmp_file)
    os.remove(remapped_file)
