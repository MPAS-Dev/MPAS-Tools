"""
Fallback correction of the MALI avgGroundingLineFlux (used for ISMIP
ligroundf) for runs made with a MALI version that divided the flux by the
grounded-cell thickness.

Old MALI (per time step, then time-averaged):
    fluxAcrossGroundingLineOnCells(cell) = sum_edges F * rho_ice / H(cell)
Current MALI:
    fluxAcrossGroundingLineOnCells(cell) = sum_edges F * dvEdge * rho_ice
                                           / areaCell(cell)
where F (m2 s-1) is the grounding-line flux per unit width.  The corrected
field is therefore approximately

    new = old * H * dvEdge / areaCell

with dvEdge taken as the mean over the cell's edges (exact for cells with
one grounding-line edge) and H a representative thickness of the grounded
cell over the averaging period: the mean of the thickness at the end of
the period (from the flux file) and from any restart files dated within the
period.  The time-averaged old field was divided by the instantaneous
thickness at each step, so the result is approximate (a few percent against
the globalStats grounding-line flux in tests).

Cells where the old field was dominated by a division by a nearly zero
thickness cannot be recovered; where |new| exceeds MAX_ABS_FLUX (the
ISMIP7 checker limit) the cell is set to 0.

A previous approach, rebuilding the flux from restart velocities, was
removed because it was about 45% below globalStats; see git history
(commit da22e7d7, gl_flux_from_restart_ismip7.py).
"""

import cftime
import numpy as np
import xarray as xr
from netCDF4 import date2num

from time_utils import get_calendar, _cftime_calendar

MAX_ABS_FLUX = 10.0  # kg m-2 s-1, ISMIP7 range limit for ligroundf


def _decode_str(values):
    return values.tobytes().decode('utf-8').strip().strip('\x00')


def _restart_days(restart_files):
    """daysSinceStart of each restart file (dict file -> days)."""
    days = {}
    start_date = None
    cal = None
    for f in restart_files:
        with xr.open_dataset(f, decode_cf=False) as ds:
            if start_date is None:
                cal = _cftime_calendar(get_calendar(ds.attrs))
                s = ds['simulationStartTime'].values
                start_date = _decode_str(
                    s[0] if s.ndim > 1 else s).split('_')[0]
            xtime = _decode_str(ds['xtime'].values[0])
            date_part, time_part = xtime.split('_')
            year, month, day = [int(v) for v in date_part.split('-')]
            hh, mm, ss = [int(float(v)) for v in time_part.split(':')]
            t = cftime.datetime(year, month, day, hh, mm, ss, calendar=cal)
            days[f] = float(date2num(
                t, units=f'days since {start_date}', calendar=cal))
    return days


def correct_gl_flux(old_flux, thickness, flux_days, input_file_mesh,
                    restart_files=None):
    """
    Correct the thickness division in avgGroundingLineFlux.

    Parameters
    ----------
    old_flux : numpy.ndarray
        avgGroundingLineFlux, (nTime, nCells), kg m-2 s-1, MALI sign.
    thickness : numpy.ndarray
        thickness at the flux records, (nTime, nCells).
    flux_days : array-like
        daysSinceStart of the flux records.
    input_file_mesh : str
        File with areaCell, edgesOnCell, nEdgesOnCell and dvEdge.
    restart_files : list of str, optional
        Restart files used to average the thickness over each period.

    Returns
    -------
    numpy.ndarray, (nTime, nCells), corrected flux (MALI sign).
    """
    with xr.open_dataset(input_file_mesh, decode_cf=False) as ds_mesh:
        for v in ['areaCell', 'edgesOnCell', 'nEdgesOnCell', 'dvEdge']:
            if v not in ds_mesh:
                raise ValueError(f"{v} not found in {input_file_mesh}")
        area = ds_mesh['areaCell'].values.astype(float)
        edges = ds_mesh['edgesOnCell'].values - 1
        n_edges = ds_mesh['nEdgesOnCell'].values
        dv_edge = ds_mesh['dvEdge'].values.astype(float)

    valid = np.arange(edges.shape[1])[None, :] < n_edges[:, None]
    dv_cell = (np.where(valid, dv_edge[np.where(valid, edges, 0)], 0.0)
               .sum(axis=1) / n_edges)
    weight = dv_cell / area

    flux_days = np.asarray(flux_days, dtype=float)
    rdays = _restart_days(restart_files) if restart_files else {}

    out = np.zeros(old_flux.shape)
    for i, d in enumerate(flux_days):
        if i == 0 and np.isclose(d, 0.0):
            continue  # start of simulation: nothing averaged yet
        d_prev = flux_days[i - 1] if i > 0 else -np.inf
        h = [thickness[i].astype(float)]
        for f, rd in rdays.items():
            if d_prev + 0.5 < rd <= d + 0.5:
                with xr.open_dataset(f, decode_cf=False) as ds:
                    h.append(ds['thickness'].values[0].astype(float))
        new = old_flux[i] * np.mean(h, axis=0) * weight
        bad = ~np.isfinite(new) | (np.abs(new) > MAX_ABS_FLUX)
        new[bad] = 0.0
        out[i] = new
        print(f"  ligroundf corrected for flux time {d:.0f} using "
              f"{len(h)} thickness field(s); {bad.sum()} cell(s) beyond "
              f"+/-{MAX_ABS_FLUX:g} set to 0")
    return out
