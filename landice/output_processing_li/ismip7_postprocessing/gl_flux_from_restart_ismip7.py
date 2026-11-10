"""
Fallback approximation of the MALI grounding-line flux (avgGroundingLineFlux,
used for ISMIP ligroundf) from restart files.

This is for runs made with a MALI version whose avgGroundingLineFlux was
divided by the grounded-cell thickness (huge values at very thin ice).  It
reproduces the corrected MALI calculation from a snapshot of the state:

    layerNormalVelocity(k) = 0.5 * (normalVelocity(k) + normalVelocity(k+1))
    layerThicknessEdge(k)  = first-order upwind of thickness*fraction(k)
    F (m2 s-1)             = sign * sum_k layerNormalVelocity*layerThicknessEdge
    flux on grounded cell  = sum_edges F * dvEdge * rho_ice / areaCell

with sign > 0 for flow from the grounded to the floating cell.  The
mean of the snapshots from all restarts within each averaging period (e.g.
a year) approximates the time-averaged flux.
"""

import cftime
import numpy as np
import xarray as xr
from netCDF4 import date2num

from time_utils import get_calendar, _cftime_calendar

RHO_ICE = 910.0  # MALI default config_ice_density, kg m-3

# cellMask bits (mpas_li_mask.F)
_ICE = 32
_DYNAMIC = 2
_FLOATING = 4


def _decode_str(values):
    return values.tobytes().decode('utf-8').strip().strip('\x00')


def _edge_normal_velocity(ds):
    """Normal velocity on edges at layer interfaces, (nEdges, nVertInterfaces).
    Uses the restart's normalVelocity if present; otherwise averages the
    reconstructed cell-center velocities to edges and projects them on the
    edge normal."""
    if 'normalVelocity' in ds:
        return ds['normalVelocity'].values[0].astype(float)
    for v in ['uReconstructX', 'uReconstructY', 'angleEdge', 'cellsOnEdge']:
        if v not in ds:
            raise ValueError(
                f"Restart lacks normalVelocity and {v}; cannot compute the "
                "grounding-line flux.")
    ux = ds['uReconstructX'].values[0].astype(float)
    uy = ds['uReconstructY'].values[0].astype(float)
    coe = ds['cellsOnEdge'].values - 1
    angle = ds['angleEdge'].values.astype(float)[:, None]
    c1 = np.where(coe[:, 0] >= 0, coe[:, 0], coe[:, 1])
    c2 = np.where(coe[:, 1] >= 0, coe[:, 1], coe[:, 0])
    ux_e = 0.5 * (ux[c1] + ux[c2])
    uy_e = 0.5 * (uy[c1] + uy[c2])
    return ux_e * np.cos(angle) + uy_e * np.sin(angle)


def compute_gl_flux(ds, fractions, rho_ice=RHO_ICE):
    """
    Grounding-line flux on cells (kg m-2 s-1, positive from grounded to
    floating ice) from one restart record.

    Parameters
    ----------
    ds : xarray.Dataset
        Restart file (decode_cf=False) with thickness, cellMask, cellsOnEdge,
        dvEdge, areaCell and normalVelocity (or uReconstructX/Y, angleEdge).
    fractions : numpy.ndarray
        layerThicknessFractions, (nVertLevels,).
    """
    thickness = ds['thickness'].values[0].astype(float)
    mask = ds['cellMask'].values[0].astype(int)
    coe = ds['cellsOnEdge'].values - 1  # 0 (no cell) becomes -1
    dv = ds['dvEdge'].values.astype(float)
    area = ds['areaCell'].values.astype(float)

    nv = _edge_normal_velocity(ds)
    if nv.shape[1] != len(fractions) + 1:
        raise ValueError(
            f"normalVelocity has {nv.shape[1]} interfaces but "
            f"layerThicknessFractions has {len(fractions)} layers.")
    lnv = 0.5 * (nv[:, :-1] + nv[:, 1:])  # edges, layers

    valid = (coe[:, 0] >= 0) & (coe[:, 1] >= 0)
    c1 = np.where(valid, coe[:, 0], 0)
    c2 = np.where(valid, coe[:, 1], 0)

    ice = (mask & _ICE) > 0
    floating = (mask & _FLOATING) > 0
    dynamic = (mask & _DYNAMIC) > 0
    grounded = ice & ~floating
    float_dyn = floating & dynamic

    g1 = grounded[c1]
    g2 = grounded[c2]
    gl = valid & ((g1 & float_dyn[c2]) | (g2 & float_dyn[c1]))
    sign = np.where(g1, 1.0, -1.0)  # + is cell 1 -> cell 2
    grounded_cell = np.where(g1, c1, c2)

    # first-order upwind layer thickness on edges
    h1 = thickness[c1][:, None] * fractions[None, :]
    h2 = thickness[c2][:, None] * fractions[None, :]
    vsign = np.sign(lnv)
    h_edge = np.maximum(vsign * h1, -vsign * h2)

    flux_edge = sign * np.sum(lnv * h_edge, axis=1)  # m2 s-1
    flux_edge = np.where(gl, flux_edge, 0.0)

    out = np.zeros(thickness.shape)
    np.add.at(out, grounded_cell[gl],
              flux_edge[gl] * dv[gl] * rho_ice / area[grounded_cell[gl]])
    return out


def gl_flux_from_restarts(restart_files, input_file_mesh, flux_days,
                          rho_ice=RHO_ICE):
    """
    Build avgGroundingLineFlux (nTime, nCells; MALI sign) from restarts.

    The flux for each flux record (time-average over the preceding period,
    ending at flux_days[i]) is approximated by the mean of the snapshot
    fluxes of all restart files dated in (flux_days[i-1], flux_days[i]].
    With one restart per year this is a single snapshot; more frequent
    restarts give a better approximation of the average.

    Parameters
    ----------
    restart_files : list of str
    input_file_mesh : str
        File with layerThicknessFractions.
    flux_days : array-like
        daysSinceStart of the flux records to be replaced.
    """
    with xr.open_dataset(input_file_mesh, decode_cf=False) as ds_mesh:
        if 'layerThicknessFractions' not in ds_mesh:
            raise ValueError(
                f"layerThicknessFractions not found in {input_file_mesh}")
        fractions = ds_mesh['layerThicknessFractions'].values.astype(float)

    flux_days = np.asarray(flux_days, dtype=float)
    restart_days = {}
    start_date = None
    cal = None
    for f in restart_files:
        with xr.open_dataset(f, decode_cf=False) as ds:
            if start_date is None:
                calendar = get_calendar(ds.attrs)
                cal = _cftime_calendar(calendar)
                s = ds['simulationStartTime'].values
                start_date = _decode_str(
                    s[0] if s.ndim > 1 else s).split('_')[0]
            xtime = _decode_str(ds['xtime'].values[0])
            date_part, time_part = xtime.split('_')
            year, month, day = [int(v) for v in date_part.split('-')]
            hh, mm, ss = [int(float(v)) for v in time_part.split(':')]
            t = cftime.datetime(year, month, day, hh, mm, ss, calendar=cal)
            restart_days[f] = float(date2num(
                t, units=f'days since {start_date}', calendar=cal))

    out = None
    missing = []
    for i, d in enumerate(flux_days):
        if i == 0 and np.isclose(d, 0.0):
            continue  # start of simulation: nothing averaged yet
        d_prev = flux_days[i - 1] if i > 0 else -np.inf
        files = [f for f, rd in restart_days.items()
                 if d_prev + 0.5 < rd <= d + 0.5]
        if not files:
            missing.append(d)
            continue
        acc = None
        for f in files:
            with xr.open_dataset(f, decode_cf=False) as ds:
                flux = compute_gl_flux(ds, fractions, rho_ice)
            acc = flux if acc is None else acc + flux
        if out is None:
            out = np.zeros((len(flux_days), acc.size))
        out[i] = acc / len(files)
        print(f"  ligroundf approximated from {len(files)} restart(s) "
              f"for flux time {d:.0f}")
    if missing:
        raise ValueError(
            "No restart file falls in the averaging period ending at flux "
            f"times (days since start): {missing}. Check --restart_dir and "
            "--ligroundf_restart_pattern.")
    if out is None:
        raise ValueError("No restart records available.")
    return out
