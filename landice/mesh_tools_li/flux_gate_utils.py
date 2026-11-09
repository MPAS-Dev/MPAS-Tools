#!/usr/bin/env python3
"""
Shared helpers for matching MALI flux-gate masks (see make_flux_gate_masks.py) to GEUS
(Mankoff et al., 2020) observed discharge, and for reading MALI fluxGatesOutput files.

Used by output_processing_li/plot_flux_gates.py and
mesh_tools_li/adjust_mu_friction_from_flux_gates.py.

Created: 2026-10-01
@author: Trevor Hillebrand
"""

import os
import sqlite3
import sys

import numpy as np
import pandas as pd
from netCDF4 import Dataset, chartostring

DEFAULT_OBS_DIR = '/global/cfs/cdirs/m4288/users/trhille/ISMIP7/gris_flux_gates'
DAYS_IN_MONTH = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]


def xtime_to_decimal_year(xtime_strs):
    """Convert MPAS xtime strings (YYYY-MM-DD_hh:mm:ss) to decimal years, assuming a 365-day calendar."""
    years = np.zeros(len(xtime_strs))
    for i, s in enumerate(xtime_strs):
        date_part, time_part = s.strip().split('_')
        y, m, d = (int(v) for v in date_part.split('-'))
        hh, mm, ss = (int(v) for v in time_part.split(':'))
        day_of_year = sum(DAYS_IN_MONTH[:m - 1]) + (d - 1) + hh / 24.0 + mm / 1440.0 + ss / 86400.0
        years[i] = y + day_of_year / 365.0
    return years


def load_mali_flux(flux_files):
    """
    Read and concatenate iceFluxThroughGates from one or more MALI fluxGatesOutput files.

    Returns
    -------
    years : (nTime,) ndarray of decimal years, sorted
    flux : (nTime, nFluxGates) ndarray, ice flux through each gate in Gt/yr
    """
    all_years = []
    all_flux = []
    for fname in flux_files:
        with Dataset(fname, 'r') as f:
            years = xtime_to_decimal_year(chartostring(f.variables['xtime'][:]))
            flux = f.variables['iceFluxThroughGates'][:, :] / 1.0e12  # kg/yr -> Gt/yr
            all_years.append(years)
            all_flux.append(flux)
    years = np.concatenate(all_years)
    flux = np.concatenate(all_flux, axis=0)
    order = np.argsort(years)
    return years[order], flux[order, :]


def load_mask_gate_ids(mask_file):
    """Read fluxGateIds from the mask file produced by make_flux_gate_masks.py."""
    with Dataset(mask_file, 'r') as f:
        if 'fluxGateIds' not in f.variables:
            sys.exit(f"ERROR: {mask_file} has no fluxGateIds variable. "
                     "Regenerate it with the current make_flux_gate_masks.py.")
        gate_ids = [int(g) for g in f.variables['fluxGateIds'][:]]
    return gate_ids


def load_mask_gate_locations(gpkg_path, gate_ids, table='gates_final', gate_field='gate',
                             x_field='mean_x', y_field='mean_y', name_field='Mouginot_2019'):
    """Read one representative (mean_x, mean_y, name) per mask gate_id from the gpkg."""
    conn = sqlite3.connect(gpkg_path)
    cursor = conn.cursor()
    locations = {}
    for gid in gate_ids:
        cursor.execute(
            f"SELECT {x_field}, {y_field}, {name_field} FROM {table} "
            f"WHERE {gate_field}=? LIMIT 1",
            (gid,)
        )
        row = cursor.fetchone()
        if row is not None:
            locations[gid] = (row[0], row[1], row[2])
    conn.close()
    return locations


def match_gates_by_location(mask_locations, meta, tolerance=1000.0):
    """
    Match each mask gate to the nearest gate_meta.csv gate by (mean_x, mean_y), since gate
    ID numbers can differ between the gpkg used to build the mask and gate_meta.csv /
    gate_D.csv (e.g. the same physical gate can be numbered 27 in one and 24 in the other).

    Returns
    -------
    dict {mask_gate_id: obs_gate_id}, only gates matched within tolerance (meters)
    """
    obs_ids = meta.index.values
    obs_xy = meta[['mean_x', 'mean_y']].values

    match = {}
    unmatched = []
    renamed = []
    for mask_gid, (x, y, name) in mask_locations.items():
        dist = np.hypot(obs_xy[:, 0] - x, obs_xy[:, 1] - y)
        i = np.argmin(dist)
        if dist[i] > tolerance:
            unmatched.append(mask_gid)
            continue
        obs_gid = int(obs_ids[i])
        match[mask_gid] = obs_gid
        obs_name = meta.loc[obs_gid, 'Mouginot_2019']
        if name and obs_name and name != obs_name:
            renamed.append((mask_gid, name, obs_gid, obs_name))

    if unmatched:
        print(f"WARNING: {len(unmatched)} mask gate(s) had no gate_meta.csv match within "
              f"{tolerance:.0f} m, excluded from observations: {unmatched}")
    if renamed:
        print(f"WARNING: {len(renamed)} gate(s) matched by location but have different names:")
        for mask_gid, name, obs_gid, obs_name in renamed:
            print(f"  mask gate {mask_gid} ('{name}') -> obs gate {obs_gid} ('{obs_name}')")

    return match


def load_observations(obs_dir, gpkg, gate_ids, tolerance=1000.0):
    """
    Load GEUS gate discharge, errors, and region metadata, matched to mask gate_ids by
    location rather than by gate ID (see match_gates_by_location).

    Returns
    -------
    region_of_gate : dict {mask_gate_id: region}
    gate_D, gate_err : DataFrames indexed by Date, columns are mask gate_ids (int)
    """
    meta = pd.read_csv(os.path.join(obs_dir, 'gate_meta.csv'), index_col='gate')
    gate_D_raw = pd.read_csv(os.path.join(obs_dir, 'gate_D.csv'), index_col='Date', parse_dates=True)
    gate_err_raw = pd.read_csv(os.path.join(obs_dir, 'gate_err.csv'), index_col='Date', parse_dates=True)
    gate_D_raw.columns = gate_D_raw.columns.astype(int)
    gate_err_raw.columns = gate_err_raw.columns.astype(int)

    mask_locations = load_mask_gate_locations(gpkg, gate_ids)
    match = match_gates_by_location(mask_locations, meta, tolerance=tolerance)

    region_of_gate = {}
    no_column = []
    for mask_gid, obs_gid in match.items():
        if obs_gid not in gate_D_raw.columns:
            no_column.append((mask_gid, obs_gid))
            continue
        region_of_gate[mask_gid] = meta.loc[obs_gid, 'region']
    if no_column:
        print(f"WARNING: {len(no_column)} matched obs gate(s) not in gate_D.csv, excluded "
              f"(mask_id, obs_id): {no_column}")

    valid_mask_gids = list(region_of_gate.keys())
    gate_D = gate_D_raw[[match[g] for g in valid_mask_gids]].copy()
    gate_D.columns = valid_mask_gids
    gate_err = gate_err_raw[[match[g] for g in valid_mask_gids]].copy()
    gate_err.columns = valid_mask_gids

    return region_of_gate, gate_D, gate_err
