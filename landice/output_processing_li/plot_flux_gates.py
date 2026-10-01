#!/usr/bin/env python3
"""
Plot MALI flux-gate ice discharge against GEUS (Mankoff et al., 2020) observations,
aggregated by Mouginot & Rignot (2019) drainage region.

Reads one or more MALI fluxGatesOutput files (iceFluxThroughGates), looks up each
gate's region from the GEUS gate_meta.csv using fluxGateIds in the mask file (see
mesh_tools_li/make_flux_gate_masks.py), and compares regional totals against GEUS
gate_D.csv / gate_err.csv, restricted to only the gates present in the mask.

Example usage:
    ./plot_flux_gates.py \\
        -f fluxGates.nc \\
        -m fluxGateMasks.nc \\
        -o flux_gates.png

Created: 2026-10-01
@author: Trevor Hillebrand
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from netCDF4 import Dataset, chartostring

REGIONS = ['NO', 'NE', 'CE', 'SE', 'SW', 'CW', 'NW']
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
        gate_ids = f.variables['fluxGateIds'][:].astype(int)
    return gate_ids


def load_observations(obs_dir, gate_ids):
    """
    Load GEUS gate discharge, errors, and region metadata, restricted to gate_ids.

    Returns
    -------
    region_of_gate : dict {gate_id: region}, only gates found in both gate_meta.csv and gate_D.csv
    gate_D, gate_err : DataFrames indexed by Date, columns are gate IDs (int)
    """
    meta = pd.read_csv(os.path.join(obs_dir, 'gate_meta.csv'), index_col='gate')
    gate_D = pd.read_csv(os.path.join(obs_dir, 'gate_D.csv'), index_col='Date', parse_dates=True)
    gate_err = pd.read_csv(os.path.join(obs_dir, 'gate_err.csv'), index_col='Date', parse_dates=True)
    gate_D.columns = gate_D.columns.astype(int)
    gate_err.columns = gate_err.columns.astype(int)

    region_of_gate = {}
    no_meta = []
    no_column = []
    for gid in gate_ids:
        if gid not in meta.index:
            no_meta.append(gid)
            continue
        if gid not in gate_D.columns:
            no_column.append(gid)
            continue
        region_of_gate[gid] = meta.loc[gid, 'region']
    if no_meta:
        print(f"WARNING: {len(no_meta)} mask gate(s) not in gate_meta.csv, excluded from "
              f"observations: {no_meta}")
    if no_column:
        print(f"WARNING: {len(no_column)} mask gate(s) not in gate_D.csv, excluded from "
              f"observations: {no_column}")

    return region_of_gate, gate_D, gate_err


def aggregate_mali(flux, gate_ids, region_of_gate):
    """Group gate column indices by region and sum MALI flux within each region."""
    region_gate_indices = {r: [] for r in REGIONS}
    for idx, gid in enumerate(gate_ids):
        region = region_of_gate.get(gid)
        if region in region_gate_indices:
            region_gate_indices[region].append(idx)

    region_totals = {
        r: flux[:, idx].sum(axis=1) if idx else None
        for r, idx in region_gate_indices.items()
    }
    gis_total = flux.sum(axis=1)
    return region_gate_indices, region_totals, gis_total


def aggregate_obs(gate_D, gate_err, region_of_gate):
    """Sum observed discharge by region, combining errors by linear sum (Mankoff et al., 2020)."""
    region_gates = {r: [gid for gid, reg in region_of_gate.items() if reg == r] for r in REGIONS}

    region_D = {}
    region_err = {}
    for r, gids in region_gates.items():
        if gids:
            region_D[r] = gate_D[gids].sum(axis=1)
            region_err[r] = gate_err[gids].sum(axis=1)
        else:
            region_D[r] = None
            region_err[r] = None

    gis_gids = list(region_of_gate.keys())
    gis_D = gate_D[gis_gids].sum(axis=1)
    gis_err = gate_err[gis_gids].sum(axis=1)
    return region_gates, region_D, region_err, gis_D, gis_err


def plot_panel(ax, mali_years, flux, gate_indices, mali_total,
               obs_years, obs_D, obs_err, title):
    """Plot per-gate MALI fluxes (gray), the MALI regional total (black), and GEUS obs (blue)."""
    for idx in gate_indices:
        ax.plot(mali_years, flux[:, idx], color='gray', linewidth=0.5, alpha=0.6)
    if mali_total is not None:
        ax.plot(mali_years, mali_total, color='black', linewidth=1.5, label='MALI total')
    if obs_D is not None:
        ax.plot(obs_years, obs_D.values, color='tab:blue', linewidth=2, label='GEUS obs.')
        if obs_err is not None:
            ax.fill_between(obs_years, obs_D.values - obs_err.values, obs_D.values + obs_err.values,
                             color='tab:blue', alpha=0.3, linewidth=0)
    ax.set_title(title)
    ax.set_ylabel('Discharge (Gt yr$^{-1}$)')


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('-f', '--flux-files', nargs='+', required=True,
                        help='MALI fluxGatesOutput file(s) (e.g. fluxGates.nc)')
    parser.add_argument('-m', '--mask-file', required=True,
                        help='Flux-gate mask file with fluxGateIds (from make_flux_gate_masks.py)')
    parser.add_argument('--obs-dir', default=DEFAULT_OBS_DIR,
                        help=f'Directory with GEUS gate_D.csv, gate_err.csv, gate_meta.csv '
                             f'(default: {DEFAULT_OBS_DIR})')
    parser.add_argument('-o', '--output', default='flux_gates.png',
                        help='Output figure filename (default: flux_gates.png)')
    parser.add_argument('--show', action='store_true',
                        help='Display the figure interactively')
    args = parser.parse_args(argv)

    print(f"Reading MALI flux output: {args.flux_files}")
    years, flux = load_mali_flux(args.flux_files)

    print(f"Reading gate IDs from mask file: {args.mask_file}")
    gate_ids = load_mask_gate_ids(args.mask_file)
    if flux.shape[1] != len(gate_ids):
        sys.exit(f"ERROR: nFluxGates mismatch between flux file(s) ({flux.shape[1]}) "
                 f"and mask file ({len(gate_ids)})")

    print(f"Reading GEUS observations from: {args.obs_dir}")
    region_of_gate, gate_D, gate_err = load_observations(args.obs_dir, gate_ids)
    obs_years = gate_D.index.year + (gate_D.index.dayofyear - 1) / 365.0

    region_gate_indices, mali_region_totals, mali_gis_total = aggregate_mali(
        flux, gate_ids, region_of_gate)
    _, obs_region_D, obs_region_err, obs_gis_D, obs_gis_err = aggregate_obs(
        gate_D, gate_err, region_of_gate)

    fig, axs = plt.subplots(2, 4, figsize=(18, 8), sharex=True)
    axs = axs.ravel()

    for ax, region in zip(axs[:7], REGIONS):
        indices = region_gate_indices[region]
        if not indices:
            ax.set_title(f'{region} (no gates in mask)')
            ax.axis('off')
            continue
        gate_word = 'gate' if len(indices) == 1 else 'gates'
        plot_panel(ax, years, flux, indices, mali_region_totals[region],
                   obs_years, obs_region_D[region], obs_region_err[region],
                   f'{region} ({len(indices)} {gate_word})')

    plot_panel(axs[7], years, flux, list(range(len(gate_ids))), mali_gis_total,
               obs_years, obs_gis_D, obs_gis_err, 'Greenland total')

    handles, labels = axs[7].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', ncol=2)
    fig.suptitle('MALI vs. GEUS (Mankoff et al., 2020) ice discharge by region')
    fig.tight_layout(rect=(0, 0.05, 1, 1))

    fig.savefig(args.output, dpi=150)
    print(f"Wrote {args.output}")
    if args.show:
        plt.show()


if __name__ == '__main__':
    main()
