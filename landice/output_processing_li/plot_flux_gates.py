#!/usr/bin/env python3
"""
Plot MALI flux-gate ice discharge against GEUS (Mankoff et al., 2020) observations,
aggregated by Mouginot & Rignot (2019) drainage region.

Reads one or more MALI fluxGatesOutput files per run (optionally several runs, for
comparison) and the gpkg used to build the mask file (see
mesh_tools_li/make_flux_gate_masks.py), then matches each mask gate to the nearest gate
in GEUS gate_meta.csv by (mean_x, mean_y) location rather than by gate ID number, since
the two datasets can assign different numbers to the same physical gate. Compares
regional totals against GEUS gate_D.csv / gate_err.csv.

Example usage (single run):
    ./plot_flux_gates.py \\
        -f fluxGates.nc \\
        -m fluxGateMasks.nc \\
        -g gates.gpkg \\
        -o flux_gates.png

Example usage (comparing two runs, each split across yearly files):
    ./plot_flux_gates.py \\
        -f "run1/output/fluxGates_*.nc" "run2/output/fluxGates_*.nc" \\
        --labels run1 run2 \\
        -m fluxGateMasks.nc \\
        -g gates.gpkg \\
        -o flux_gates.png

Created: 2026-10-01
@author: Trevor Hillebrand
"""

import argparse
import glob
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'mesh_tools_li'))
from flux_gate_utils import (  # noqa: E402
    DEFAULT_OBS_DIR,
    load_mali_flux,
    load_mask_gate_ids,
    load_observations,
)

REGIONS = ['NO', 'NE', 'CE', 'SE', 'SW', 'CW', 'NW']
# Colors for run totals/per-gate lines when comparing multiple runs; tab:blue is
# reserved for GEUS observations.
RUN_COLORS = ['tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown',
              'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']


def expand_run(pattern):
    """Expand a single run specifier (a file path or a glob pattern) into its file list."""
    matches = sorted(glob.glob(pattern))
    return matches if matches else [pattern]


def get_run_styles(n_runs):
    """Return (total_color, gate_color) per run; black/gray for a single run, else RUN_COLORS."""
    if n_runs == 1:
        return [('black', 'gray')]
    return [(RUN_COLORS[i % len(RUN_COLORS)], RUN_COLORS[i % len(RUN_COLORS)]) for i in range(n_runs)]


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


def plot_panel(ax, mali_runs, gate_indices, obs_years, obs_gate_D, obs_D, obs_err, title):
    """
    Plot per-gate MALI fluxes, each run's regional/GIS total, and GEUS obs (blue).

    mali_runs : list of dict, each with keys 'label', 'years', 'flux', 'total',
        'color', 'gate_color'. 'total' is the region/GIS total for that run, or None.
    """
    for run in mali_runs:
        for idx in gate_indices:
            ax.plot(run['years'], run['flux'][:, idx], color=run['gate_color'],
                     linewidth=0.5, alpha=0.5)
        if run['total'] is not None:
            ax.plot(run['years'], run['total'], color=run['color'], linewidth=1.5,
                     label=run['label'])
    if obs_gate_D is not None:
        for gid in obs_gate_D.columns:
            ax.plot(obs_years, obs_gate_D[gid].values, color='tab:blue', linewidth=0.5, alpha=0.4)
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
                        help='One run specifier per run, space-separated: a single '
                             'fluxGates.nc file, or a quoted glob pattern matching the '
                             'multiple files of one run (e.g. "run1/output/fluxGates_*.nc"). '
                             'Give multiple specifiers to compare runs.')
    parser.add_argument('--labels', nargs='+',
                        help='Label for each run, in the order given by -f (default: '
                             '"MALI total" for a single run, or "Run 1", "Run 2", ... '
                             'for multiple runs)')
    parser.add_argument('-m', '--mask-file', required=True,
                        help='Flux-gate mask file with fluxGateIds (from make_flux_gate_masks.py)')
    parser.add_argument('-g', '--gpkg', required=True,
                        help='GeoPackage used to build the mask file, for matching mask gates '
                             'to GEUS gate_meta.csv by location (gate ID numbers can differ '
                             'between the two)')
    parser.add_argument('--match-tolerance', type=float, default=1000.0,
                        help='Max distance (m) between mask and gate_meta.csv gate centroids '
                             'to consider them the same physical gate (default: 1000)')
    parser.add_argument('--obs-dir', default=DEFAULT_OBS_DIR,
                        help=f'Directory with GEUS gate_D.csv, gate_err.csv, gate_meta.csv '
                             f'(default: {DEFAULT_OBS_DIR})')
    parser.add_argument('-o', '--output', default='flux_gates.png',
                        help='Output figure filename (default: flux_gates.png)')
    parser.add_argument('--show', action='store_true',
                        help='Display the figure interactively')
    args = parser.parse_args(argv)

    run_file_lists = [expand_run(pattern) for pattern in args.flux_files]
    if args.labels:
        if len(args.labels) != len(run_file_lists):
            parser.error(f"--labels count ({len(args.labels)}) must match number of runs "
                         f"({len(run_file_lists)}, from -f arguments)")
        labels = args.labels
    elif len(run_file_lists) == 1:
        labels = ['MALI total']
    else:
        labels = [f'Run {i + 1}' for i in range(len(run_file_lists))]

    print(f"Reading gate IDs from mask file: {args.mask_file}")
    gate_ids = load_mask_gate_ids(args.mask_file)

    mali_runs_data = []
    for label, run_files in zip(labels, run_file_lists):
        print(f"Reading MALI flux output for '{label}': {run_files}")
        run_years, run_flux = load_mali_flux(run_files)
        if run_flux.shape[1] != len(gate_ids):
            sys.exit(f"ERROR: nFluxGates mismatch for run '{label}' ({run_flux.shape[1]}) "
                     f"and mask file ({len(gate_ids)})")
        mali_runs_data.append({'label': label, 'years': run_years, 'flux': run_flux})

    for run, (total_color, gate_color) in zip(mali_runs_data, get_run_styles(len(mali_runs_data))):
        run['color'] = total_color
        run['gate_color'] = gate_color

    print(f"Reading GEUS observations from: {args.obs_dir}")
    region_of_gate, gate_D, gate_err = load_observations(
        args.obs_dir, args.gpkg, gate_ids, tolerance=args.match_tolerance)
    obs_years = gate_D.index.year + (gate_D.index.dayofyear - 1) / 365.0

    region_gate_indices = None
    for run in mali_runs_data:
        run_indices, run_region_totals, run_gis_total = aggregate_mali(
            run['flux'], gate_ids, region_of_gate)
        if region_gate_indices is None:
            region_gate_indices = run_indices
        run['region_totals'] = run_region_totals
        run['gis_total'] = run_gis_total

    region_gates, obs_region_D, obs_region_err, obs_gis_D, obs_gis_err = aggregate_obs(
        gate_D, gate_err, region_of_gate)
    obs_region_gate_D = {
        r: gate_D[gids] if gids else None for r, gids in region_gates.items()
    }
    obs_gis_gate_D = gate_D[list(region_of_gate.keys())]

    fig, axs = plt.subplots(2, 4, figsize=(18, 8), sharex=True)
    axs = axs.ravel()

    for ax, region in zip(axs[:7], REGIONS):
        indices = region_gate_indices[region]
        if not indices:
            ax.set_title(f'{region} (no gates in mask)')
            ax.axis('off')
            continue
        gate_word = 'gate' if len(indices) == 1 else 'gates'
        region_runs = [
            {**run, 'total': run['region_totals'][region]} for run in mali_runs_data
        ]
        plot_panel(ax, region_runs, indices,
                   obs_years, obs_region_gate_D[region], obs_region_D[region], obs_region_err[region],
                   f'{region} ({len(indices)} {gate_word})')

    gis_runs = [{**run, 'total': run['gis_total']} for run in mali_runs_data]
    plot_panel(axs[7], gis_runs, list(range(len(gate_ids))),
               obs_years, obs_gis_gate_D, obs_gis_D, obs_gis_err, 'Greenland total')

    handles, labels_ = axs[7].get_legend_handles_labels()
    fig.legend(handles, labels_, loc='lower center', ncol=min(len(mali_runs_data) + 1, 4))
    fig.suptitle('MALI vs. GEUS (Mankoff et al., 2020) ice discharge by region')
    fig.tight_layout(rect=(0, 0.05, 1, 1))

    fig.savefig(args.output, dpi=150)
    print(f"Wrote {args.output}")
    if args.show:
        plt.show()


if __name__ == '__main__':
    main()
