#!/usr/bin/env python
"""
Calibrate the damage threshold for MALI's crevasse-depth calving law.

For each Antarctic ice-shelf region (the ``ice_shelf_extended_*`` features in
geometric_features, tag ``ais_ice_shelf_calving_calibration``) choose the
damage threshold that best preserves the present calving front, given the
``damage`` field along that front in a MALI output file.

Calving-front cells are dynamic floating cells that either neighbor ice-free
ocean, or neighbor a non-dynamic floating cell that itself neighbors ice-free
ocean.  A cell calves when its damage >= threshold, so a threshold T is judged
by simulating the retreat it would cause starting from the observed front:
front cells with damage >= T are removed, which exposes their dynamic floating
neighbors, which are removed in turn if damage >= T, and so on.  (This captures
thin low-damage margins that protect a high-damage interior.)  The calibrated
threshold is the smallest T for which the floating-area loss stays within
--max-area-loss (default 5%), plus an optional --growth-margin (default 0).

Output: a netCDF file with ``damageCalvingThreshold(Time, nCells)``, constant
within each region.  Where regions overlap the largest value is used (least
retreat).  Cells outside all regions get the median of the calibrated values.

Requires numpy, scipy, netCDF4, matplotlib, shapely (>=2), pyproj and
geometric_features.  Example:
  calibrate_damage_threshold_for_cd_calving.py -i output.nc -m mesh.nc \\
      -o damageCalvingThreshold.nc
"""
import argparse
import glob
import json
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import netCDF4
import numpy as np
import pyproj
import shapely
from shapely.geometry import shape
from shapely.ops import transform

# MALI cellMask bits (mpas_li_mask.F)
MASK_ICE = 32
MASK_DYNAMIC = 2
MASK_FLOATING = 4


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('-i', '--input', required=True,
                   help='MALI output file with cellMask, damage, xCell, yCell, '
                        'dcEdge and cellsOnVertex')
    p.add_argument('-m', '--mesh', required=True,
                   help='MALI mesh file providing bedTopography (same cells as input)')
    p.add_argument('-o', '--output', default='damageCalvingThreshold.nc')
    p.add_argument('--region-tag', default='ais_ice_shelf_calving_calibration',
                   help='geometric_features tag selecting the landice regions to '
                        'calibrate (default ais_ice_shelf_calving_calibration)')
    p.add_argument('-g', '--geometric-features', default=None,
                   help='optional path to a local geometric_features checkout; if '
                        'given, regions are read from its geometric_data/landice/'
                        'region directory (names matching --region-glob) instead '
                        'of via the geometric_features package')
    p.add_argument('--region-glob', default='*',
                   help='region-name pattern used with --geometric-features')
    p.add_argument('--max-area-loss', type=float, default=0.05,
                   help='floating-area fraction that must stay strictly below this '
                        '(default 0.05)')
    p.add_argument('--growth-margin', type=float, default=0.0,
                   help='damage added to the minimal threshold to favor minor '
                        'growth (default 0)')
    p.add_argument('--crs', default='EPSG:3031',
                   help='projection of xCell/yCell (default EPSG:3031)')
    p.add_argument('--time-index', type=int, default=0)
    p.add_argument('--fig-dir', default='figures',
                   help='directory for per-region figures (default figures)')
    return p.parse_args()


def neighbors_from_vertices(cells_on_vertex, ncells):
    """Return a sparse symmetric adjacency (CSR) of cells sharing a vertex pair
    (cellsOnVertex is 1-based, 0 = missing)."""
    from scipy.sparse import coo_matrix
    cov = cells_on_vertex - 1
    pairs = []
    for a, b in ((0, 1), (1, 2), (2, 0)):
        ok = (cov[:, a] >= 0) & (cov[:, b] >= 0) & (cov[:, a] != cov[:, b])
        pairs.append(np.stack([cov[ok, a], cov[ok, b]], axis=1))
    pairs = np.concatenate(pairs)
    adj = coo_matrix((np.ones(len(pairs), dtype=np.int8), (pairs[:, 0], pairs[:, 1])),
                     shape=(ncells, ncells)).tocsr()
    adj.data[:] = 1
    adj = adj.maximum(adj.T)  # duplicates collapse to 1
    return adj


def load_regions(gf_path, pattern, tag, crs):
    """Return {name: shapely geometry in projected coordinates}."""
    to_proj = pyproj.Transformer.from_crs('EPSG:4326', crs, always_xy=True).transform
    features = {}
    if gf_path is None:
        from geometric_features import GeometricFeatures
        fc = GeometricFeatures().read('landice', 'region', tags=[tag])
        for feat in fc.features:
            features[feat['properties']['name']] = [feat]
    else:
        files = sorted(glob.glob(os.path.join(gf_path, 'geometric_data', 'landice',
                                              'region', pattern, 'region.geojson')))
        for f in files:
            with open(f) as fp:
                gj = json.load(fp)
            features[os.path.basename(os.path.dirname(f))] = \
                gj['features'] if gj['type'] == 'FeatureCollection' else [gj]
    if not features:
        raise SystemExit('No regions found; check --region-tag / --geometric-features')
    regions = {}
    for name in sorted(features):
        geoms = [shape(ft['geometry']) for ft in features[name]]
        geom = shapely.union_all([shapely.segmentize(g, 0.05) for g in geoms])
        regions[name] = transform(to_proj, geom)
    return regions


def plot_region(name, geom, thr, loss, removed, in_reg, front, ice, damage, x, y, dc, fig_dir):
    """Map of damage on ice cells near a region, with the chosen-threshold contour."""
    # zoom to the region's front cells (the region itself extends far over ocean)
    fz = np.where(front & in_reg)[0]
    x0, x1, y0, y1 = x[fz].min(), x[fz].max(), y[fz].min(), y[fz].max()
    pad = 0.1 * max(x1 - x0, y1 - y0) + 20 * dc
    box = ice & (x > x0 - pad) & (x < x1 + pad) & (y > y0 - pad) & (y < y1 + pad)
    idx = np.where(box)[0]
    fig, ax = plt.subplots(figsize=(8, 7))
    if len(idx) > 3:
        tri = mtri.Triangulation(x[idx], y[idx])
        xt, yt = x[idx][tri.triangles], y[idx][tri.triangles]
        edge = np.max(np.hypot(xt - np.roll(xt, 1, axis=1), yt - np.roll(yt, 1, axis=1)), axis=1)
        tri.set_mask(edge > 1.8 * dc)
        d = damage[idx]
        pc = ax.tripcolor(tri, d, vmin=0, vmax=1, cmap='viridis', shading='gouraud')
        fig.colorbar(pc, ax=ax, label='damage', shrink=0.8)
        if d.min() < thr < d.max():
            ax.tricontour(tri, d, levels=[thr], colors='red', linewidths=1.5)
    fi = np.where(front & in_reg)[0]
    ri = np.where(removed)[0]
    ax.scatter(x[ri], y[ri], s=6, marker='s', c='red', alpha=0.6, linewidths=0)
    ax.scatter(x[fi], y[fi], s=10, c='none', edgecolors='white', linewidths=0.6)
    gs = [geom] if geom.geom_type == 'Polygon' else list(geom.geoms)
    for g in gs:
        ax.plot(*g.exterior.xy, 'k--', lw=0.8)
    ax.set_autoscale_on(False)
    ax.plot([], [], 'r-', label=f'damage = {thr:.3f} (chosen threshold)')
    ax.scatter([], [], s=10, c='none', edgecolors='k', label='calving-front cells (outlined)')
    ax.scatter([], [], s=12, marker='s', c='red', alpha=0.6,
               label=f'cells removed (connected to terminus): {100 * loss:.2f}% of floating area')
    ax.legend(loc='lower left', fontsize=8)
    ax.set_aspect('equal')
    ax.set_xlim(x0 - pad, x1 + pad)
    ax.set_ylim(y0 - pad, y1 + pad)
    ax.set_title(f'{name}\nchosen damageCalvingThreshold = {thr:.3f}  (area removed {100 * loss:.2f}%)')
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    fig.tight_layout()
    fig.savefig(os.path.join(fig_dir, f'{name}.png'), dpi=150)
    plt.close(fig)


def simulate_retreat(thr, dyn_float, damage, area, adj, seeds_exposed, region):
    """Remove region dynamic-floating cells with damage >= thr, iteratively
    starting from cells exposed to ocean.  Returns the boolean removed mask."""
    removable = region & dyn_float & (damage >= thr)
    removed = np.zeros_like(region)
    exposed = seeds_exposed.copy()
    frontier = removable & exposed
    while frontier.any():
        removed |= frontier
        # neighbors of newly removed cells become exposed
        newly = np.asarray(adj[frontier].sum(axis=0)).ravel() > 0
        frontier = removable & newly & ~removed
    return removed


def main():
    args = parse_args()
    ds = netCDF4.Dataset(args.input)
    t = args.time_index
    ncells = len(ds.dimensions['nCells'])
    mask = ds['cellMask'][t, :].astype(int)
    damage = np.asarray(ds['damage'][t, :], dtype=float)
    with netCDF4.Dataset(args.mesh) as mesh:
        bed = np.asarray(mesh['bedTopography'][0, :], dtype=float)
        assert bed.shape[0] == ncells, 'mesh file and input have different nCells'
    x = np.asarray(ds['xCell'][:])
    y = np.asarray(ds['yCell'][:])
    dc = float(np.median(ds['dcEdge'][:]))
    # hexagonal Voronoi cell area from the median cell spacing (no areaCell in file)
    area = np.full(ncells, np.sqrt(3.0) / 2.0 * dc ** 2)
    adj = neighbors_from_vertices(np.asarray(ds['cellsOnVertex'][:]), ncells)

    ice = (mask & MASK_ICE) > 0
    floating = ice & ((mask & MASK_FLOATING) > 0)
    dyn_float = floating & ((mask & MASK_DYNAMIC) > 0)
    nondyn_float = floating & ~dyn_float
    # ice-free ocean: no ice and bed below sea level
    ocean = ~ice & (bed < 0.0)

    def has_neighbor(flag):
        return np.asarray(adj @ flag.astype(np.int8)).ravel() > 0

    nbr_ocean = has_neighbor(ocean)
    nondyn_near_ocean = nondyn_float & nbr_ocean
    front = dyn_float & (nbr_ocean | has_neighbor(nondyn_near_ocean))
    # cells whose removal/exposure starts retreat: next to ocean or to
    # (non-dynamic) floating ice that borders ocean
    exposed0 = nbr_ocean | has_neighbor(nondyn_near_ocean)
    print(f'{front.sum()} calving-front cells, {dyn_float.sum()} dynamic floating cells, '
          f'cell area {area[0] / 1e6:.2f} km^2')

    regions = load_regions(args.geometric_features, args.region_glob, args.region_tag,
                           args.crs)
    os.makedirs(args.fig_dir, exist_ok=True)
    result = {}
    hdr = f"{'region':45s} {'front':>6s} {'d_p10':>6s} {'d_p50':>6s} {'d_p90':>6s} " \
          f"{'thr':>6s} {'loss%':>6s}"
    print(hdr)
    for name, geom in regions.items():
        shapely.prepare(geom)
        in_reg = shapely.contains_xy(geom, x, y)
        fl_area = area[in_reg & floating].sum()
        fr = in_reg & front
        if fl_area == 0 or fr.sum() == 0:
            print(f'{name:45s} no floating/front cells -> default later')
            result[name] = (in_reg, np.nan)
            continue
        cand = np.unique(np.concatenate([damage[fr], [0.0, 1.0 + 1e-9]]))
        # loss is non-increasing in threshold; binary search the smallest feasible one
        lo, hi = 0, len(cand) - 1
        while lo < hi:
            mid = (lo + hi) // 2
            loss = area[simulate_retreat(cand[mid], dyn_float, damage, area, adj,
                                         exposed0, in_reg)].sum() / fl_area
            if loss < args.max_area_loss:
                hi = mid
            else:
                lo = mid + 1
        thr_min = cand[lo]
        thr = min(thr_min + args.growth_margin, 1.0)
        removed = simulate_retreat(thr, dyn_float, damage, area, adj, exposed0, in_reg)
        loss = area[removed].sum() / fl_area
        p10, p50, p90 = np.percentile(damage[fr], [10, 50, 90])
        print(f'{name:45s} {fr.sum():6d} {p10:6.3f} {p50:6.3f} {p90:6.3f} '
              f'{thr:6.3f} {100 * loss:6.2f}')
        result[name] = (in_reg, thr)
        plot_region(name, geom, thr, loss, removed, in_reg, front, ice, damage, x, y, dc,
                    args.fig_dir)

    vals = np.array([v for _, v in result.values() if np.isfinite(v)])
    default = float(np.median(vals))
    print(f'Default (median over regions) for cells outside regions / uncalibrated: {default:.3f}')

    out = np.full(ncells, np.nan)
    for in_reg, thr in result.values():
        v = thr if np.isfinite(thr) else default
        out[in_reg] = np.fmax(out[in_reg], v)  # overlap -> largest value
    out[np.isnan(out)] = default

    with netCDF4.Dataset(args.output, 'w') as o:
        o.createDimension('Time', None)
        o.createDimension('nCells', ncells)
        var = o.createVariable('damageCalvingThreshold', 'f8', ('Time', 'nCells'))
        var.long_name = 'damage threshold for crevasse-depth calving'
        var[0, :] = out
        for nm in ('xCell', 'yCell', 'zCell'):
            v = o.createVariable(nm, 'f8', ('nCells',))
            v[:] = ds[nm][:]
        o.createDimension('nVertices', len(ds.dimensions['nVertices']))
        o.createDimension('vertexDegree', len(ds.dimensions['vertexDegree']))
        v = o.createVariable('cellsOnVertex', 'i4', ('nVertices', 'vertexDegree'))
        v[:] = ds['cellsOnVertex'][:]
        o.history = f'Created by calibrate_damage_threshold.py from {args.input}; ' \
                    f'max_area_loss={args.max_area_loss}, growth_margin={args.growth_margin}'
    print(f'Wrote {args.output}')


if __name__ == '__main__':
    main()
