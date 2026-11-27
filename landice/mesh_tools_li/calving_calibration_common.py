"""
Helpers shared by the MALI calving-law calibration scripts
(``calibrate_damage_threshold_for_cd_calving.py`` and
``calibrate_eigencalving_parameter.py``): mesh adjacency, cellMask decoding and
calving-front identification, region loading, combining per-region values into
a cell field, writing the output file, and the base map used in the figures.
"""
import glob
import json
import os

import matplotlib.tri as mtri
import netCDF4
import numpy as np
import pyproj
import shapely
from matplotlib.colors import ListedColormap
from shapely.geometry import shape
from shapely.ops import transform

# MALI cellMask bits (mpas_li_mask.F)
MASK_ICE = 32
MASK_DYNAMIC = 2
MASK_FLOATING = 4
MASK_MARGIN = 8


def add_region_arguments(p):
    """Add the region-selection and projection arguments shared by the scripts."""
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
    p.add_argument('--crs', default='EPSG:3031',
                   help='projection of xCell/yCell (default EPSG:3031)')


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


def estimate_cell_area(dc, ncells):
    """Hexagonal Voronoi cell area from the median cell spacing dc (m), for files
    without areaCell."""
    return np.full(ncells, np.sqrt(3.0) / 2.0 * dc ** 2)


def read_cell_area(ds, mesh_path, dc, t=0):
    """areaCell from the input file or the mesh file if either has it, else the
    hexagonal estimate from the median cell spacing dc."""
    ncells = len(ds.dimensions['nCells'])
    if 'areaCell' in ds.variables:
        return np.asarray(ds['areaCell'][:], dtype=float)
    if mesh_path is not None:
        with netCDF4.Dataset(mesh_path) as mesh:
            if 'areaCell' in mesh.variables:
                return np.asarray(mesh['areaCell'][:], dtype=float)
    return estimate_cell_area(dc, ncells)


def read_bed_topography(ds, mesh_path, t=0):
    """bedTopography from the input file if it has one, else from the mesh file."""
    ncells = len(ds.dimensions['nCells'])
    if 'bedTopography' in ds.variables:
        return np.asarray(ds['bedTopography'][t, :], dtype=float)
    if mesh_path is None:
        raise SystemExit('bedTopography is not in the input file; give a mesh file '
                         'that has it with --mesh')
    with netCDF4.Dataset(mesh_path) as mesh:
        bed = np.asarray(mesh['bedTopography'][0, :], dtype=float)
    assert bed.shape[0] == ncells, 'mesh file and input have different nCells'
    return bed


def identify_front_cells(mask, bed, adj):
    """Decode cellMask and find the calving front.

    The front is made of dynamic floating cells that neighbor ice-free ocean
    (no ice and bed below sea level), or neighbor a non-dynamic floating cell
    that itself neighbors ice-free ocean.  Grounded marine margin cells, which
    MALI's rate-based calving laws also act on, are those dynamic grounded cells
    with bed below sea level that neighbor ice-free ocean or an inactive
    (non-dynamic) floating margin cell that borders ocean.  Returns a dict of
    boolean arrays: ice, floating, dyn_float, front, grounded_margin and
    exposed0 (cells next to ocean, or next to non-dynamic floating ice bordering
    ocean, where a retreat would start).
    """
    ice = (mask & MASK_ICE) > 0
    floating = ice & ((mask & MASK_FLOATING) > 0)
    dynamic = ice & ((mask & MASK_DYNAMIC) > 0)
    dyn_float = floating & dynamic
    nondyn_float = floating & ~dyn_float
    grounded = ice & ~floating
    ocean = ~ice & (bed < 0.0)

    def has_neighbor(flag):
        return np.asarray(adj @ flag.astype(np.int8)).ravel() > 0

    nbr_ocean = has_neighbor(ocean)
    nondyn_near_ocean = nondyn_float & nbr_ocean
    exposed0 = nbr_ocean | has_neighbor(nondyn_near_ocean)
    front = dyn_float & exposed0
    inactive_margin = nondyn_float & ((mask & MASK_MARGIN) > 0)
    grounded_margin = (grounded & dynamic & (bed < 0.0)
                       & (nbr_ocean | has_neighbor(inactive_margin & nbr_ocean)))
    return dict(ice=ice, floating=floating, dyn_float=dyn_float, front=front,
                grounded_margin=grounded_margin, exposed0=exposed0)


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


def region_mask(geom, x, y):
    """Boolean array of cell centers inside a (projected) region geometry."""
    shapely.prepare(geom)
    return shapely.contains_xy(geom, x, y)


def combine_region_values(results, ncells, overlap='max'):
    """Combine per-region constants into one cell field.

    results: {name: (in_region_mask, value)}; a NaN value means the region could
    not be calibrated.  Where regions overlap the largest ('max') or smallest
    ('min') value is used; cells outside all regions and uncalibrated regions get
    the median of the calibrated values.  Returns (field, default).
    """
    vals = np.array([v for _, v in results.values() if np.isfinite(v)])
    default = float(np.median(vals))
    pick = np.fmax if overlap == 'max' else np.fmin
    out = np.full(ncells, np.nan)
    for in_reg, value in results.values():
        v = value if np.isfinite(value) else default
        out[in_reg] = pick(out[in_reg], v)
    out[np.isnan(out)] = default
    return out, default


def write_calibration_file(path, varname, long_name, units, values, ds, history):
    """Write the calibrated cell field plus the coordinates and connectivity
    needed to visualize it."""
    ncells = len(values)
    with netCDF4.Dataset(path, 'w') as o:
        o.createDimension('Time', None)
        o.createDimension('nCells', ncells)
        var = o.createVariable(varname, 'f8', ('Time', 'nCells'))
        var.long_name = long_name
        if units:
            var.units = units
        var[0, :] = values
        for nm in ('xCell', 'yCell', 'zCell'):
            v = o.createVariable(nm, 'f8', ('nCells',))
            v[:] = ds[nm][:]
        o.createDimension('nVertices', len(ds.dimensions['nVertices']))
        o.createDimension('vertexDegree', len(ds.dimensions['vertexDegree']))
        v = o.createVariable('cellsOnVertex', 'i4', ('nVertices', 'vertexDegree'))
        v[:] = ds['cellsOnVertex'][:]
        o.history = history


def draw_base_map(ax, geom, in_reg, front, ice, floating, x, y, dc):
    """Zoom to the region's front cells and draw grounded ice in gray.

    Returns (tri, idx, extent): tri is a Triangulation of the ice cells idx in
    the view with triangles that touch non-floating cells (or span a gap)
    masked, ready for coloring floating ice; extent = (x0, x1, y0, y1, pad).
    """
    fz = np.where(front & in_reg)[0]
    x0, x1, y0, y1 = x[fz].min(), x[fz].max(), y[fz].min(), y[fz].max()
    pad = 0.1 * max(x1 - x0, y1 - y0) + 20 * dc
    box = ice & (x > x0 - pad) & (x < x1 + pad) & (y > y0 - pad) & (y < y1 + pad)
    idx = np.where(box)[0]
    tri = None
    if len(idx) > 3:
        tri = mtri.Triangulation(x[idx], y[idx])
        xt, yt = x[idx][tri.triangles], y[idx][tri.triangles]
        edge = np.max(np.hypot(xt - np.roll(xt, 1, axis=1), yt - np.roll(yt, 1, axis=1)),
                      axis=1)
        long_edge = edge > 1.8 * dc
        # grounded ice in gray underneath; the caller colors floating ice
        # (separate Triangulation objects: a mask is read at draw time)
        tri_all = mtri.Triangulation(x[idx], y[idx], triangles=tri.triangles,
                                     mask=long_edge)
        ax.tripcolor(tri_all, np.zeros(len(idx)), vmin=0, vmax=1,
                     cmap=ListedColormap(['0.6']), shading='gouraud')
        tri.set_mask(long_edge | ~np.all(floating[idx][tri.triangles], axis=1))
    return tri, idx, (x0, x1, y0, y1, pad)


def finish_map(ax, geom, extent):
    """Draw the region outline and set the view to the extent from draw_base_map."""
    x0, x1, y0, y1, pad = extent
    gs = [geom] if geom.geom_type == 'Polygon' else list(geom.geoms)
    for g in gs:
        ax.plot(*g.exterior.xy, 'k--', lw=0.8)
    ax.set_autoscale_on(False)
    ax.set_aspect('equal')
    ax.set_xlim(x0 - pad, x1 + pad)
    ax.set_ylim(y0 - pad, y1 + pad)
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
