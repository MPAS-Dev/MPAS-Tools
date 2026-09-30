#!/usr/bin/env python3
"""
Generate MALI flux-gate mask files from GeoPackage gate geometry.

This script converts ice-discharge gate geometry (e.g. from Mankoff et al., GEUS
https://doi.org/10.22008/promice/data/ice_discharge/gates/v02) into an MPAS/MALI
edge mask file for use with the flux-gates analysis member.

MALI's fluxGatesInput stream expects a NetCDF file with integer variables
`fluxGateEdgeMasks` and `fluxGateEdgeMasksSigns` dimensioned (nEdges, nFluxGates).
The script produces this from a GeoPackage by:
  1. Reading per-pixel gate coordinates from the gpkg table (via sqlite3)
  2. Transforming coordinates to EPSG:4326 and building LineString transects
  3. Calling mpas_tools.mesh.mask.compute_mpas_transect_masks to compute edge masks
  4. Culling ambiguous edges (mask==1 & sign==0)
  5. Renaming transect variables → fluxGate variables and writing the output

Dependencies (all in the compass pixi environment):
  - Python stdlib: sqlite3
  - pyproj, shapely (>=2.0), geometric_features, mpas_tools, xarray, netcdf4

Example usage:
    ./make_flux_gate_masks.py \\
        -g gates.gpkg \\
        -m Humboldt_1to10km.nc \\
        -o Humboldt_flux_gates_mask.nc \\
        --gates 28 \\
        --geojson Humboldt_flux_gates.geojson

Created: 2026-09-30
@author: Trevor Hillebrand
"""

import argparse
import os
import sqlite3
import sys
import tempfile
from datetime import datetime

import numpy as np
import xarray as xr
from geometric_features import FeatureCollection, read_feature_collection
from mpas_tools.io import write_netcdf
from mpas_tools.logging import LoggingContext
from mpas_tools.mesh.mask import compute_mpas_transect_masks
from mpas_tools.parallel import create_pool
from pyproj import Transformer
from shapely import LineString
from shapely.geometry import mapping


def read_gate_coordinates(gpkg_path, gate_ids, table="gates_final", gate_field="gate",
                          x_field="x", y_field="y"):
    """
    Extract ordered (x, y) coordinates for each gate from a GeoPackage.

    Parameters
    ----------
    gpkg_path : str
        Path to the GeoPackage file (SQLite database).
    gate_ids : list of int
        Gate IDs to extract.
    table : str
        Feature table name in the gpkg.
    gate_field : str
        Column name identifying gate IDs.
    x_field, y_field : str
        Column names holding projected coordinates.

    Returns
    -------
    dict
        {gate_id: [(x1, y1), (x2, y2), ...]}
    """
    conn = sqlite3.connect(gpkg_path)
    cursor = conn.cursor()

    # Validate columns exist
    cursor.execute(f"PRAGMA table_info({table})")
    columns = [row[1] for row in cursor.fetchall()]
    missing = [f for f in [gate_field, x_field, y_field] if f not in columns]
    if missing:
        conn.close()
        sys.exit(f"ERROR: columns {missing} not found in table '{table}'. "
                 f"Available columns: {columns}")

    # Validate gate IDs exist
    cursor.execute(f"SELECT DISTINCT {gate_field} FROM {table}")
    available_gates = sorted([row[0] for row in cursor.fetchall()])
    missing_gates = [gid for gid in gate_ids if gid not in available_gates]
    if missing_gates:
        conn.close()
        sys.exit(f"ERROR: gate IDs {missing_gates} not found in table '{table}'. "
                 f"Available gates: {available_gates}")

    gate_coords = {}
    skipped = []
    for gate_id in gate_ids:
        cursor.execute(
            f"SELECT {x_field}, {y_field} FROM {table} "
            f"WHERE {gate_field}=? ORDER BY fid",
            (gate_id,)
        )
        coords = cursor.fetchall()
        if len(coords) < 2:
            skipped.append((gate_id, len(coords)))
        else:
            gate_coords[gate_id] = coords

    conn.close()

    if skipped:
        print(f"WARNING: Skipped {len(skipped)} gate(s) with insufficient points:")
        for gate_id, n_points in skipped:
            print(f"  Gate {gate_id}: {n_points} point(s) (need at least 2)")

    if not gate_coords:
        sys.exit("ERROR: No valid gates found (all had fewer than 2 points)")

    return gate_coords


def make_geojson_features(gate_coords, source_epsg, name_prefix="gate"):
    """
    Build GeoJSON features for geometric_features from gate coordinates.

    Parameters
    ----------
    gate_coords : dict
        {gate_id: [(x1, y1), ...]} in source CRS.
    source_epsg : int
        EPSG code of the source coordinates.
    name_prefix : str
        Prefix for feature names (e.g. "gate28").

    Returns
    -------
    list of dict
        GeoJSON Feature dicts ready for geometric_features.FeatureCollection.
    """
    # pyproj Transformer from source → EPSG:4326
    # Default axis order for 4326 is (lat, lon), so transform returns (lat, lon)
    transformer = Transformer.from_crs(source_epsg, 4326, always_xy=False)

    features = []
    for gate_id, coords in sorted(gate_coords.items()):
        # Transform each point: trans(x, y) → (lat, lon)
        # GeoJSON LineString coordinates are [lon, lat], so swap
        lonlat_coords = []
        for x, y in coords:
            lat, lon = transformer.transform(x, y)
            lonlat_coords.append([lon, lat])

        line = LineString(lonlat_coords)
        geometry = mapping(line)

        properties = {
            "name": f"{name_prefix}{gate_id}",
            "component": "landice",
            "object": "transect",
            "tags": "",
            "constituents": "",
        }

        features.append({
            "type": "Feature",
            "geometry": geometry,
            "properties": properties,
        })

    return features


def cull_ambiguous_edges(ds_masks):
    """
    Zero out transectEdgeMasks where mask==1 and sign==0 (ambiguous edges).

    Parameters
    ----------
    ds_masks : xarray.Dataset
        Dataset with transectEdgeMasks and transectEdgeMaskSigns.

    Returns
    -------
    int
        Number of edges culled.
    """
    masks = ds_masks["transectEdgeMasks"]
    signs = ds_masks["transectEdgeMaskSigns"]
    ambiguous = (masks == 1) & (signs == 0)
    n_culled = int(ambiguous.sum())
    ds_masks["transectEdgeMasks"] = masks.where(~ambiguous, 0)
    return n_culled


def convert_to_compact_format(ds_masks):
    """
    Convert sparse 2D transect masks to compact 1D per-edge gate IDs.

    MALI's fluxGatesInput stream reads:
      - fluxGateEdgeID(nEdges) - gate number (1..nFluxGates) per edge, 0 if none
      - fluxGateEdgeSign(nEdges) - flux sign (-1/+1) per edge, 0 if none
      - fluxGateNames(nFluxGates) - gate names

    mpas_tools produces:
      - transectEdgeMasks(nEdges, nTransects) - sparse 2D one-hot
      - transectEdgeMaskSigns(nEdges, nTransects) - sparse 2D signs

    For each edge, find which gate (if any) it belongs to and the corresponding sign.
    Warns if an edge belongs to multiple gates (takes the first).

    Parameters
    ----------
    ds_masks : xarray.Dataset
        Dataset with transectEdgeMasks, transectEdgeMaskSigns, transectNames.

    Returns
    -------
    xarray.Dataset
        Compact format dataset with fluxGateEdgeID, fluxGateEdgeSign, fluxGateNames.
    """
    masks = ds_masks["transectEdgeMasks"].values  # (nEdges, nTransects)
    signs = ds_masks["transectEdgeMaskSigns"].values
    n_edges, n_gates = masks.shape

    # Find which gate each edge belongs to (argmax gives 0-based index)
    # edges with no gates will have mask sum = 0
    edge_gate_id = np.zeros(n_edges, dtype=np.int32)
    edge_gate_sign = np.zeros(n_edges, dtype=np.int32)

    n_overlap = 0
    for i_edge in range(n_edges):
        gate_indices = np.where(masks[i_edge, :] == 1)[0]
        if len(gate_indices) > 0:
            if len(gate_indices) > 1:
                n_overlap += 1
            # Use first gate (1-based indexing for MALI)
            i_gate = gate_indices[0]
            edge_gate_id[i_edge] = i_gate + 1  # 1-based
            edge_gate_sign[i_edge] = signs[i_edge, i_gate]

    if n_overlap > 0:
        print(f"WARNING: {n_overlap} edges belong to multiple gates (using first gate)")

    # Build compact dataset with proper dimensions
    ds_compact = xr.Dataset(
        {
            "fluxGateEdgeID": (("nEdges",), edge_gate_id),
            "fluxGateEdgeSign": (("nEdges",), edge_gate_sign),
        },
        coords={"nEdges": np.arange(n_edges)},
    )

    # Add gate names if present
    if "transectNames" in ds_masks:
        gate_names = ds_masks["transectNames"].values
        ds_compact["fluxGateNames"] = (("nFluxGates",), gate_names)
        ds_compact.coords["nFluxGates"] = np.arange(n_gates)

    return ds_compact


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("-g", "--gpkg", required=True,
                        help="Input GeoPackage (or fiona-readable) file with gate geometry")
    parser.add_argument("-m", "--mesh", required=True,
                        help="MPAS/MALI mesh file")
    parser.add_argument("-o", "--output", default="fluxGateMasks.nc",
                        help="Output mask file for MALI fluxGatesInput stream "
                             "(default: fluxGateMasks.nc)")
    parser.add_argument("--gates", nargs="+", type=int,
                        help="Gate IDs to include (space-separated)")
    parser.add_argument("--all-gates", action="store_true",
                        help="Process all gates found in the GeoPackage")
    parser.add_argument("--table", default="gates_final",
                        help="GeoPackage feature table name (default: gates_final)")
    parser.add_argument("--gate-field", default="gate",
                        help="Column name for gate IDs (default: gate)")
    parser.add_argument("--x-field", default="x",
                        help="Column name for x coordinates (default: x)")
    parser.add_argument("--y-field", default="y",
                        help="Column name for y coordinates (default: y)")
    parser.add_argument("--source-epsg", type=int, default=3413,
                        help="EPSG code of the gpkg coordinates (default: 3413)")
    parser.add_argument("--name-prefix", default="gate",
                        help="Prefix for feature names (default: gate)")
    parser.add_argument("--geojson",
                        help="Path to save intermediate geojson (optional; temp file if omitted)")
    parser.add_argument("--subdivision", type=float, default=1000.0,
                        help="Transect subdivision resolution in meters (default: 1000)")
    parser.add_argument("--add-edge-sign", dest="add_edge_sign", action="store_true",
                        default=True, help="Compute edge signs (default: True)")
    parser.add_argument("--no-edge-sign", dest="add_edge_sign", action="store_false",
                        help="Skip edge sign computation (not recommended for flux gates)")
    parser.add_argument("--process-count", type=int,
                        help="Number of parallel processes (default: all available cores)")
    parser.add_argument("--show-progress", action="store_true",
                        help="Show progress bar during mask computation")
    parser.add_argument("--overwrite", action="store_true",
                        help="Overwrite existing output file")

    args = parser.parse_args(argv)

    # Validate inputs
    if not args.gates and not args.all_gates:
        parser.error("Either --gates or --all-gates must be specified")
    if args.gates and args.all_gates:
        parser.error("Cannot specify both --gates and --all-gates")
    if not os.path.exists(args.gpkg):
        parser.error(f"GeoPackage file not found: {args.gpkg}")
    if not os.path.exists(args.mesh):
        parser.error(f"Mesh file not found: {args.mesh}")
    if os.path.exists(args.output) and not args.overwrite:
        parser.error(f"Output file exists: {args.output} (use --overwrite to replace)")

    # Determine which gates to process
    if args.all_gates:
        print(f"Querying all gates from {args.gpkg} ...")
        conn = sqlite3.connect(args.gpkg)
        cursor = conn.cursor()
        cursor.execute(f"SELECT DISTINCT {args.gate_field} FROM {args.table} ORDER BY {args.gate_field}")
        gate_ids = [row[0] for row in cursor.fetchall()]
        conn.close()
        print(f"  Found {len(gate_ids)} gates")
    else:
        gate_ids = args.gates

    print(f"Reading gate coordinates from {args.gpkg} ...")
    gate_coords = read_gate_coordinates(
        args.gpkg, gate_ids, args.table, args.gate_field, args.x_field, args.y_field
    )
    for gate_id, coords in gate_coords.items():
        print(f"  Gate {gate_id}: {len(coords)} points")

    print(f"Building GeoJSON features (EPSG:{args.source_epsg} → 4326) ...")
    features = make_geojson_features(gate_coords, args.source_epsg, args.name_prefix)
    fc = FeatureCollection(features)

    # Write geojson (temp or user-specified)
    if args.geojson:
        geojson_path = args.geojson
        fc.to_geojson(geojson_path)
        print(f"  Wrote intermediate geojson: {geojson_path}")
    else:
        # Use temp file
        fd, geojson_path = tempfile.mkstemp(suffix=".geojson", text=True)
        os.close(fd)
        fc.to_geojson(geojson_path)

    try:
        print(f"Loading mesh: {args.mesh} ...")
        ds_mesh = xr.open_dataset(args.mesh, decode_cf=False, decode_times=False)

        print("Computing transect edge masks ...")
        # Import constants here (some mpas_tools versions have it at module level)
        try:
            from mpas_tools.cime.constants import constants
            earth_radius = constants["SHR_CONST_REARTH"]
        except ImportError:
            # Fallback for older mpas_tools
            earth_radius = 6.37122e6  # meters

        pool = create_pool(process_count=args.process_count, method="forkserver")
        fc_mask = read_feature_collection(geojson_path)

        with LoggingContext("make_flux_gate_masks") as logger:
            ds_masks = compute_mpas_transect_masks(
                dsMesh=ds_mesh,
                fcMask=fc_mask,
                earthRadius=earth_radius,
                maskTypes=("edge",),
                logger=logger,
                pool=pool,
                chunkSize=1000,
                showProgress=args.show_progress,
                subdivisionResolution=args.subdivision,
                addEdgeSign=args.add_edge_sign,
            )

        print("Culling ambiguous edges (mask==1 & sign==0) ...")
        n_culled = cull_ambiguous_edges(ds_masks)
        print(f"  {n_culled} edges culled from the transect mask")

        print("Converting to compact format (per-edge gate IDs) ...")
        ds_mali = convert_to_compact_format(ds_masks)

        print(f"Writing output: {args.output} ...")
        with LoggingContext("write_netcdf") as logger:
            write_netcdf(ds_mali, args.output, logger=logger)

        print(f"\nSuccess! Output file: {args.output}")
        print(f"  Variables: fluxGateEdgeID, fluxGateEdgeSign, fluxGateNames")
        print(f"  Dimensions: nEdges={ds_mali.sizes['nEdges']}, "
              f"nFluxGates={ds_mali.sizes['nFluxGates']}")

        # Report file size savings
        n_edges = ds_mali.sizes['nEdges']
        n_gates = ds_mali.sizes['nFluxGates']
        compact_size_mb = (2 * n_edges * 4) / (1024**2)  # 2 int32 arrays
        onehot_size_mb = (2 * n_edges * n_gates * 4) / (1024**2)  # what it would have been
        print(f"  File size: ~{compact_size_mb:.1f} MB (vs {onehot_size_mb:.1f} MB one-hot, "
              f"{onehot_size_mb/compact_size_mb:.0f}x reduction)")

    finally:
        # Clean up temp geojson if we created it
        if not args.geojson and os.path.exists(geojson_path):
            os.remove(geojson_path)


if __name__ == "__main__":
    main()
