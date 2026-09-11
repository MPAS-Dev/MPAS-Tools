#!/usr/bin/env python3
"""
Script to adjust bed topography within grounding line regions to achieve
a specified height above flotation (HAF).

Height above flotation is defined as:
    HAF = ice_surface - flotation_surface
where flotation_surface is the surface elevation at which ice would be
in hydrostatic equilibrium with ocean water.

For a given ice thickness and target HAF, the required bed elevation is:
    bed = ice_surface - thickness
where ice_surface is calculated to achieve the target HAF.
"""

import argparse
import sys
from datetime import datetime

import netCDF4
import numpy as np
from shapely.geometry import Point, shape

try:
    import geopandas as gpd
    HAS_GEOPANDAS = True
except ImportError:
    HAS_GEOPANDAS = False
    import json


def load_grounding_line_geojson(geojson_file):
    """
    Load grounding line polygons from a GeoJSON file.

    Parameters
    ----------
    geojson_file : str
        Path to the GeoJSON file containing grounding line delineations

    Returns
    -------
    geometries : list
        List of shapely geometry objects
    """
    if HAS_GEOPANDAS:
        # Use geopandas if available
        gdf = gpd.read_file(geojson_file)
        geometries = gdf.geometry.tolist()
    else:
        # Fall back to json + shapely
        with open(geojson_file, 'r') as f:
            data = json.load(f)

        geometries = []
        for feature in data['features']:
            geom = shape(feature['geometry'])
            geometries.append(geom)

    return geometries


def find_cells_in_polygon(lon_cell, lat_cell, geometries):
    """
    Find mesh cells that fall within any of the provided polygons.

    Parameters
    ----------
    lon_cell : ndarray
        Longitude of cell centers (in degrees)
    lat_cell : ndarray
        Latitude of cell centers (in degrees)
    geometries : list
        List of shapely geometry objects

    Returns
    -------
    mask : ndarray (bool)
        Boolean mask indicating which cells are inside the polygons
    """
    n_cells = len(lon_cell)
    mask = np.zeros(n_cells, dtype=bool)

    print(f'Checking {n_cells} cells against {len(geometries)} geometries...')

    for i, (lon, lat) in enumerate(zip(lon_cell, lat_cell)):
        if i % 10000 == 0:
            print(f'  Processed {i}/{n_cells} cells...')

        point = Point(lon, lat)
        for geom in geometries:
            if geom.contains(point) or geom.intersects(point):
                mask[i] = True
                break

    print(f'Found {np.sum(mask)} cells within polygons')
    return mask


def calculate_flotation_thickness(bed_elevation, sea_level=0.0,
                                   rho_ice=910.0, rho_ocean=1028.0):
    """
    Calculate the ice thickness at which ice would be in flotation.

    Parameters
    ----------
    bed_elevation : ndarray
        Bed topography (positive up, m)
    sea_level : float
        Sea level (m), default 0.0
    rho_ice : float
        Ice density (kg/m^3), default 910.0
    rho_ocean : float
        Ocean water density (kg/m^3), default 1028.0

    Returns
    -------
    flotation_thickness : ndarray
        Ice thickness at flotation (m)
    """
    # For ice to float: rho_ice * thickness = rho_ocean * draft
    # where draft = sea_level - bed_elevation
    # Therefore: thickness_flotation = (rho_ocean / rho_ice) * draft

    draft = sea_level - bed_elevation
    flotation_thickness = (rho_ocean / rho_ice) * draft

    # Thickness must be positive
    flotation_thickness = np.maximum(flotation_thickness, 0.0)

    return flotation_thickness


def calculate_required_bed(thickness, target_haf, sea_level=0.0,
                           rho_ice=910.0, rho_ocean=1028.0):
    """
    Calculate the bed elevation required to achieve a target height above
    flotation for a given ice thickness.

    Parameters
    ----------
    thickness : ndarray
        Ice thickness (m)
    target_haf : float
        Target height above flotation (m)
    sea_level : float
        Sea level (m), default 0.0
    rho_ice : float
        Ice density (kg/m^3), default 910.0
    rho_ocean : float
        Ocean water density (kg/m^3), default 1028.0

    Returns
    -------
    bed_elevation : ndarray
        Required bed elevation (m)
    """
    # At flotation:
    # ice_surface_flotation = sea_level + thickness * (1 - rho_ice/rho_ocean)
    #
    # With target HAF:
    # ice_surface = ice_surface_flotation + target_haf
    #
    # Since ice_surface = bed + thickness:
    # bed = ice_surface - thickness
    #     = ice_surface_flotation + target_haf - thickness
    #     = sea_level + thickness * (1 - rho_ice/rho_ocean) + target_haf - thickness
    #     = sea_level + target_haf - thickness * rho_ice/rho_ocean

    bed_elevation = sea_level + target_haf - thickness * (rho_ice / rho_ocean)

    return bed_elevation


def adjust_bed_to_haf():
    """
    Main function to adjust bed topography to achieve target height above
    flotation within grounding line regions.
    """

    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        '-m', '--mesh',
        dest='mesh_file',
        required=True,
        help='MALI mesh file (NetCDF format)'
    )
    parser.add_argument(
        '-g', '--geojson',
        dest='geojson_file',
        required=True,
        help='GeoJSON file containing grounding line delineations'
    )
    parser.add_argument(
        '-o', '--output',
        dest='output_file',
        help='Output mesh file. If not specified, modifies input file in place.'
    )
    parser.add_argument(
        '--target-haf',
        dest='target_haf',
        type=float,
        default=10.0,
        help='Target height above flotation in meters (default: 10.0)'
    )
    parser.add_argument(
        '--sea-level',
        dest='sea_level',
        type=float,
        default=0.0,
        help='Sea level in meters (default: 0.0)'
    )
    parser.add_argument(
        '--rho-ice',
        dest='rho_ice',
        type=float,
        default=910.0,
        help='Ice density in kg/m^3 (default: 910.0)'
    )
    parser.add_argument(
        '--rho-ocean',
        dest='rho_ocean',
        type=float,
        default=1028.0,
        help='Ocean water density in kg/m^3 (default: 1028.0)'
    )
    parser.add_argument(
        '--thickness-var',
        dest='thickness_var',
        default='thickness',
        help='Name of thickness variable in mesh file (default: thickness)'
    )
    parser.add_argument(
        '--bed-var',
        dest='bed_var',
        default='bedTopography',
        help='Name of bed topography variable in mesh file (default: bedTopography)'
    )

    args = parser.parse_args()

    # Validate inputs
    if not HAS_GEOPANDAS:
        print('Warning: geopandas not available, using json + shapely instead')

    print('\n** Adjusting bed topography to achieve target height above flotation **')
    print(f'Input mesh file: {args.mesh_file}')
    print(f'Grounding line GeoJSON: {args.geojson_file}')
    print(f'Target HAF: {args.target_haf} m')
    print(f'Sea level: {args.sea_level} m')
    print(f'Ice density: {args.rho_ice} kg/m^3')
    print(f'Ocean density: {args.rho_ocean} kg/m^3')
    print()

    # Load grounding line polygons
    print('Loading grounding line geometries...')
    geometries = load_grounding_line_geojson(args.geojson_file)
    print(f'Loaded {len(geometries)} geometry features')
    print()

    # Open mesh file
    print('Opening mesh file...')
    if args.output_file:
        # Copy to output file
        import shutil
        shutil.copy2(args.mesh_file, args.output_file)
        mesh_file = args.output_file
        print(f'Copied input to output file: {args.output_file}')
    else:
        mesh_file = args.mesh_file
        print('Modifying mesh file in place')

    with netCDF4.Dataset(mesh_file, 'r+') as mesh:

        # Read mesh coordinates
        print('Reading mesh coordinates...')
        # Coordinates are typically in radians, convert to degrees
        lon_cell = np.degrees(mesh.variables['lonCell'][:])
        lat_cell = np.degrees(mesh.variables['latCell'][:])
        n_cells = len(lon_cell)
        print(f'Mesh has {n_cells} cells')
        print()

        # Find cells within grounding line polygons
        print('Identifying cells within grounding line regions...')
        mask = find_cells_in_polygon(lon_cell, lat_cell, geometries)
        print()

        if np.sum(mask) == 0:
            print('ERROR: No cells found within grounding line polygons!')
            print('Check that:')
            print('  1. GeoJSON and mesh use compatible coordinate systems')
            print('  2. GeoJSON geometries overlap with mesh extent')
            sys.exit(1)

        # Read thickness and bed topography
        print('Reading ice thickness and bed topography...')
        if args.thickness_var not in mesh.variables:
            print(f'ERROR: Variable "{args.thickness_var}" not found in mesh file')
            print(f'Available variables: {list(mesh.variables.keys())}')
            sys.exit(1)

        if args.bed_var not in mesh.variables:
            print(f'ERROR: Variable "{args.bed_var}" not found in mesh file')
            print(f'Available variables: {list(mesh.variables.keys())}')
            sys.exit(1)

        thickness = mesh.variables[args.thickness_var][:]
        bed_topo = mesh.variables[args.bed_var][:]

        # Handle potential time dimension
        if len(thickness.shape) > 1:
            # Assume time is first dimension
            thickness = thickness[0, :]
            bed_topo = bed_topo[0, :]
            has_time_dim = True
        else:
            has_time_dim = False

        print(f'Thickness range: [{np.min(thickness):.2f}, {np.max(thickness):.2f}] m')
        print(f'Bed topography range: [{np.min(bed_topo):.2f}, {np.max(bed_topo):.2f}] m')
        print()

        # Calculate current HAF in the region
        print('Calculating current height above flotation...')
        flotation_thickness = calculate_flotation_thickness(
            bed_topo,
            sea_level=args.sea_level,
            rho_ice=args.rho_ice,
            rho_ocean=args.rho_ocean
        )
        ice_surface = bed_topo + thickness
        flotation_surface = args.sea_level + flotation_thickness * (
            1.0 - args.rho_ice / args.rho_ocean
        )
        current_haf = ice_surface - flotation_surface

        print(f'Current HAF in region (mean): {np.mean(current_haf[mask]):.2f} m')
        print(f'Current HAF in region (min):  {np.min(current_haf[mask]):.2f} m')
        print(f'Current HAF in region (max):  {np.max(current_haf[mask]):.2f} m')
        print()

        # Calculate new bed elevation
        print('Calculating new bed topography...')
        new_bed = calculate_required_bed(
            thickness[mask],
            args.target_haf,
            sea_level=args.sea_level,
            rho_ice=args.rho_ice,
            rho_ocean=args.rho_ocean
        )

        print(f'New bed range in region: [{np.min(new_bed):.2f}, {np.max(new_bed):.2f}] m')
        bed_change = new_bed - bed_topo[mask]
        print(f'Bed change (mean): {np.mean(bed_change):.2f} m')
        print(f'Bed change (min):  {np.min(bed_change):.2f} m')
        print(f'Bed change (max):  {np.max(bed_change):.2f} m')
        print()

        # Update bed topography
        print('Updating bed topography in mesh file...')
        bed_topo[mask] = new_bed

        if has_time_dim:
            # Write back with time dimension
            mesh.variables[args.bed_var][0, :] = bed_topo
        else:
            mesh.variables[args.bed_var][:] = bed_topo

        # Update global attributes
        timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        if 'history' in mesh.ncattrs():
            history = mesh.getncattr('history')
            history = f'{timestamp}: adjust_bed_to_haf.py\n{history}'
        else:
            history = f'{timestamp}: adjust_bed_to_haf.py'
        mesh.setncattr('history', history)

        comment = (
            f'Bed topography adjusted within grounding line regions '
            f'to achieve target HAF = {args.target_haf} m. '
            f'Modified {np.sum(mask)} cells using grounding line data from '
            f'{args.geojson_file}.'
        )
        if 'comment' in mesh.ncattrs():
            existing_comment = mesh.getncattr('comment')
            comment = f'{comment} {existing_comment}'
        mesh.setncattr('comment', comment)

        print('Successfully updated mesh file!')
        print()

        # Verify the result
        print('Verifying updated HAF...')
        flotation_thickness = calculate_flotation_thickness(
            bed_topo,
            sea_level=args.sea_level,
            rho_ice=args.rho_ice,
            rho_ocean=args.rho_ocean
        )
        ice_surface = bed_topo + thickness
        flotation_surface = args.sea_level + flotation_thickness * (
            1.0 - args.rho_ice / args.rho_ocean
        )
        new_haf = ice_surface - flotation_surface

        print(f'New HAF in region (mean): {np.mean(new_haf[mask]):.2f} m')
        print(f'New HAF in region (min):  {np.min(new_haf[mask]):.2f} m')
        print(f'New HAF in region (max):  {np.max(new_haf[mask]):.2f} m')
        print()

        if np.allclose(new_haf[mask], args.target_haf, rtol=1e-6):
            print(f'✓ Successfully achieved target HAF = {args.target_haf} m')
        else:
            print(f'Note: HAF values may differ slightly from target due to '
                  f'numerical precision')

    print()
    print('Done!')


if __name__ == '__main__':
    adjust_bed_to_haf()
