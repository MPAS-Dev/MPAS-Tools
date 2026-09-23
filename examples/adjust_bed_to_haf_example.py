#!/usr/bin/env python3
"""
Example script demonstrating how to use the adjust_bed_to_haf tool
programmatically (as a Python module) rather than from the command line.

This is useful if you want to integrate the functionality into a larger
workflow or perform additional processing.
"""

import sys
import numpy as np
import netCDF4
from shapely.geometry import Point, shape

# Import the functions from the landice module
# (assuming mpas_tools is installed)
try:
    from mpas_tools.landice.adjust_bed_to_haf import (
        load_grounding_line_geojson,
        find_cells_in_polygon,
        calculate_flotation_thickness,
        calculate_required_bed
    )
except ImportError:
    print("Error: mpas_tools not installed or not in PYTHONPATH")
    print("Install with: cd conda_package && pip install -e .")
    sys.exit(1)


def example_workflow():
    """
    Example workflow for adjusting bed topography to achieve target HAF.
    """

    # Define input files
    mesh_file = "/Users/trhille/Documents/ISMIP7/Antarctica/mesh/4to20km/AIS_4to20km_r03_20260910_ASE_extracted.nc"
    geojson_file = "/Users/trhille/Documents/ISMIP7/Antarctica/wild_grounding_lines/Thwaites_GLs_2014_201920/Thwaites_GL_2014_pinning_points.geojson"
    output_file = "output_mesh_programmatic.nc"

    # Define parameters
    target_haf = 12.0  # meters
    rho_ice = 910.0  # kg/m^3
    rho_ocean = 1028.0  # kg/m^3
    sea_level = 0.0  # meters

    print("="*60)
    print("Adjust Bed to HAF - Programmatic Example")
    print("="*60)
    print()

    # Step 1: Load grounding line geometries
    print("Step 1: Loading grounding line geometries...")
    geometries = load_grounding_line_geojson(geojson_file)
    print(f"  Loaded {len(geometries)} geometry features")
    print()

    # Step 2: Open mesh file and get coordinates
    print("Step 2: Opening mesh file...")
    with netCDF4.Dataset(mesh_file, 'r') as mesh:
        lon_cell = np.degrees(mesh.variables['lonCell'][:])
        lat_cell = np.degrees(mesh.variables['latCell'][:])
        thickness = mesh.variables['thickness'][:]
        bed_topo = mesh.variables['bedTopography'][:]

        # Handle potential time dimension
        if len(thickness.shape) > 1:
            thickness = thickness[0, :]
            bed_topo = bed_topo[0, :]

    print(f"  Mesh has {len(lon_cell)} cells")
    print(f"  Thickness range: [{np.min(thickness):.1f}, {np.max(thickness):.1f}] m")
    print(f"  Bed range: [{np.min(bed_topo):.1f}, {np.max(bed_topo):.1f}] m")
    print()

    # Step 3: Find cells within grounding line polygons
    print("Step 3: Finding cells within grounding line regions...")
    mask = find_cells_in_polygon(lon_cell, lat_cell, geometries)
    n_cells_in_region = np.sum(mask)
    print(f"  Found {n_cells_in_region} cells within grounding line")
    print(f"  ({100.0 * n_cells_in_region / len(mask):.2f}% of total mesh)")
    print()

    if n_cells_in_region == 0:
        print("ERROR: No cells found in region!")
        return

    # Step 4: Calculate current HAF
    print("Step 4: Calculating current height above flotation...")
    flotation_thickness = calculate_flotation_thickness(
        bed_topo,
        sea_level=sea_level,
        rho_ice=rho_ice,
        rho_ocean=rho_ocean
    )
    ice_surface = bed_topo + thickness
    flotation_surface = sea_level + flotation_thickness * (1.0 - rho_ice / rho_ocean)
    current_haf = ice_surface - flotation_surface

    print(f"  Current HAF in region:")
    print(f"    Mean: {np.mean(current_haf[mask]):.2f} m")
    print(f"    Min:  {np.min(current_haf[mask]):.2f} m")
    print(f"    Max:  {np.max(current_haf[mask]):.2f} m")
    print(f"    Std:  {np.std(current_haf[mask]):.2f} m")
    print()

    # Step 5: Calculate new bed elevation
    print(f"Step 5: Calculating new bed for target HAF = {target_haf} m...")
    new_bed = calculate_required_bed(
        thickness[mask],
        target_haf,
        sea_level=sea_level,
        rho_ice=rho_ice,
        rho_ocean=rho_ocean
    )

    bed_change = new_bed - bed_topo[mask]
    print(f"  Required bed change:")
    print(f"    Mean: {np.mean(bed_change):.2f} m")
    print(f"    Min:  {np.min(bed_change):.2f} m")
    print(f"    Max:  {np.max(bed_change):.2f} m")
    print(f"    Std:  {np.std(bed_change):.2f} m")
    print()

    # Step 6: Update mesh and save
    print("Step 6: Updating mesh file...")
    import shutil
    shutil.copy2(mesh_file, output_file)

    with netCDF4.Dataset(output_file, 'r+') as mesh:
        bed_var = mesh.variables['bedTopography']

        # Read current bed
        if len(bed_var.shape) > 1:
            current_bed = bed_var[0, :]
        else:
            current_bed = bed_var[:]

        # Update bed in region
        current_bed[mask] = new_bed

        # Write back
        if len(bed_var.shape) > 1:
            bed_var[0, :] = current_bed
        else:
            bed_var[:] = current_bed

        # Update metadata
        from datetime import datetime
        timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        history = f'{timestamp}: Programmatic adjust_bed_to_haf example\n'
        if 'history' in mesh.ncattrs():
            history += mesh.getncattr('history')
        mesh.setncattr('history', history)

    print(f"  Saved to: {output_file}")
    print()

    # Step 7: Verify result
    print("Step 7: Verifying result...")
    with netCDF4.Dataset(output_file, 'r') as mesh:
        new_bed_topo = mesh.variables['bedTopography'][:]
        if len(new_bed_topo.shape) > 1:
            new_bed_topo = new_bed_topo[0, :]

    new_flotation_thickness = calculate_flotation_thickness(
        new_bed_topo,
        sea_level=sea_level,
        rho_ice=rho_ice,
        rho_ocean=rho_ocean
    )
    new_ice_surface = new_bed_topo + thickness
    new_flotation_surface = sea_level + new_flotation_thickness * (1.0 - rho_ice / rho_ocean)
    verified_haf = new_ice_surface - new_flotation_surface

    print(f"  Verified HAF in region:")
    print(f"    Mean: {np.mean(verified_haf[mask]):.2f} m")
    print(f"    Min:  {np.min(verified_haf[mask]):.2f} m")
    print(f"    Max:  {np.max(verified_haf[mask]):.2f} m")
    print()

    if np.allclose(verified_haf[mask], target_haf, rtol=1e-6):
        print("✓ SUCCESS: Target HAF achieved!")
    else:
        diff = np.abs(verified_haf[mask] - target_haf)
        print(f"  Max difference from target: {np.max(diff):.6f} m")
        print("  (small differences expected due to numerical precision)")

    print()
    print("="*60)
    print("Example completed successfully!")
    print("="*60)


if __name__ == '__main__':
    example_workflow()
