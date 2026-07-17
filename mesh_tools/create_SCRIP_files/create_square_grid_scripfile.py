#!/usr/bin/env python3
"""
Create a SCRIP-format file for a square regular lat-lon grid.
This grid is used for FastIsostasy coupling in MALI.

FastIsostasy requires a PERFECT SQUARE grid: nx = ny.

Usage:
    python create_square_grid_scripfile.py \\
        --mali_scripfile <path_to_mali.scripfile.nc> \\
        --nx 512 \\
        --output fastisostasy_grid.scripfile.nc

The script will:
1. Read the MALI scripfile to determine domain extent
2. Create a square regular grid (nx x nx) covering that domain
3. Write a SCRIP-format NetCDF file for use with ESMF_RegridWeightGen
"""

import sys
import argparse
import numpy as np

try:
    import netCDF4 as nc
except ImportError:
    print("ERROR: netCDF4 module not found. Install with: pip install netCDF4")
    sys.exit(1)


def read_mali_domain(mali_scripfile):
    """Read MALI scripfile and return domain bounds in degrees."""
    print(f"Reading MALI scripfile: {mali_scripfile}")
    ds = nc.Dataset(mali_scripfile, 'r')

    lat_rad = ds.variables['grid_center_lat'][:]
    lon_rad = ds.variables['grid_center_lon'][:]

    lat_deg = np.degrees(lat_rad)
    lon_deg = np.degrees(lon_rad)

    # Handle longitude wrapping
    lon_deg = np.where(lon_deg > 180, lon_deg - 360, lon_deg)

    bounds = {
        'lat_min': float(lat_deg.min()),
        'lat_max': float(lat_deg.max()),
        'lon_min': float(lon_deg.min()),
        'lon_max': float(lon_deg.max())
    }

    ds.close()

    print(f"  Domain: lat [{bounds['lat_min']:.2f}, {bounds['lat_max']:.2f}]°")
    print(f"          lon [{bounds['lon_min']:.2f}, {bounds['lon_max']:.2f}]°")

    return bounds


def create_square_grid_scripfile(bounds, nx, output_file, padding_km=0):
    """
    Create a SCRIP-format file for a square regular grid.

    Parameters:
    -----------
    bounds : dict
        Domain bounds with keys: lat_min, lat_max, lon_min, lon_max (degrees)
    nx : int
        Number of grid points in each dimension (must be square: nx x nx)
    output_file : str
        Output SCRIP NetCDF filename
    padding_km : float
        Padding around domain in km (converted to degrees approximately)
    """
    ny = nx  # Square grid requirement

    # Apply padding (rough conversion: 1° ≈ 111 km at mid-latitudes)
    if padding_km > 0:
        padding_deg = padding_km / 111.0
        bounds['lat_min'] -= padding_deg
        bounds['lat_max'] += padding_deg
        bounds['lon_min'] -= padding_deg
        bounds['lon_max'] += padding_deg
        print(f"  Applying {padding_km} km padding (~{padding_deg:.2f}°)")

    # Create 1D coordinate arrays
    lats = np.linspace(bounds['lat_min'], bounds['lat_max'], ny)
    lons = np.linspace(bounds['lon_min'], bounds['lon_max'], nx)

    dlat = lats[1] - lats[0]
    dlon = lons[1] - lons[0]

    print(f"Creating {nx}x{ny} square grid:")
    print(f"  Resolution: dlat={dlat:.4f}°, dlon={dlon:.4f}°")
    print(f"  Total cells: {nx*ny}")

    # Create 2D meshgrid
    lon_2d, lat_2d = np.meshgrid(lons, lats)

    # Flatten to 1D (SCRIP format)
    grid_size = nx * ny
    grid_center_lat = lat_2d.flatten()
    grid_center_lon = lon_2d.flatten()

    # Create corner arrays (4 corners for regular grid)
    grid_corner_lat = np.zeros((grid_size, 4))
    grid_corner_lon = np.zeros((grid_size, 4))

    # Corner offsets (counter-clockwise from bottom-left)
    corner_dlat = np.array([-dlat/2, -dlat/2, dlat/2, dlat/2])
    corner_dlon = np.array([-dlon/2, dlon/2, dlon/2, -dlon/2])

    for i in range(grid_size):
        grid_corner_lat[i, :] = grid_center_lat[i] + corner_dlat
        grid_corner_lon[i, :] = grid_center_lon[i] + corner_dlon

    # Convert to radians
    grid_center_lat_rad = np.radians(grid_center_lat)
    grid_center_lon_rad = np.radians(grid_center_lon)
    grid_corner_lat_rad = np.radians(grid_corner_lat)
    grid_corner_lon_rad = np.radians(grid_corner_lon)

    # Calculate grid cell areas (approximate spherical)
    R_earth = 6.371e6  # Earth radius in meters
    grid_area = np.zeros(grid_size)
    for i in range(grid_size):
        lat1 = grid_corner_lat_rad[i, 0]
        lat2 = grid_corner_lat_rad[i, 2]
        lon1 = grid_corner_lon_rad[i, 0]
        lon2 = grid_corner_lon_rad[i, 1]
        # Area = R^2 * |sin(lat2) - sin(lat1)| * |lon2 - lon1|
        grid_area[i] = R_earth**2 * abs(np.sin(lat2) - np.sin(lat1)) * abs(lon2 - lon1)

    # Normalize to steradians
    grid_area_rad = grid_area / R_earth**2

    # All cells are active (imask = 1)
    grid_imask = np.ones(grid_size, dtype=np.int32)

    # Write SCRIP file
    print(f"Writing SCRIP file: {output_file}")
    ds_out = nc.Dataset(output_file, 'w', format='NETCDF3_64BIT_OFFSET')

    # Dimensions
    ds_out.createDimension('grid_size', grid_size)
    ds_out.createDimension('grid_corners', 4)
    ds_out.createDimension('grid_rank', 2)

    # Variables
    var_lat = ds_out.createVariable('grid_center_lat', 'f8', ('grid_size',))
    var_lat.units = 'radians'
    var_lat[:] = grid_center_lat_rad

    var_lon = ds_out.createVariable('grid_center_lon', 'f8', ('grid_size',))
    var_lon.units = 'radians'
    var_lon[:] = grid_center_lon_rad

    var_clat = ds_out.createVariable('grid_corner_lat', 'f8', ('grid_size', 'grid_corners'))
    var_clat.units = 'radians'
    var_clat[:, :] = grid_corner_lat_rad

    var_clon = ds_out.createVariable('grid_corner_lon', 'f8', ('grid_size', 'grid_corners'))
    var_clon.units = 'radians'
    var_clon[:, :] = grid_corner_lon_rad

    var_area = ds_out.createVariable('grid_area', 'f8', ('grid_size',))
    var_area.units = 'radian^2'
    var_area[:] = grid_area_rad

    var_mask = ds_out.createVariable('grid_imask', 'i4', ('grid_size',))
    var_mask.units = 'unitless'
    var_mask[:] = grid_imask

    var_dims = ds_out.createVariable('grid_dims', 'i4', ('grid_rank',))
    var_dims[:] = [nx, ny]

    # Global attributes
    ds_out.title = f'FastIsostasy square grid {nx}x{ny}'
    ds_out.history = 'Created by create_square_grid_scripfile.py'

    ds_out.close()
    print(f"SUCCESS: Created {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description='Create SCRIP file for FastIsostasy square grid',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument('--mali_scripfile', required=True,
                        help='Path to MALI mesh SCRIP file')
    parser.add_argument('--nx', type=int, required=True,
                        help='Grid dimension (creates nx x nx square grid)')
    parser.add_argument('--output', required=True,
                        help='Output SCRIP filename')
    parser.add_argument('--padding', type=float, default=0,
                        help='Padding around domain in km (default: 0)')

    args = parser.parse_args()

    # Validate square grid
    if args.nx <= 0:
        print(f"ERROR: nx must be positive, got {args.nx}")
        sys.exit(1)

    # Read MALI domain
    bounds = read_mali_domain(args.mali_scripfile)

    # Create square grid scripfile
    create_square_grid_scripfile(bounds, args.nx, args.output, args.padding)

    print("\nNext steps:")
    print(f"1. Generate weight files with ESMF_RegridWeightGen:")
    print(f"   # MALI to FastIsostasy")
    print(f"   ESMF_RegridWeightGen -s {args.mali_scripfile} -d {args.output} \\")
    print(f"       -w mpas_to_fastisostasy.nc -i --64bit_offset --src_regional")
    print(f"   # FastIsostasy to MALI")
    print(f"   ESMF_RegridWeightGen -s {args.output} -d {args.mali_scripfile} \\")
    print(f"       -w fastisostasy_to_mpas.nc -i --64bit_offset --dst_regional")


if __name__ == '__main__':
    main()
