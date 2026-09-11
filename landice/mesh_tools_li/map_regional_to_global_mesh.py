#!/usr/bin/env python
'''
Map cells from a regional mesh back to a global mesh based on exact coordinate
matching.

This script is designed for workflows where a regional mesh has been extracted
from a larger global mesh (e.g., using compass subdomain extractor), modified,
and then needs to be mapped back to the global mesh. Since the regional mesh
cells are a subset of the global mesh cells with identical coordinates, this
script creates an exact coordinate-based mapping rather than using interpolation.

The script:
1. Loads both regional and global meshes
2. Finds cells in the global mesh that match regional mesh coordinates exactly
3. Copies specified variables from regional to global mesh for matched cells
4. Leaves all other global mesh cells unchanged

Trevor Hillebrand, 2026
'''

import sys
import numpy as np
import xarray as xr
from argparse import ArgumentParser
from datetime import datetime


def parse_args():
    parser = ArgumentParser(description=__doc__,
                           formatter_class=lambda prog: ArgumentParser.
                           RawDescriptionHelpFormatter(prog, max_help_position=30))
    parser.add_argument('-r', '--regional', dest='regional_file', required=True,
                       metavar='FILENAME',
                       help='Regional mesh file (NetCDF format)')
    parser.add_argument('-g', '--global', dest='global_file', required=True,
                       metavar='FILENAME',
                       help='Global mesh file (NetCDF format)')
    parser.add_argument('-o', '--output', dest='output_file', required=True,
                       metavar='FILENAME',
                       help='Output mesh file (must differ from global file)')
    parser.add_argument('-v', '--vars', dest='variables', required=True,
                       nargs='+',
                       help='Variable(s) to copy from regional to global mesh')
    parser.add_argument('--coord-type', dest='coord_type',
                       choices=['spherical', 'planar'], default='spherical',
                       help='Coordinate type: spherical (lon/lat) or planar '
                            '(x/y) (default: spherical)')
    parser.add_argument('--tolerance', dest='tolerance', type=float,
                       default=1e-10,
                       help='Coordinate matching tolerance (default: 1e-10)')
    parser.add_argument('--verify-only', dest='verify_only', action='store_true',
                       help='Only verify mapping without copying data')

    return parser.parse_args()


def build_coordinate_mapping(regional_coords, global_coords, tolerance=1e-10):
    '''
    Build a mapping from regional mesh cell indices to global mesh cell indices
    based on exact coordinate matching.

    Parameters
    ----------
    regional_coords : ndarray, shape (n_regional, 2)
        Regional mesh coordinates (x, y) or (lon, lat)
    global_coords : ndarray, shape (n_global, 2)
        Global mesh coordinates (x, y) or (lon, lat)
    tolerance : float
        Maximum distance for considering coordinates as matching

    Returns
    -------
    mapping : dict
        Dictionary mapping regional cell index to global cell index
    unmatched_regional : list
        List of regional cell indices that have no match in global mesh
    '''
    n_regional = regional_coords.shape[0]
    n_global = global_coords.shape[0]

    print(f'Building coordinate mapping...')
    print(f'  Regional mesh: {n_regional} cells')
    print(f'  Global mesh: {n_global} cells')
    print(f'  Tolerance: {tolerance}')

    mapping = {}
    unmatched_regional = []

    # For each regional cell, find matching global cell
    for i in range(n_regional):
        if i % 1000 == 0:
            print(f'  Processed {i}/{n_regional} regional cells...')

        regional_coord = regional_coords[i, :]

        # Calculate distances to all global cells
        distances = np.sqrt(np.sum((global_coords - regional_coord)**2, axis=1))

        # Find closest match
        min_idx = np.argmin(distances)
        min_dist = distances[min_idx]

        if min_dist <= tolerance:
            mapping[i] = min_idx
        else:
            unmatched_regional.append(i)

    print(f'  Completed mapping!')
    print(f'  Matched cells: {len(mapping)} / {n_regional}')
    if unmatched_regional:
        print(f'  WARNING: {len(unmatched_regional)} regional cells have no '
              f'match in global mesh')

    return mapping, unmatched_regional


def verify_mapping(mapping, regional_coords, global_coords, tolerance):
    '''
    Verify that the mapping is correct by checking coordinate differences.

    Parameters
    ----------
    mapping : dict
        Regional to global cell index mapping
    regional_coords : ndarray
        Regional mesh coordinates
    global_coords : ndarray
        Global mesh coordinates
    tolerance : float
        Matching tolerance

    Returns
    -------
    max_error : float
        Maximum coordinate error
    '''
    print('Verifying mapping...')

    errors = []
    for regional_idx, global_idx in mapping.items():
        regional_coord = regional_coords[regional_idx, :]
        global_coord = global_coords[global_idx, :]
        error = np.sqrt(np.sum((regional_coord - global_coord)**2))
        errors.append(error)

    errors = np.array(errors)
    max_error = np.max(errors)
    mean_error = np.mean(errors)

    print(f'  Mean coordinate error: {mean_error:.2e}')
    print(f'  Max coordinate error: {max_error:.2e}')
    print(f'  Tolerance: {tolerance:.2e}')

    if max_error > tolerance:
        print('  WARNING: Some matches exceed tolerance!')
        return False
    else:
        print('  ✓ All matches within tolerance')
        return True


def copy_variables(regional_ds, global_ds, mapping, variables):
    '''
    Copy specified variables from regional to global mesh using the mapping.

    Parameters
    ----------
    regional_ds : xarray.Dataset
        Regional mesh dataset
    global_ds : xarray.Dataset
        Global mesh dataset
    mapping : dict
        Regional to global cell index mapping
    variables : list
        List of variable names to copy

    Returns
    -------
    global_ds : xarray.Dataset
        Modified global mesh dataset
    '''
    print(f'Copying variables: {", ".join(variables)}')

    for var in variables:
        print(f'  Processing {var}...')

        if var not in regional_ds.variables:
            print(f'    ERROR: Variable "{var}" not found in regional mesh')
            continue

        if var not in global_ds.variables:
            print(f'    ERROR: Variable "{var}" not found in global mesh')
            continue

        regional_var = regional_ds[var].values
        global_var = global_ds[var].values

        # Get variable shape
        var_shape = regional_var.shape
        print(f'    Regional shape: {var_shape}')
        print(f'    Global shape: {global_var.shape}')

        # Handle different dimensionalities
        if len(var_shape) == 1:
            # 1D variable (nCells)
            for regional_idx, global_idx in mapping.items():
                global_var[global_idx] = regional_var[regional_idx]

        elif len(var_shape) == 2:
            # 2D variable, could be (Time, nCells) or (nCells, nVertLevels)
            dim_names = regional_ds[var].dims

            if 'nCells' in dim_names:
                cell_dim = dim_names.index('nCells')

                if cell_dim == 0:
                    # (nCells, nVertLevels)
                    for regional_idx, global_idx in mapping.items():
                        global_var[global_idx, :] = regional_var[regional_idx, :]
                else:
                    # (Time, nCells) or (nVertLevels, nCells)
                    for regional_idx, global_idx in mapping.items():
                        global_var[:, global_idx] = regional_var[:, regional_idx]
            else:
                print(f'    WARNING: Cannot determine cell dimension for {var}')
                continue

        elif len(var_shape) == 3:
            # 3D variable (Time, nCells, nVertLevels)
            dim_names = regional_ds[var].dims

            if 'nCells' in dim_names:
                cell_dim = dim_names.index('nCells')

                if cell_dim == 1:
                    # (Time, nCells, nVertLevels)
                    for regional_idx, global_idx in mapping.items():
                        global_var[:, global_idx, :] = \
                            regional_var[:, regional_idx, :]
                else:
                    print(f'    WARNING: Unexpected cell dimension position '
                          f'for {var}')
                    continue
            else:
                print(f'    WARNING: Cannot determine cell dimension for {var}')
                continue

        else:
            print(f'    WARNING: Variable {var} has {len(var_shape)} dimensions,'
                  f' not currently supported')
            continue

        # Update the global dataset
        global_ds[var].values = global_var
        print(f'    ✓ Copied {var}')

    return global_ds


def main():
    '''
    Main function to map regional mesh data to global mesh.
    '''

    args = parse_args()

    # Validate that output file differs from global file
    if args.output_file == args.global_file:
        print('ERROR: Output file must differ from global file')
        print('       (to avoid overwriting the original global mesh)')
        sys.exit(1)

    print('\n** Mapping regional mesh to global mesh **')
    print(f'Regional mesh: {args.regional_file}')
    print(f'Global mesh: {args.global_file}')
    print(f'Output file: {args.output_file}')
    print(f'Variables to copy: {", ".join(args.variables)}')
    print(f'Coordinate type: {args.coord_type}')
    print()

    # Load meshes
    print('Loading regional mesh...')
    regional_ds = xr.open_dataset(args.regional_file)
    n_regional = regional_ds.dims['nCells']
    print(f'  Regional mesh has {n_regional} cells')

    print('Loading global mesh...')
    global_ds = xr.open_dataset(args.global_file)
    n_global = global_ds.dims['nCells']
    print(f'  Global mesh has {n_global} cells')
    print()

    # Get coordinates based on type
    if args.coord_type == 'spherical':
        print('Using spherical coordinates (lonCell, latCell)...')
        regional_coords = np.column_stack([
            regional_ds['lonCell'].values,
            regional_ds['latCell'].values
        ])
        global_coords = np.column_stack([
            global_ds['lonCell'].values,
            global_ds['latCell'].values
        ])
    else:
        print('Using planar coordinates (xCell, yCell)...')
        regional_coords = np.column_stack([
            regional_ds['xCell'].values,
            regional_ds['yCell'].values
        ])
        global_coords = np.column_stack([
            global_ds['xCell'].values,
            global_ds['yCell'].values
        ])

    print()

    # Build coordinate mapping
    mapping, unmatched = build_coordinate_mapping(
        regional_coords, global_coords, args.tolerance
    )

    if not mapping:
        print('ERROR: No matching cells found!')
        sys.exit(1)

    print()

    # Verify mapping
    mapping_ok = verify_mapping(mapping, regional_coords, global_coords,
                                args.tolerance)
    if not mapping_ok:
        print('WARNING: Mapping verification failed, but continuing...')

    print()

    # Report unmatched cells if any
    if unmatched:
        print(f'WARNING: {len(unmatched)} regional cells have no match in '
              f'global mesh')
        print('These cells will be skipped.')
        print()

    # If verify-only mode, stop here
    if args.verify_only:
        print('Verify-only mode: stopping before copying data')
        return

    # Copy variables
    print('Copying variables from regional to global mesh...')
    global_ds = copy_variables(regional_ds, global_ds, mapping, args.variables)
    print()

    # Update global attributes
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    if 'history' in global_ds.attrs:
        history = global_ds.attrs['history']
        history = f'{timestamp}: map_regional_to_global_mesh.py\n{history}'
    else:
        history = f'{timestamp}: map_regional_to_global_mesh.py'
    global_ds.attrs['history'] = history

    comment = (
        f'Variables {", ".join(args.variables)} updated from regional mesh '
        f'{args.regional_file}. {len(mapping)} cells modified.'
    )
    if 'comment' in global_ds.attrs:
        existing_comment = global_ds.attrs['comment']
        comment = f'{comment} {existing_comment}'
    global_ds.attrs['comment'] = comment

    # Write output
    print(f'Writing output to {args.output_file}...')
    # Use encoding to preserve data types and compression
    encoding = {var: {'_FillValue': None} for var in global_ds.data_vars}
    global_ds.to_netcdf(args.output_file, encoding=encoding)

    print('Done!')
    print()
    print(f'Summary:')
    print(f'  Regional cells: {n_regional}')
    print(f'  Global cells: {n_global}')
    print(f'  Matched cells: {len(mapping)}')
    print(f'  Unmatched regional cells: {len(unmatched)}')
    print(f'  Variables copied: {len(args.variables)}')
    print(f'  Output: {args.output_file}')

    # Clean up
    regional_ds.close()
    global_ds.close()


if __name__ == '__main__':
    main()
