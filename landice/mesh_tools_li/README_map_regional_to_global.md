# map_regional_to_global_mesh.py

Map cells from a regional mesh back to a global mesh based on exact coordinate matching.

## Purpose

This script is designed for workflows where:
1. A regional mesh was extracted from a larger global mesh (e.g., using compass subdomain extractor)
2. Variables in the regional mesh were modified
3. The modified values need to be copied back to the corresponding cells in the global mesh
4. All other global mesh cells must remain unchanged

The script uses **exact coordinate matching** (not interpolation) since the regional mesh cells are a subset of the global mesh cells with identical coordinates.

## Quick Start

```bash
python map_regional_to_global_mesh.py \
    --regional regional_mesh.nc \
    --global global_mesh.nc \
    --output updated_global_mesh.nc \
    --vars bedTopography thickness
```

## Usage

```
python map_regional_to_global_mesh.py [-h] -r FILENAME -g FILENAME -o FILENAME
                                      -v VARIABLES [VARIABLES ...]
                                      [--coord-type {spherical,planar}]
                                      [--tolerance TOLERANCE]
                                      [--verify-only]
```

### Required Arguments

- `-r`, `--regional FILENAME`: Regional mesh file (NetCDF format)
- `-g`, `--global FILENAME`: Global mesh file (NetCDF format)
- `-o`, `--output FILENAME`: Output file (must differ from global file)
- `-v`, `--vars VARIABLES`: Variable(s) to copy from regional to global mesh

### Optional Arguments

- `--coord-type {spherical,planar}`: Coordinate type (default: spherical)
  - `spherical`: Use lonCell/latCell (for spherical meshes)
  - `planar`: Use xCell/yCell (for planar meshes)
- `--tolerance TOLERANCE`: Coordinate matching tolerance (default: 1e-10)
- `--verify-only`: Only verify mapping without copying data

## Example Workflow

### 1. Extract regional mesh from global mesh
```bash
# Using compass subdomain extractor or similar tool
compass extract-subdomain \
    --mesh global_mesh.nc \
    --output regional_mesh.nc \
    --region amundsen_sea
```

### 2. Modify regional mesh
```bash
# Make modifications to the regional mesh
python adjust_bed_to_haf.py \
    --mesh regional_mesh.nc \
    --geojson grounding_line.geojson \
    --target-haf 10.0
```

### 3. Map back to global mesh
```bash
# Copy modified variables back to global mesh
python map_regional_to_global_mesh.py \
    --regional regional_mesh.nc \
    --global global_mesh.nc \
    --output global_mesh_updated.nc \
    --vars bedTopography
```

## Use Cases

### Copy modified bed topography
```bash
python map_regional_to_global_mesh.py \
    -r ASE_mesh_modified.nc \
    -g AIS_mesh.nc \
    -o AIS_mesh_updated.nc \
    -v bedTopography
```

### Copy multiple variables
```bash
python map_regional_to_global_mesh.py \
    -r ASE_mesh_modified.nc \
    -g AIS_mesh.nc \
    -o AIS_mesh_updated.nc \
    -v bedTopography thickness surfaceSpeed temperature
```

### Verify mapping without copying
```bash
python map_regional_to_global_mesh.py \
    -r ASE_mesh.nc \
    -g AIS_mesh.nc \
    -o dummy_output.nc \
    -v bedTopography \
    --verify-only
```

### Use planar coordinates
```bash
python map_regional_to_global_mesh.py \
    -r regional_planar.nc \
    -g global_planar.nc \
    -o global_updated.nc \
    -v thickness \
    --coord-type planar
```

## How It Works

1. **Load both meshes**: Opens regional and global mesh files
2. **Extract coordinates**: Gets lonCell/latCell (spherical) or xCell/yCell (planar)
3. **Build mapping**: For each regional cell, finds the global cell with matching coordinates
4. **Verify mapping**: Checks that all matches are within tolerance
5. **Copy variables**: Updates specified variables in global mesh for matched cells only
6. **Write output**: Saves updated global mesh to new file

## Coordinate Matching

The script uses exact coordinate matching with a small tolerance (default 1e-10) to account for floating-point precision:

- **Spherical**: Matches based on (lonCell, latCell) in radians
- **Planar**: Matches based on (xCell, yCell) in meters

For each regional cell, the script finds the global cell with minimum distance. If the distance exceeds the tolerance, the regional cell is marked as unmatched.

## Variable Handling

The script automatically handles variables with different dimensions:

- **1D** `(nCells)`: Direct cell-to-cell copy
- **2D** `(nCells, nVertLevels)`: Copies all vertical levels
- **2D** `(Time, nCells)`: Copies all time levels
- **3D** `(Time, nCells, nVertLevels)`: Copies all time and vertical levels

## Error Checking

The script includes several safety checks:

- Verifies output file differs from global file (prevents accidental overwrite)
- Checks that all requested variables exist in both meshes
- Reports unmatched regional cells
- Verifies coordinate matching is within tolerance
- Updates file metadata with history and comments

## Output

The script reports:
- Number of cells in each mesh
- Number of matched cells
- Number of unmatched regional cells (if any)
- Coordinate matching errors (mean and max)
- Variables successfully copied

## Notes

- The output file must differ from the global input file (safety feature)
- Unmatched regional cells are skipped (with warning)
- The global mesh is loaded into memory, so ensure sufficient RAM for large meshes
- Original global file is never modified (read-only)
- All cells in global mesh not matching regional cells remain unchanged

## Dependencies

Required Python packages:
- xarray
- numpy

## Common Issues

### "No matching cells found"
- Check that regional mesh was actually extracted from this global mesh
- Verify coordinate type (spherical vs planar)
- Try increasing tolerance

### "WARNING: Some matches exceed tolerance"
- Regional and global meshes may have been modified with different precision
- Try increasing tolerance with `--tolerance 1e-6`

### Memory issues
- For very large global meshes, ensure sufficient RAM
- Consider processing variables one at a time

## See Also

- `adjust_bed_to_haf.py`: Tool for modifying bed topography
- `compass extract-subdomain`: Tool for extracting regional meshes
