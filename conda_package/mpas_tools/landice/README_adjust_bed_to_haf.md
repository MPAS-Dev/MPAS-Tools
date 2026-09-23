# Adjust Bed to Height Above Flotation

This tool adjusts bed topography within grounding line regions to achieve a specified height above flotation (HAF).

## Quick Start

```bash
adjust_bed_to_haf \
    --mesh input_mesh.nc \
    --geojson grounding_line.geojson \
    --output output_mesh.nc \
    --target-haf 10.0
```

## What is Height Above Flotation?

Height above flotation (HAF) measures how far ice is from hydrostatic equilibrium with ocean water:

- **HAF = 0**: Ice is exactly at flotation (neither grounded nor floating freely)
- **HAF > 0**: Ice surface is above the flotation level (typical for grounded ice)
- **HAF < 0**: Ice surface is below flotation level (should not occur in steady state)

## How It Works

The tool:

1. Loads grounding line polygons from a GeoJSON file
2. Identifies mesh cells within these polygons
3. Calculates required bed elevation to achieve target HAF
4. Updates the `bedTopography` field in the mesh file
5. Preserves ice thickness (only bed is modified)

## Physics

For a given ice thickness H and target HAF, the required bed elevation is:

```
bed = sea_level + HAF_target - H × (ρ_ice / ρ_ocean)
```

Where:
- ρ_ice = 910 kg/m³ (default)
- ρ_ocean = 1028 kg/m³ (default)

## Input File Requirements

### Mesh File (NetCDF)
Must contain:
- `lonCell`, `latCell`: Cell center coordinates (in radians)
- `thickness`: Ice thickness (m)
- `bedTopography`: Bed elevation (m, positive up)

### GeoJSON File
Must contain:
- Polygon or MultiLineString geometries
- Coordinates in WGS 84 (longitude, latitude in degrees)

## Example: Thwaites Glacier

```bash
adjust_bed_to_haf \
    --mesh /path/to/AIS_4to20km_mesh.nc \
    --geojson /path/to/Thwaites_GL_2014_pinning_points.geojson \
    --output AIS_4to20km_mesh_haf15m.nc \
    --target-haf 15.0 \
    --rho-ice 918.0
```

## Common Options

- `--target-haf 10.0`: Set target height above flotation (meters)
- `--sea-level 0.0`: Set sea level (meters)
- `--rho-ice 910.0`: Set ice density (kg/m³)
- `--rho-ocean 1028.0`: Set ocean water density (kg/m³)
- `--output file.nc`: Save to new file (otherwise modifies in place)

## Troubleshooting

### "No cells found within grounding line polygons"

This usually means:
1. Coordinate system mismatch (mesh in radians vs. GeoJSON in degrees)
2. Grounding line doesn't overlap with mesh domain
3. GeoJSON uses wrong projection

**Solution**: Verify that your GeoJSON is in WGS 84 (EPSG:4326) coordinates.

### Variable not found

The tool looks for `thickness` and `bedTopography` by default. If your mesh uses different names:

```bash
adjust_bed_to_haf \
    --thickness-var myThickness \
    --bed-var myBedTopo \
    ...
```

## Use Cases

### 1. Initialize grounding line regions
Set grounding line areas to a specific HAF for model initialization:

```bash
adjust_bed_to_haf --mesh init.nc --geojson gl.geojson --target-haf 10.0
```

### 2. Sensitivity experiments
Test model response to different HAF values:

```bash
for haf in 5 10 15 20; do
    adjust_bed_to_haf \
        --mesh base.nc \
        --geojson gl.geojson \
        --output mesh_haf${haf}m.nc \
        --target-haf $haf
done
```

### 3. Adjust pinning points
Modify bed topography at pinning points to control ice sheet stability:

```bash
adjust_bed_to_haf \
    --mesh mesh.nc \
    --geojson pinning_points.geojson \
    --target-haf 5.0 \
    --output mesh_weak_pinning.nc
```

## Dependencies

**Required:**
- netCDF4
- numpy
- shapely

**Optional:**
- geopandas (recommended for better performance)

## See Also

- [Full documentation](../../docs/landice/adjust_bed_to_haf.rst)
- [MPAS-Tools landice module](https://mpas-dev.github.io/MPAS-Tools/stable/landice.html)
