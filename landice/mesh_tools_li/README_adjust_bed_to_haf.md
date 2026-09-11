# adjust_bed_to_haf.py

Adjust bed topography within grounding line regions to achieve a specified height above flotation (HAF).

## Quick Start

```bash
python adjust_bed_to_haf.py \
    --mesh input_mesh.nc \
    --geojson grounding_line.geojson \
    --output output_mesh.nc \
    --target-haf 10.0
```

## Usage

```
python adjust_bed_to_haf.py [-h] -m FILENAME -g FILENAME [-o FILENAME]
                            [--target-haf TARGET_HAF] [--sea-level SEA_LEVEL]
                            [--rho-ice RHO_ICE] [--rho-ocean RHO_OCEAN]
                            [--thickness-var THICKNESS_VAR] [--bed-var BED_VAR]
```

### Required Arguments

- `-m`, `--mesh FILENAME`: MALI mesh file (NetCDF format)
- `-g`, `--geojson FILENAME`: GeoJSON file containing grounding line delineations

### Optional Arguments

- `-o`, `--output FILENAME`: Output mesh file (default: modifies input in place)
- `--target-haf TARGET_HAF`: Target height above flotation in meters (default: 10.0)
- `--sea-level SEA_LEVEL`: Sea level in meters (default: 0.0)
- `--rho-ice RHO_ICE`: Ice density in kg/m³ (default: 910.0)
- `--rho-ocean RHO_OCEAN`: Ocean water density in kg/m³ (default: 1028.0)
- `--thickness-var THICKNESS_VAR`: Name of thickness variable (default: 'thickness')
- `--bed-var BED_VAR`: Name of bed topography variable (default: 'bedTopography')

## Example

```bash
python adjust_bed_to_haf.py \
    --mesh AIS_4to20km_mesh.nc \
    --geojson Thwaites_GL_2014_pinning_points.geojson \
    --output AIS_4to20km_mesh_adjusted.nc \
    --target-haf 15.0
```

## What is Height Above Flotation?

Height above flotation (HAF) measures how far ice is from hydrostatic equilibrium:
- **HAF = 0**: Ice at flotation (neutral buoyancy)
- **HAF > 0**: Ice grounded above flotation level
- **HAF < 0**: Ice below flotation (uncommon)

## Physics

The tool calculates required bed elevation using:

```
bed = sea_level + HAF_target - thickness × (ρ_ice / ρ_ocean)
```

This ensures the ice surface achieves the target HAF given the existing thickness.

## Input Requirements

### Mesh File
Must contain:
- `lonCell`, `latCell`: Cell coordinates (in radians)
- `thickness`: Ice thickness (m)
- `bedTopography`: Bed elevation (m)

### GeoJSON File
Must contain:
- Polygon or MultiLineString geometries
- WGS 84 coordinates (EPSG:4326, longitude/latitude in degrees)

## Dependencies

**Required:**
- xarray
- numpy
- shapely

**Optional:**
- geopandas (recommended for better performance)

## Notes

- The tool converts mesh coordinates from radians to degrees for comparison with GeoJSON
- Ice thickness is preserved; only bed topography is adjusted
- The output file includes updated `history` and `comment` attributes
- If no output file is specified, the input file is modified in place

## References

Wild, C. T., et al. (2022). Thwaites Glacier 2014 and 2019/20 grounding line positions.
U.S. Antarctic Program (USAP) Data Center. doi: 10.15784/601499
