#!/bin/bash
# Test script for adjust_bed_to_haf.py

# Test with the provided files using the standalone script
python3 landice/mesh_tools_li/adjust_bed_to_haf.py \
    --mesh /Users/trhille/Documents/ISMIP7/Antarctica/mesh/4to20km/AIS_4to20km_r03_20260910_ASE_extracted.nc \
    --geojson /Users/trhille/Documents/ISMIP7/Antarctica/wild_grounding_lines/Thwaites_GLs_2014_201920/Thwaites_GL_2014_pinning_points.geojson \
    --projection ais-bedmap2 \
    --output test_output.nc \
    --target-haf 10.0

echo ""
echo "Test completed. Output saved to test_output.nc"
echo "You can compare the original and modified bed topography using ncview or similar tools."
