# ======== DEFINE PROJECTIONS =============
# Create empty dictionary to store projection definitions:
projections = dict()
# add more as needed:

# CISM's projection is as follows, with the vertical datum as EIGEN-GL04C
# geoid. Datum is actually EIGEN-GL04C but that is not an option in Proj.
# Therefore using EGM08 which should be within ~1m everywhere (and 10-20
# cm in most places)
# NOTE!!!!!!  egm08_25.gtx can be downloaded from:
# http://download.osgeo.org/proj/vdatum/egm08_25/egm08_25.gtx  and the
# path in the projection specification line should point to it!

# projections['gis-bamber'] = \
#     '+proj=stere +lat_ts=71.0 +lat_0=90 +lon_0=321.0 +k_0=1.0 '
#     '+x_0=800000.0 +y_0=3400000.0 +ellps=WGS84 ' \
#     '+geoidgrids=./egm08_25.gtx'

# This version ignores the vertical datum shift, which should be a very
# small error for horizontal-only positions
projections['gis-bamber'] = (
    '+proj=stere +lat_ts=71.0 +lat_0=90 +lon_0=321.0 +k_0=1.0 '
    '+x_0=800000.0 +y_0=3400000.0 +ellps=WGS84'
)

# GIMP projection: This is also polar stereographic but with different
# standard parallel and using the WGS84 ellipsoid.
projections['gis-gimp'] = (
    '+proj=stere +lat_ts=70.0 +lat_0=90 +lon_0=315.0 +k_0=1.0 +x_0=0.0 '
    '+y_0=0.0 +ellps=WGS84'
)

# BEDMAP2 projection
# Note: BEDMAP2 elevations use EIGEN-GL04C geoid
projections['ais-bedmap2'] = (
    '+proj=stere +lat_ts=-71.0 +lat_0=-90 +lon_0=0.0 +k_0=1.0 +x_0=0.0 '
    '+y_0=0.0 +ellps=WGS84'
)

# BEDMAP2 projection of sphere. This projection must be used to adjust MALI
# mesh when performing coupled MALI-SeaLevelModel simulations. Otherwise, ice
# mass won't be conserved between the MALI planar mesh and the spherical
# sea-level model grid during the post-processing (output analysis) step.
projections['ais-bedmap2-sphere'] = (
    '+proj=stere +lat_ts=-71.0 +lat_0=-90 +lon_0=0.0 +k_0=1.0 +x_0=0.0 '
    '+y_0=0.0 +ellps=sphere'
)

# Standard Lat/Long
projections['latlon'] = '+proj=longlat +ellps=WGS84'
