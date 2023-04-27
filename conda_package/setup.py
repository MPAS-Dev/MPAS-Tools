#!/usr/bin/env python

import os
import re
from setuptools import setup, find_packages
import shutil


install_requires = [
    'cartopy',
    'cmocean',
    'dask',
    'inpoly',
    'matplotlib',
    'netcdf4',
    'numpy',
    'progressbar2',
    'pyamg',
    'pyevtk',
    'pyproj',
    'igraph',
    'scikit-image',
    'scipy',
    'shapely >=2.0,<3.0',
    'xarray']

here = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(here, 'mpas_tools', '__init__.py')) as f:
    init_file = f.read()

version = re.search(r'{}\s*=\s*[(]([^)]*)[)]'.format('__version_info__'),
                    init_file).group(1).replace(', ', '.')

os.chdir(here)

for path in ['ocean', 'landice', 'visualization', 'mesh_tools']:
    if not os.path.exists(path):
        shutil.copytree('../{}'.format(path), './{}'.format(path))

setup(name='mpas_tools',
      version=version,
      description='A set of tools for creating and manipulating meshes for the'
                  ' climate components based on the Model for Prediction '
                  'Across Scales (MPAS) framework',
      url='https://github.com/MPAS-Dev/MPAS-Tools',
      author='MPAS-Tools Developers',
      author_email='mpas-developers@googlegroups.com',
      license='BSD',
      classifiers=[
          'Development Status :: 5 - Production/Stable',
          'License :: OSI Approved :: BSD License',
          'Operating System :: OS Independent',
          'Intended Audience :: Science/Research',
          'Programming Language :: Python',
          'Programming Language :: Python :: 3',
          'Programming Language :: Python :: 3.8',
          'Programming Language :: Python :: 3.9',
          'Programming Language :: Python :: 3.10',
          'Programming Language :: Python :: 3.11',
          'Topic :: Scientific/Engineering',
      ],
      packages=find_packages(include=['mpas_tools', 'mpas_tools.*']),
      package_data={'mpas_tools.viz': ['SciVisColorColormaps/*.xml']},
      scripts=['mesh_tools/mesh_conversion_tools/mark_horns_for_culling.py',
               'mesh_tools/planar_grid_transformations/set_lat_lon_fields_in_planar_grid.py',
               'mesh_tools/create_SCRIP_files/create_SCRIP_file_from_MPAS_mesh.py',
               'mesh_tools/create_SCRIP_files/create_SCRIP_file_from_planar_rectangular_grid.py',
               'landice/mesh_tools_li/create_landice_grid_from_generic_MPAS_grid.py',
               'landice/mesh_tools_li/define_cullMask.py',
               'landice/mesh_tools_li/interpolate_to_mpasli_grid.py',
               'landice/mesh_tools_li/mark_domain_boundaries_dirichlet.py',
               'ocean/coastline_alteration/add_land_locked_cells_to_mask.py',
               'ocean/coastline_alteration/widen_transect_edge_masks.py',
               'ocean/coastline_alteration/add_critical_land_blockages_to_mask.py',
               'ocean/moc_southern_boundary_extractor/moc_southern_boundary_extractor.py',
               'visualization/paraview_vtk_field_extractor/paraview_vtk_field_extractor.py'],
      install_requires=install_requires,
      entry_points={'console_scripts':
                    ['planar_hex = mpas_tools.planar_hex:main',
                     'translate_planar_grid = mpas_tools.translate:main',
                     'merge_grids = mpas_tools.merge_grids:main',
                     'split_grids = mpas_tools.split_grids:main',
                     'inject_bathymetry = mpas_tools.ocean.inject_bathymetry:main',
                     'inject_preserve_floodplain = mpas_tools.ocean.inject_preserve_floodplain:main',
                     'ocean_add_depth = mpas_tools.ocean.depth:main_add_depth',
                     'ocean_add_zmid = mpas_tools.ocean.depth:main_add_zmid',
                     'ocean_write_time_varying_zmid = mpas_tools.ocean.depth:main_write_time_varying_zmid',
                     'plot_ocean_transects = mpas_tools.ocean.viz.transects:main',
                     'mpas_to_triangle = mpas_tools.mesh.creation.mpas_to_triangle:main',
                     'triangle_to_netcdf = mpas_tools.mesh.creation.triangle_to_netcdf:main',
                     'jigsaw_to_netcdf = mpas_tools.mesh.creation.jigsaw_to_netcdf:main',
                     'scrip_from_mpas = mpas_tools.scrip.from_mpas:main',
                     'compute_mpas_region_masks = mpas_tools.mesh.mask:entry_point_compute_mpas_region_masks',
                     'compute_mpas_transect_masks = mpas_tools.mesh.mask:entry_point_compute_mpas_transect_masks',
                     'compute_mpas_flood_fill_mask = mpas_tools.mesh.mask:entry_point_compute_mpas_flood_fill_mask',
                     'compute_lon_lat_region_masks = mpas_tools.mesh.mask:entry_point_compute_lon_lat_region_masks',
                     'compute_projection_region_masks = mpas_tools.mesh.mask:entry_point_compute_projection_region_masks',
                     'prepare_seaice_partitions = mpas_tools.seaice.partition:prepare_partitions',
                     'create_seaice_partitions = mpas_tools.seaice.partition:create_partitions',
                     'simple_seaice_partitions = mpas_tools.seaice.partition:simple_partitions']})
