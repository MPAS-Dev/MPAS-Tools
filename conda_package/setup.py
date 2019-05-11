#!/usr/bin/env python

from setuptools import setup, find_packages

version = '0.0.2'

setup(name='mpas_tools',
      version=version,
      description='A set of tools for creating and manipulating meshes for the'
                  ' climate components based on the Model for Prediction '
                  'Across Scales (MPAS) framework',
      url='https://github.com/MPAS-Dev/MPAS-Tools',
      author='MPAS-Analysis Developers',
      author_email='mpas-developers@googlegroups.com',
      license='BSD',
      classifiers=[
          'Development Status :: 3 - Alpha',
          'License :: OSI Approved :: BSD License',
          'Operating System :: OS Independent',
          'Intended Audience :: Science/Research',
          'Programming Language :: Python',
          'Programming Language :: Python :: 2',
          'Programming Language :: Python :: 2.7',
          'Programming Language :: Python :: 3',
          'Programming Language :: Python :: 3.6',
          'Programming Language :: Python :: 3.7',
          'Topic :: Scientific/Engineering',
      ],
      packages=find_packages(),
      package_data={},
      scripts=['mesh_tools/mesh_conversion_tools/mark_horns_for_culling.py',
               'landice/mesh_tools_li/create_landice_grid_from_generic_MPAS_grid.py',
               'landice/mesh_tools_li/define_cullMask.py',
               'landice/mesh_tools_li/interpolate_to_mpasli_grid.py',
               'landice/mesh_tools_li/mark_domain_boundaries_dirichlet.py',
               'ocean/coastline_alteration/add_land_locked_cells_to_mask.py',
               'ocean/coastline_alteration/widen_transect_edge_masks.py',
               'ocean/coastline_alteration/add_critical_land_blockages_to_mask.py',
               'ocean/coastline_alteration/add_critical_land_blockages_to_mask.py',
               'visualization/paraview_vtk_field_extractor/paraview_vtk_field_extractor.py'],
      install_requires=['numpy', 'xarray', 'netCDF4', 'pyevtk'],
      entry_points={'console_scripts':
                    ['planar_hex = mpas_tools.planar_hex:main',
                     'translate_planar_grid = mpas_tools.translate:main']})
