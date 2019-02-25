#!/usr/bin/env python

from setuptools import setup, find_packages

version = '0.0.1'

setup(name='mpas_mesh_tools',
      version=version,
      description='A set of tools for creating and manipulating meshes for the'
                  ' climate components based on the Model for Prediction '
                  'Across Scales (MPAS) framework',
      url='https://github.com/MPAS-Dev/MPAS-Tools/tree/master/mesh_tools',
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
      install_requires=['numpy', 'xarray', 'netCDF4'],
      entry_points={'console_scripts':
                    ['planar_hex = mpas_mesh_tools.planar_hex:main']})
