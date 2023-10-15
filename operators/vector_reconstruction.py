#!/usr/bin/env python

"""
Extract Cartesian (X, Y, Z), zonal and meridional components of an MPAS vector
field, given the field on edge normals.

This tool requires that the field 'coeffs_reconstruct' has been saved to a
NetCDF file.  The simplest way to do this is to include the following stream
in a forward run:

<stream name="vector_reconstruction"
        clobber_mode="truncate"
        type="output"
        output_interval="initial_only"
        filename_template="vector_reconstruction.nc">

    <var name="coeffs_reconstruct"/>
</stream>

and run the model for one time step.

"""
from mpas_tools.vector.reconstruct import main


if __name__ == '__main__':
    main()
