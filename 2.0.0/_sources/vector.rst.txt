.. _vector:

*****************
Vector Operations
*****************

MPAS-Tools has a ``Vector`` class and an unrelated tool for performing
vector reconstructions at cell centers from fields defined at edge normals.

.. _vector_class:

Vector Class
============

The :py:class:`mpas_tools.vector.Vector` class defines a single vector (with
components that are floats) or a vector field (with components that are
:py:class:`numpy.ndarray` objects).  See the API documentation for the
individual methods to find out more.


.. _vector_reconstruct:

Vector Reconstructions
======================

The command-line tool ``vector_reconstruct`` and the function
:py:func:`mpas_tools.vector.reconstruct.reconstruct_variable()` are used to
reconstruct Cartesian (X, Y, Z), zonal and meridional components of an MPAS
vector field at cell centers, given the field on edge normals.

This tool requires that the field ``coeffs_reconstruct`` has been saved to a
NetCDF file.  The simplest way to do this is to include the following
stream in a forward run:

.. code-block:: xml

    <stream name="vector_reconstruction"
            clobber_mode="truncate"
            type="output"
            output_interval="initial_only"
            filename_template="vector_reconstruction.nc">

        <var name="coeffs_reconstruct"/>
    </stream>

and run the model for one time step.

Then, ``vector_reconstruct`` is called with:

.. code-block::

    $ vector_reconstruct --help
    usage: vector_reconstruct [-h] [-m MESH_FILENAME] [-w WEIGHTS_FILENAME] -i
                              IN_FILENAME -v VARIABLES [VARIABLES ...]
                              [--out_variables OUT_VARIABLES [OUT_VARIABLES ...]]
                              -o OUT_FILENAME

You must supply input and output files and a list of one or more variables on
edge normals to reconstruct.  You can optionally supply a separate file that
contains the MPAS mesh if it is not part of the input file, a file with
``coeffs_reconstruct`` if it is not in the input file, and a list of variable
prefixes corresponding to each variable supplied that should be prepended to
the Cartesian, zonal and meridional reconstructions for that variable.
