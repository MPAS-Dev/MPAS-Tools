.. _logging:

*******
Logging
*******

MPAS-Tools includes capabilities for logging output from many of its function
calls to either ``stdout``/``stderr`` or to a log file.

Using logging
=============

Logging is performed by creating a
:py:class:`mpas_tools.logging.LoggingContext` object inside a ``with``
statement:

.. code-block:: python

    with LoggingContext(__name__) as logger:

        ...
        logger.info('Step 1. Generate mesh with JIGSAW')
        ...
        logger.info('Step 2. Convert triangles from jigsaw format to netcdf')
        ...
        logger.info('Step 3. Convert from triangles to MPAS mesh')
        ...

A ``LoggingContext`` is given a unique name (typically the name of the module,
``__name__``).  If no other arguments are supplied, the logger will write to
``stdout`` with calls to ``logger.info()`` and ``logger.debug()``, and to
``stderr`` for calls to ``logger.error()``.

It is convenient to create this kind of a ``logger`` in contexts where logging
might be to a file or might be to ``stdout``/``stderr``, depending on the
calling code.

Often, a function should support output to a ``logger`` but you do not
necessarily want the calling code to have to crate one if the calling code just
wants output to ``stdout``/``stderr``.  In such contexts, the function can take
``logger=None`` as an optional argument.  Then, in the function itself, a new
``LoggingContext`` can be created that will just use the existing ``logger``
if there is one or create a new ong for ``stdout``/``stderr`` if none is
provided.  As an example, here is a snippet from
:py:func:`mpas_tools.mesh.creation.build_mesh.build_spherical_mesh()`:

.. code-block:: python

    def build_spherical_mesh(cellWidth, lon, lat, earth_radius,
                             out_filename='base_mesh.nc', plot_cellWidth=True,
                             dir='./', logger=None):
        ...
        with LoggingContext(__name__, logger=logger) as logger:

            ...
            logger.info('Step 1. Generate mesh with JIGSAW')
            jigsaw_driver(cellWidth, lon, lat, on_sphere=True,
                          earth_radius=earth_radius, logger=logger)

            logger.info('Step 2. Convert triangles from jigsaw format to netcdf')
            jigsaw_to_netcdf(msh_filename='mesh-MESH.msh',
                             output_name='mesh_triangles.nc', on_sphere=True,
                             sphere_radius=earth_radius)

            logger.info('Step 3. Convert from triangles to MPAS mesh')
            write_netcdf(convert(xarray.open_dataset('mesh_triangles.nc'), dir=dir,
                                 logger=logger),
                         out_filename)

The optional argument ``logger`` is passed to a new ``LoggingContext``.  That
way ``logger`` is guaranteed not to be ``None`` (so calls to ``logger.info``
work properly).  The ``logger`` is also passed on to
:py:func:`mpas_tools.mesh.creation.jigsaw_driver.jigsaw_driver()` and
:py:func:`mpas_tools.io.write_netcdf()` for output from those functions.

To log to a new log file, simply pass a file name to the ``log_filename``
argument:

.. code-block:: python

    with LoggingContext(name=__name__, log_filename='output.log') as logger:
        ...
        logger.info('Step 1. Generate mesh with JIGSAW')


If both the ``logger`` and ``log_filename`` arguments are not ``None``, the
``log_filename`` is ignored and the existing ``logger`` is simply used for
logging.

Logging subprocess calls
========================

You can also run subprocesses and capture the output to a ``logger``.  This is
accomplished with the function :py:func:`mpas_tools.logging.check_call`, which
acts a lot like :py:func:`subprocess.check_call()` but with output going to
the logger, as in this example from
:py:func:`mpas_tools.mesh.creation.jigsaw_driver.jigsaw_driver()`

.. code-block:: python

    from mpas_tools.logging import check_call


    def jigsaw_driver(cellWidth, x, y, on_sphere=True, earth_radius=6371.0e3,
                      geom_points=None, geom_edges=None, logger=None):

        ...
        opts.jcfg_file = 'mesh.jig'
        ...
        check_call(['jigsaw', opts.jcfg_file], logger=logger)
