.. _dev_making_changes:

****************************
Making Changes to mpas_tools
****************************

New python functions and modules (``.py`` files) can be added within the
``conda_package/mpas_tools``.  These will automatically be part of the
``mpas_tools`` package.  New directories with python modules should include an
``__init__.py`` file (which can be empty) to indicate that they are also part of
the package.

Entry Points
============

The best way to add new "scripts" to the package is to add a function without
any arguments somewhere in the package, and then to add it as an "entry point"
both in ``conda_package/setup.py`` and ``conda_package/recipe/meta.yaml``.

As an example, the entry point ``planar_hex`` is defined in ``setup.py`` as:

.. code-block:: python

  setup(name='mpas_tools',
  ...
        entry_points={'console_scripts':
                    ['planar_hex = mpas_tools.planar_hex:main',
  ...

and in ``meta.yaml`` as:

.. code-block::

  build:
    number: 0
    entry_points:
      - planar_hex = mpas_tools.planar_hex:main

When the package is installed in a conda environment, a stub script
``planar_hex`` will be in the user's path that will call the function ``main()``
in the module ``mpas_tools.planar_hex``:

.. code-block:: python

  def main():

      parser = argparse.ArgumentParser(
          description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
      parser.add_argument('--nx', dest='nx', type=int, required=True,
                          help='Cells in x direction')
      ...
      args = parser.parse_args()

      make_planar_hex_mesh(args.nx, args.ny, args.dc,
                           args.nonperiodic_x, args.nonperiodic_y,
                           args.outFileName)

As you can see, the function pointed to by the entry point can used to parse
command-line arguments, just as a "normal" python script would do

By convention, entry points do not typically include the ``.py`` extension.

Dependencies
============

If you changes introduce new dependencies, these need to be added to the recipe
for the conda package in ``conda_package/recipe/meta.yaml``

Add these changes to the end of the ``run`` section of ``requirements``:

.. code-block::

  requirements:
  ...
    run:
      - python
      - netcdf4
      ...
      - affine

These requirements *must* be on the ``conda-forge`` anaconda channel.  If you
need help with this, please contact the developers.

