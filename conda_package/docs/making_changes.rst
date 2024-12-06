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
both in ``conda_package/pyproject.toml`` and ``conda_package/recipe/meta.yaml``.

As an example, the entry point ``planar_hex`` is defined in ``pyproject.toml`` as:

.. code-block:: python

  [project.scripts]
  ...
  planar_hex = "mpas_tools.planar_hex:main"
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

If you changes introduce new dependencies, these need to be added to both
the recipe for the conda package in ``conda_package/recipe/meta.yaml`` and
to the text file describing the development environment,
``conda_package/dev-spec.txt``.

In ``meta.yaml``, add these changes in alphabetical order to the ``run``
section of ``requirements``:

.. code-block:: yaml

  requirements:
  ...
    run:
      - python
      - affine
      ...

These requirements *must* be on the ``conda-forge`` anaconda channel.  If you
need help with this, please contact the developers.

Add the new dependencies in alphabetical order to ``dev-speck.txt``
under the ``# Base`` comment:

.. code-block:: none

    ...
    # This file may be used to create an environment using:
    # $ conda create --name <env> --file <this file>

    # Base
    python>=3.9
    cartopy
    ...

Updating the Version
====================

Before a release of the package, the version of ``mpas_tools`` needs to be
updated in 3 places.  First, in ``conda_package/mpas_tools/__init__.py``:

.. code-block:: python

  __version_info__ = (0, 6, 0)
  __version__ = '.'.join(str(vi) for vi in __version_info__)

Increment ``__version_info__`` (major, minor or micro version, depending on
what makes sense).

Second, the version in the conda recipe (``conda_package/recipe/meta.yaml``)
needs to match:

.. code-block::

  {% set name = "mpas_tools" %}
  {% set version = "0.6.0" %}

Third, Add the new version to the :ref:`versions` in the documentation.

.. code-block::

    `v0.6.0`_         `0.6.0`_
    ================ ===============

    ...

    .. _`v0.6.0`: ../0.6.0/index.html
    .. _`0.6.0`: https://github.com/MPAS-Dev/MPAS-Analysis/tree/0.6.0


The new links won't be valid until a new release is made and Azure Pipelines
has generated the associated documentation.  Eventually, it should be possible
to do this automatically but that has not yet been implemented.
