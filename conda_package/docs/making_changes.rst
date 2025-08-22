.. _dev_making_changes:

****************************
Making Changes to mpas_tools
****************************

New python functions and modules (``.py`` files) can be added within the
``conda_package/mpas_tools``.  These will automatically be part of the
``mpas_tools`` package.  New directories with python modules should include an
``__init__.py`` file (which can be empty) to indicate that they are also part of
the package.

Code Styling and Linting
========================
``mpas_tools`` utilizes ``pre-commit`` to lint incoming code when you make a
commit (as long as you have your environment set up correctly), and on GitHub
whenever you make a pull request to the ``mpas_tools`` repository. Linting makes sure
your code follows the formatting guidelines of PEP8, and cleans up additional
things like whitespace at the end of lines.

The first time you set up the ``mpas_tools_dev`` environment, you will need to set up
``pre-commit``. This is done by running:

```bash
pre-commit install
```

You only need to do this once when you create the ``mpas_tools_dev``
environment. If you create a new version of ``mpas_tools_dev``, then you will
need to run it again.

When you run ``git commit <filename>``, ``pre-commit`` will automatically lint
your code before committing. Some formatting will be updated by ``pre-commit``
automatically, in which case it will terminate the commit and inform you of the
change. Then you can run ``git commit <filename>`` again to continue the
linting process until your commit is successful. Some changes need to be made
manually, such as a line being too long. When this happens, you must update the
file to ``pre-commit``'s standards, and then attempt to re-commit the file.

Internally, ``pre-commit``  uses `ruff <https://docs.astral.sh/ruff/>` to check
PEP8 compliance, as well as sort, check and format imports,
`flynt <https://github.com/ikamensh/flynt>` to change any format strings to
f-strings, and `mypy <https://mypy-lang.org/>` to check for consistent variable
types. An example error might be:

```bash
example.py:77:1: E302 expected 2 blank lines, found 1
```

For this example, we would just add an additional blank line after line 77 and
try the commit again to make sure we've resolved the issue.

You may also find it useful to use an IDE with a PEP8 style checker built in,
such as `VS Code <https://code.visualstudio.com/>`. See
`Formatting Python in VS Code <https://code.visualstudio.com/docs/python/formatting>`
for some tips on checking code style in VS Code.

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
    python>=3.10
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
