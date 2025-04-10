.. _config_mod:

Config files
============

The ``mpas_tools.config`` module includes the
:py:class:`mpas_tools.config.MpasConfigParser` class reading, getting, setting,
and writing config options and config files.

The :py:meth:`mpas_tools.config.MpasConfigParser.add_from_package()` method can
be used to add the contents of a config file within a package to the config
options.

Here is an example from `compass <https://github.com/MPAS-Dev/compass>`_
.. code-block:: python

    self.config.add_from_package(
       'compass.ocean.tests.global_ocean.make_diagnostics_files',
       'make_diagnostics_files.cfg', exception=True)

The first and second arguments are the name of a package containing the config
file and the name of the config file itself, respectively.  You can see that
the file is in the path ``compass/ocean/tests/global_ocean/make_diagnostics_files``
(replacing the ``.`` in the module name with ``/``).  In this case, we know
that the config file should always exist, so we would like the code to raise
an exception (``exception=True``) if the file is not found.  This is the
default behavior.  In some cases, you would like the code to add the config
options if the config file exists and do nothing if it does not
(``exception=False``).

Ihe ``MpasConfigParser`` class also includes methods for adding a user
config file, :py:meth:`mpas_tools.config.MpasConfigParser.add_user_config()`,
and other config files by file name,
:py:meth:`mpas_tools.config.MpasConfigParser.add_from_file()`.

The :py:meth:`mpas_tools.config.MpasConfigParser.copy()` method can be used to
make a deep copy of the config parser.  This is useful in cases where config
options should be added or modified without affecting the original config
object.  For example, this feature is used in MPAS-Analysis to set a reference
year as the start year in some analysis without affecting the start year in
other analysis.

The :py:meth:`mpas_tools.config.MpasConfigParser.set()` method has some
optional arguments not present in :py:class:`configparser.ConfigParser.set()`.
The ``comment`` argument can be used to add a comment that will be written
out above the config option.  The comment can cover multiple lines by including
a ``\n`` character.  The comment should not include the ``#`` comment
character, as this is added automatically.  The argument ``user=True`` can be
used to set "user" config options similar to reading a user config file with
:py:meth:`mpas_tools.config.MpasConfigParser.add_user_config()`.

Other methods for the ``MpasConfigParser`` are similar to those for
:py:class:`configparser.ConfigParser`.  In addition to ``get()``,
``getinteger()``, ``getfloat()`` and ``getboolean()`` methods, this class
implements :py:meth:`mpas_tools.config.MpasConfigParser.getlist()`, which
can be used to parse a config value separated by spaces and/or commas into
a list of strings, floats, integers, booleans, etc.  Another useful method
is :py:meth:`mpas_tools.config.MpasConfigParser.getexpression()`, which can
be used to get python dictionaries, lists and tuples as well as a small set
of functions (``range()``, :py:meth:`numpy.linspace()`,
:py:meth:`numpy.arange()`, and :py:meth:`numpy.array()`)

Currently, ``MpasConfigParser`` supports accessing a config section using
section names as keys, e.g.:

.. code-block:: python

    section = self.config['enthalpy_benchmark_viz']
    display_image = section.getboolean('display_image')
    ...

But it does not allow assignment of a section or many of the other
dictionary-like features supported by :py:class:`configparser.ConfigParser`.

.. _config_comments:

Comments in config files
------------------------

One of the main advantages of :py:class:`mpas_tools.config.MpasConfigParser`
over :py:class:`configparser.ConfigParser` is that it keeps track of comments
that are associated with config sections and options.  There are a few "rules"
that make this possible.

Comments must begin with the ``#`` character.  They must be placed *before* the
config section or option in question (preferably without blank lines between).
The comments can be any number of lines.

.. note::

    Inline comments (after a config option on the same line) are not allowed
    and will be parsed as part of the config option itself.

