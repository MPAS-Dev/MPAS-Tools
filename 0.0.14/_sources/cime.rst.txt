.. _cime_mod:

CIME Constants
==============

The module :py:mod:`mpas_tools.cime.constants` contains constants that are in
sync with `CIME <https://github.com/ESMCI/cime>`_, which provides infrastructure
and utilities for Earth System Models such at E3SM.  Currently, the only
constant being synched with CIME is the radius of the Earth.  Constants are
checked against their values on CIME's master branch during tests of the
conda build.  See
:py:func:`mpas_tools.tests.test_cime_constants.test_cime_constants`.
