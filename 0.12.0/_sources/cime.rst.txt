.. _cime_mod:

CIME Constants
==============

The module :py:mod:`mpas_tools.cime.constants` contains constants that are in
sync with `CIME <https://github.com/ESMCI/cime>`_, which provides infrastructure
and utilities for Earth System Models such at E3SM.  Currently, we sync only
those constants given numerical values in CIME, not those that are derivied
from other constants.  Constants are checked against their values on CIME's
master branch during tests of the conda build.  See
:py:func:`mpas_tools.tests.test_cime_constants.test_cime_constants`.

Some of the constants most likely to be useful in MPAS-Tools, COMPASS and other
related projects are:

* ``SHR_CONST_CDAY`` - sec in calendar day (s)
* ``SHR_CONST_REARTH`` - radius of Earth (m)
* ``SHR_CONST_G`` - acceleration of gravity (m/s^2)
* ``SHR_CONST_RHOFW`` - density of fresh water (kg/m^3)
* ``SHR_CONST_RHOSW`` - density of sea water (kg/m^3)
* ``SHR_CONST_RHOICE`` - density of ice (kg/m^3)
* ``SHR_CONST_CPFW`` - specific heat of fresh water (J/kg/K)
* ``SHR_CONST_CPSW`` - specific heat of sea water (J/kg/K)
* ``SHR_CONST_CPICE`` - specific heat of fresh ice (J/kg/K)
* ``SHR_CONST_LATICE`` - latent heat of fusion (J/kg)
* ``SHR_CONST_LATVAP`` - latent heat of evaporation (J/kg)
