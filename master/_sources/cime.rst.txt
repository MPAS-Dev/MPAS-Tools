.. _cime_mod:

Physical Constants
==================

E3SM takes its physical constants from the
`Physical Constants Dictionary (PCD) <https://icams-us.github.io/PhysicalConstantsDictionary/>`_,
a YAML file found at ``share/pcd.yaml`` in E3SM.  MPAS-Tools packages a copy of
this file in :py:mod:`mpas_tools.constants`.  During tests of the rattler-build
package build, the copy is checked against the one on E3SM's ``master``
branch.

.. _pcd_constants:

PCD constants
-------------

Use :py:func:`mpas_tools.constants.get_constant` to get a constant by its PCD
name:

.. code-block:: python

    from mpas_tools.constants import get_constant

    earth_radius = get_constant('mean_radius')

Some of the constants are given by a ``formula`` in terms of other constants,
rather than a ``value``.  These are evaluated in double precision, in the same
way as in the compiled models.
:py:func:`mpas_tools.constants.get_pcd_version` gives the version of the
packaged PCD.

.. _cime_constants:

CIME constants
--------------

The dictionary ``constants`` in :py:mod:`mpas_tools.cime.constants` contains
constants with the names and values of those in E3SM's
``share/util/shr_const_mod.F90``.  E3SM sets most of these from the PCD, and
MPAS-Tools does the same.  A few, such as ``SHR_CONST_SDAY`` and
``SHR_CONST_KARMAN``, are still set in ``shr_const_mod`` itself.  The constants
that ``shr_const_mod`` derives from other constants (e.g. ``SHR_CONST_RGAS``)
are not included.  During tests of the rattler-build package build, each
constant is checked against ``shr_const_mod`` on E3SM's ``master`` branch.

Some of the constants most likely to be useful in MPAS-Tools, Polaris and
other related projects are (with their PCD names):

* ``SHR_CONST_CDAY`` - sec in calendar day (s), ``seconds_per_day``
* ``SHR_CONST_REARTH`` - radius of Earth (m), ``mean_radius``
* ``SHR_CONST_G`` - acceleration of gravity (m/s^2),
  ``standard_acceleration_of_gravity``
* ``SHR_CONST_RHOFW`` - density of fresh water (kg/m^3),
  ``pure_water_density_reference``
* ``SHR_CONST_RHOSW`` - density of sea water (kg/m^3),
  ``seawater_density_reference``
* ``SHR_CONST_RHOICE`` - density of ice (kg/m^3),
  ``sea_ice_density_reference``
* ``SHR_CONST_CPFW`` - specific heat of fresh water (J/kg/K),
  ``pure_water_specific_heat_capacity_reference``
* ``SHR_CONST_CPSW`` - specific heat of sea water (J/kg/K),
  ``seawater_specific_heat_capacity_reference``
* ``SHR_CONST_CPICE`` - specific heat of fresh ice (J/kg/K),
  ``sea_ice_specific_heat_capacity_reference``
* ``SHR_CONST_LATICE`` - latent heat of fusion (J/kg),
  ``latent_heat_of_fusion_reference``
* ``SHR_CONST_LATVAP`` - latent heat of evaporation (J/kg),
  ``latent_heat_of_vaporization_reference``

Updating the PCD
----------------

When E3SM updates its PCD, copy E3SM's ``share/pcd.yaml`` to
``conda_package/mpas_tools/constants/pcd.yaml``.  If E3SM has also changed
which PCD constants ``shr_const_mod.F90`` uses, update
``conda_package/mpas_tools/cime/constants.py`` to match.
