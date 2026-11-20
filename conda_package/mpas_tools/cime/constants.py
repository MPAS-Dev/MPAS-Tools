"""
Constants with the names and values of those in E3SM's ``shr_const_mod``

E3SM sets most of the constants in ``share/util/shr_const_mod.F90`` from the
Physical Constants Dictionary (PCD).  Here, they are taken from the same PCD
in :py:mod:`mpas_tools.constants.pcd`, in the same way, so that their values
match E3SM's.  New code should prefer
:py:func:`mpas_tools.constants.get_constant` with PCD names.
"""

from mpas_tools.constants.pcd import get_constant

# the PCD name of each constant that E3SM takes directly from the PCD
_PCD_NAMES = {
    'SHR_CONST_PI': 'pi',
    'SHR_CONST_CDAY': 'seconds_per_day',
    'SHR_CONST_OMEGA': 'angular_velocity',
    'SHR_CONST_REARTH': 'mean_radius',
    'SHR_CONST_G': 'standard_acceleration_of_gravity',
    'SHR_CONST_STEBOL': 'stefan_boltzmann_constant',
    'SHR_CONST_BOLTZ': 'boltzmann_constant',
    'SHR_CONST_PSTD': 'standard_atmosphere',
    'SHR_CONST_PDB': 'pee_dee_belemnite_ratio_13c_to_12c',
    'SHR_CONST_TKTRIP': 'water_triple_point_temperature',
    'SHR_CONST_TKFRZ': 'pure_water_freezing_temperature_reference',
    'SHR_CONST_RHOFW': 'pure_water_density_reference',
    'SHR_CONST_RHOSW': 'seawater_density_reference',
    'SHR_CONST_RHOICE': 'sea_ice_density_reference',
    'SHR_CONST_CPDAIR': 'dry_air_specific_heat_capacity_reference',
    'SHR_CONST_CPWV': 'water_vapor_specific_heat_capacity_reference',
    'SHR_CONST_CPFW': 'pure_water_specific_heat_capacity_reference',
    'SHR_CONST_CPSW': 'seawater_specific_heat_capacity_reference',
    'SHR_CONST_CPICE': 'sea_ice_specific_heat_capacity_reference',
    'SHR_CONST_LATICE': 'latent_heat_of_fusion_reference',
    'SHR_CONST_LATVAP': 'latent_heat_of_vaporization_reference',
    'SHR_CONST_CONDICE': 'sea_ice_thermal_conductivity_reference',
    'SHR_CONST_OCN_REF_SAL': 'ocean_reference_salinity',
    'SHR_CONST_ICE_REF_SAL': 'sea_ice_reference_salinity',
    'SHR_CONST_VSMOW_18O': 'vsmow_ratio_18o_to_16o',
    'SHR_CONST_VSMOW_17O': 'vsmow_ratio_17o_to_16o',
    'SHR_CONST_VSMOW_16O': 'vsmow_fraction_16o_of_total_oxygen',
    'SHR_CONST_VSMOW_D': 'vsmow_ratio_2h_to_1h',
    'SHR_CONST_VSMOW_T': 'vsmow_ratio_3h_to_1h',
    'SHR_CONST_VSMOW_H': 'vsmow_fraction_1h_of_total_hydrogen',
}

# the PCD name of each constant that E3SM converts from per mole in the PCD
# to per kilomole
_PCD_NAMES_PER_KMOLE = {
    'SHR_CONST_AVOGAD': 'avogadro_constant',
    'SHR_CONST_MWDAIR': 'dry_air_molar_mass_reference',
    'SHR_CONST_MWWV': 'water_vapor_molar_mass_reference',
}

constants = {
    name: get_constant(pcd_name) for name, pcd_name in _PCD_NAMES.items()
}

for name, pcd_name in _PCD_NAMES_PER_KMOLE.items():
    constants[name] = 1000.0 * get_constant(pcd_name)

# constants that E3SM still sets in shr_const_mod rather than from the PCD
constants.update(
    {
        'SHR_CONST_SDAY': 86164.0,
        'SHR_CONST_KARMAN': 0.4,
        'SHR_CONST_ZSRFLYR': 3.0,
        'SHR_CONST_TF0': 6.22e-2,
        'SHR_CONST_DTF_DP': -7.43e-8,
        'SHR_CONST_DTF_DS': -5.63e-2,
        'SHR_CONST_DTF_DPDS': -1.74e-10,
        'SHR_CONST_SPVAL': 1.0e30,
        'SHR_CONST_SPVAL_TOLMIN': 0.99 * 1.0e30,
        'SHR_CONST_SPVAL_TOLMAX': 1.01 * 1.0e30,
        'SHR_CONST_SPVAL_AERODEP': 1.0e29,
        'SHR_CONST_RSTD_H2ODEV': 1.0,
    }
)
