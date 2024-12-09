import numpy as np


constants = \
    {'SHR_CONST_PI': 3.14159265358979323846,
     'SHR_CONST_CDAY': 86400.0,
     'SHR_CONST_SDAY': 86164.0,
     'SHR_CONST_REARTH': 6.37122e6,
     'SHR_CONST_G': 9.80616,
     'SHR_CONST_STEBOL': 5.67e-8,
     'SHR_CONST_BOLTZ': 1.38065e-23,
     'SHR_CONST_AVOGAD': 6.02214e26,
     'SHR_CONST_MWDAIR': 28.966,
     'SHR_CONST_MWWV': 18.016,
     'SHR_CONST_KARMAN': 0.4,
     'SHR_CONST_PSTD': 101325.0,
     'SHR_CONST_PDB': 0.0112372,
     'SHR_CONST_TKTRIP': 273.16,
     'SHR_CONST_TKFRZ': 273.15,
     'SHR_CONST_ZSRFLYR': 3.0,
     'SHR_CONST_RHOFW': 1.000e3,
     'SHR_CONST_RHOSW': 1.026e3,
     'SHR_CONST_RHOICE': 0.917e3,
     'SHR_CONST_CPDAIR': 1.00464e3,
     'SHR_CONST_CPWV': 1.810e3,
     'SHR_CONST_CPFW': 4.188e3,
     'SHR_CONST_CPSW': 3.996e3,
     'SHR_CONST_CPICE': 2.11727e3,
     'SHR_CONST_LATICE': 3.337e5,
     'SHR_CONST_LATVAP': 2.501e6,
     'SHR_CONST_CONDICE': 2.1,
     'SHR_CONST_TF0': 6.22e-2,
     'SHR_CONST_DTF_DP': -7.43e-8,
     'SHR_CONST_DTF_DS': -5.63e-2,
     'SHR_CONST_DTF_DPDS': -1.74e-10,
     'SHR_CONST_OCN_REF_SAL': 34.7,
     'SHR_CONST_ICE_REF_SAL':  4.0,
     'SHR_CONST_SPVAL': 1.0e30,
     'SHR_CONST_SPVAL_TOLMIN': 0.99 * 1.0e30,
     'SHR_CONST_SPVAL_TOLMAX': 1.01 * 1.0e30,
     'SHR_CONST_SPVAL_AERODEP': 1.e29,
     'SHR_CONST_VSMOW_18O': 2005.2e-6,
     'SHR_CONST_VSMOW_17O': 379.e-6,
     'SHR_CONST_VSMOW_16O': 0.997628,
     'SHR_CONST_VSMOW_D': 155.76e-6,
     'SHR_CONST_VSMOW_T': 1.85e-6,
     'SHR_CONST_VSMOW_H': 0.99984426,
     'SHR_CONST_RSTD_H2ODEV': 1.0}

constants['SHR_CONST_OMEGA'] = 2.0 * np.pi / constants['SHR_CONST_SDAY']
