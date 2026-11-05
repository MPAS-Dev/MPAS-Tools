"""Small helpers for writing ISMIP7 NetCDF files."""

import numpy as np

# values at or beyond this magnitude are treated as missing (ESMF/ncremap
# and MALI fill values are 1e34 to 1e37)
_FILL_THRESHOLD = 1.0e30


def to_masked(values):
    """
    Return ``values`` as a masked array in which NaN, infinite and fill-like
    (|x| > 1e30) entries are masked.  Assigning a masked array to a netCDF4
    variable writes the variable's _FillValue in the masked cells, as the
    ISMIP7 compliance checker requires for missing data (NaN is not
    accepted).
    """
    values = np.ma.asarray(values)
    bad = ~np.isfinite(values.filled(0.0)) | \
        (np.abs(values.filled(0.0)) > _FILL_THRESHOLD)
    return np.ma.masked_where(bad | np.ma.getmaskarray(values), values)


def zero_where_invalid(values):
    """
    Return ``values`` with NaN, infinite and fill-like (|x| > 1e30) entries
    replaced by 0.0.  The ISMIP7 data request does not permit missing values
    for some variables (e.g., ice thickness, area fractions and fluxes): where
    there is no ice the value is 0, not missing.
    """
    values = np.ma.filled(np.ma.asarray(values, dtype=float), np.nan)
    bad = ~np.isfinite(values) | (np.abs(np.nan_to_num(values)) > _FILL_THRESHOLD)
    return np.where(bad, 0.0, values)
