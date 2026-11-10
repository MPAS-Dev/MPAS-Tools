"""
Calendar-aware time helpers for MALI output.

MALI daysSinceStart is relative to simulationStartTime and uses the calendar
given by the config_calendar_type global attribute ('gregorian' or 'noleap').
"""

import cftime
import numpy as np
from netCDF4 import date2num, num2date

VALID_CALENDARS = ('gregorian', 'noleap')
DEFAULT_CALENDAR = 'gregorian'
# tolerance (days) for time levels that land just before Jan. 1
_TOL_DAYS = 1.0e-6


def get_calendar(attrs):
    """
    Return the calendar from a mapping/object of global attributes.

    Parameters
    ----------
    attrs : dict-like or netCDF4.Dataset
        Source of the config_calendar_type global attribute.
    """
    if hasattr(attrs, 'ncattrs'):
        attrs = {k: attrs.getncattr(k) for k in attrs.ncattrs()}
    calendar = attrs.get('config_calendar_type', DEFAULT_CALENDAR)
    if calendar not in VALID_CALENDARS:
        raise ValueError(
            f"config_calendar_type must be one of {VALID_CALENDARS}, "
            f"got '{calendar}'")
    return calendar


def year_start_days(years, start_date, calendar):
    """
    Days since start_date ('YYYY-MM-DD') of Jan. 1 of each year in years.
    """
    units = f'days since {start_date}'
    years = np.atleast_1d(years)
    return np.array([
        date2num(cftime.datetime(int(y), 1, 1, calendar=calendar),
                 units=units, calendar=calendar)
        for y in years
    ], dtype=float)


def decimal_years(days_since_start, start_date, calendar):
    """
    Decimal year for each daysSinceStart value, equal to an integer exactly on
    Jan. 1, with the fraction of the year elapsed (using the real length of
    each year in the given calendar) otherwise.
    """
    days = np.atleast_1d(np.asarray(days_since_start, dtype=float))
    ref_year = int(start_date[0:4])
    n_years = int(np.ceil(np.nanmax(days) / 365.0)) + 2
    bounds = year_start_days(ref_year + np.arange(n_years + 1),
                             start_date, calendar)
    idx = np.searchsorted(bounds, days + _TOL_DAYS, side='right') - 1
    frac = (days - bounds[idx]) / (bounds[idx + 1] - bounds[idx])
    return ref_year + idx + frac


def to_days_since_1850(days, start_date, calendar):
    """
    Convert days since start_date in the given calendar to days since
    1850-01-01 in the standard (Gregorian) calendar.
    """
    return date2num(
        num2date(days, units=f'days since {start_date}', calendar=calendar),
        units='days since 1850-01-01',
        calendar='standard',
    )
