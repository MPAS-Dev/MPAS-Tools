"""
Calendar-aware time helpers for the ISMIP7 post-processing scripts.

MALI's ``daysSinceStart`` counts days of the model calendar
(``config_calendar_type`` = 'gregorian' or 'noleap') since
``simulationStartTime``.  The ismip7_run compass test cases use 'gregorian',
so years are not all 365 days long and conversions must use the calendar of
the run rather than assuming 365-day years.
"""

import numpy as np
from netCDF4 import num2date, date2num
import cftime

VALID_CALENDARS = ('noleap', 'gregorian')


def get_calendar(attrs, default='gregorian'):
    """
    Return the model calendar from the global attributes of a MALI output file
    (the ``config_calendar_type`` attribute).  ``attrs`` is a mapping such as
    ``xarray.Dataset.attrs``.  Falls back to ``default`` if the attribute is
    missing.
    """
    calendar = attrs.get('config_calendar_type', default)
    if calendar not in VALID_CALENDARS:
        raise ValueError(
            f"config_calendar_type must be one of {VALID_CALENDARS}, "
            f"got '{calendar}'")
    return calendar


def _cftime_calendar(calendar):
    # netCDF4/cftime name for the model calendar
    return 'noleap' if calendar == 'noleap' else 'proleptic_gregorian'


def decimal_years(days_since_start, sim_start_date, calendar):
    """
    Convert days since ``sim_start_date`` ('YYYY-MM-DD') to decimal years,
    where an integer value means Jan. 1 of that year at 00:00 and the
    fractional part is the fraction of that calendar year elapsed.
    """
    days_since_start = np.atleast_1d(np.asarray(days_since_start, dtype=float))
    cal = _cftime_calendar(calendar)
    dates = num2date(days_since_start, units=f'days since {sim_start_date}',
                     calendar=cal)
    out = np.empty(days_since_start.shape)
    for i, date in enumerate(dates):
        year = date.year
        start = cftime.datetime(year, 1, 1, calendar=cal)
        end = cftime.datetime(year + 1, 1, 1, calendar=cal)
        out[i] = year + (date - start).total_seconds() / \
            (end - start).total_seconds()
    return out


def days_at_jan1(year, sim_start_date, calendar):
    """
    Days since ``sim_start_date`` of Jan. 1 00:00 of ``year`` in the model
    calendar.  ``year`` may be a scalar or array of integers.
    """
    cal = _cftime_calendar(calendar)
    units = f'days since {sim_start_date}'
    year = np.atleast_1d(np.asarray(year))
    dates = [cftime.datetime(int(y), 1, 1, calendar=cal) for y in year]
    return np.asarray(date2num(dates, units=units, calendar=cal), dtype=float)


def days_to_ismip7_time(days_since_start, sim_start_date, calendar):
    """
    Convert days since ``sim_start_date`` on the model calendar to days since
    1850-01-01 on the standard (Gregorian) calendar, the time axis used by
    ISMIP7 output files.
    """
    cal = _cftime_calendar(calendar)
    return date2num(
        num2date(days_since_start, units=f'days since {sim_start_date}',
                 calendar=cal),
        units='days since 1850-01-01',
        calendar='standard',
    )
