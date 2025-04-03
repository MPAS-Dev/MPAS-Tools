from datetime import datetime
from typing import Optional

import numpy as np
import xarray as xr


def _set_time(ds: xr.Dataset, xtime_var: Optional[str]):
    """
    Set the time variable in the dataset to be a DataArray of seconds since
    the first entry.

    Parameters
    ----------
    ds : xr.Dataset
        The dataset containing the time variable.
    xtime_var : str
        The name of the time variable in the dataset.
    """
    if 'Time' in ds.dims:
        if xtime_var is not None:
            if xtime_var not in ds.data_vars:
                raise ValueError(
                    f"xtime variable '{xtime_var}' not found in dataset."
                )
            ds['Time'] = _xtime_to_seconds(ds[xtime_var])
        else:
            ds['Time'] = xr.DataArray(np.arange(ds.sizes['Time']), dims='Time')


def _xtime_to_seconds(xtime: xr.DataArray) -> xr.DataArray:
    """
    Convert an xarray DataArray of xtime strings to seconds since the first
    entry.

    Parameters
    ----------
    xtime : xr.DataArray
        An array of strings representing time in the format
        'YYYY-MM-DD_HH:MM:SS.sss'.

    Returns
    -------
    xr.DataArray
        An array of seconds since the first entry in `xtime`.
    """
    # Convert xtime strings to datetime objects using datetime.strptime
    timestamps = [
        datetime.strptime(time_str, '%Y-%m-%d_%H:%M:%S.%f')
        for time_str in xtime.values.astype(str)
    ]
    # Calculate seconds since the first timestamp
    seconds_since_start = [
        (ts - timestamps[0]).total_seconds() for ts in timestamps
    ]
    # Return as a DataArray
    return xr.DataArray(
        seconds_since_start, dims=xtime.dims, coords=xtime.coords
    )
