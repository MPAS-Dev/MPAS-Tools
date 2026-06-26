"""
Generic validation utilities for sets of MALI output files.
"""

import os
import xarray as xr


def validate_mali_files(files, required_vars, label=''):
    """
    Validate a sorted list of MALI output files before processing.

    Checks that:
    - The list is not empty
    - Each file exists
    - Each file contains all required variables
    - simulationStartTime is consistent across all files
    - No time overlaps (daysSinceStart) exist between consecutive files
    - No unexpectedly large time gaps (> 366 days) exist between consecutive files

    Parameters
    ----------
    files : list of str
        Sorted list of file paths.
    required_vars : list of str
        Variable names that must be present in every file.
    label : str, optional
        Human-readable label for the file type, used in messages.

    Raises
    ------
    ValueError
        If any validation check fails.
    FileNotFoundError
        If any file does not exist.
    """
    tag = f' ({label})' if label else ''

    if len(files) == 0:
        raise ValueError(f"No files provided{tag}.")

    for f in files:
        if not os.path.exists(f):
            raise FileNotFoundError(f"File not found{tag}: {f}")

    # Check required variables in each file
    for f in files:
        with xr.open_dataset(f, decode_cf=False) as ds:
            missing = [v for v in required_vars if v not in ds]
        if missing:
            raise ValueError(
                f"File '{f}'{tag} is missing expected variables: {missing}")

    # Check simulationStartTime consistency across all files
    if len(files) > 1:
        start_times = []
        for f in files:
            with xr.open_dataset(f, decode_cf=False) as ds:
                start_times.append(
                    ds['simulationStartTime'].values.tobytes()
                    .decode('utf-8').strip().strip('\x00'))
        if len(set(start_times)) > 1:
            raise ValueError(
                f"Inconsistent simulationStartTime across files{tag}: "
                f"{set(start_times)}")

    # Check for time overlaps or large gaps between consecutive files
    for i in range(len(files) - 1):
        with xr.open_dataset(files[i], decode_cf=False) as ds_a:
            end_a = float(ds_a['daysSinceStart'].values[-1])
        with xr.open_dataset(files[i + 1], decode_cf=False) as ds_b:
            start_b = float(ds_b['daysSinceStart'].values[0])
        if start_b <= end_a:
            raise ValueError(
                f"Time overlap detected between files{tag}:\n"
                f"  {files[i]} (ends at day {end_a})\n"
                f"  {files[i + 1]} (starts at day {start_b})")
        gap_days = start_b - end_a
        if gap_days > 366:
            print(f"WARNING: Gap of {gap_days:.1f} days between files{tag}:\n"
                  f"  {files[i]}\n  {files[i + 1]}")

    print(f"Validated {len(files)} file(s){tag}.")
