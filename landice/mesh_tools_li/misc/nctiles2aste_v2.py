from __future__ import annotations

from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np
import xarray as xr


def _get_struct_value(struct: Any, key: str) -> Any:
    """Get a field from dict-like or attribute-style MATLAB-loaded structs."""
    if isinstance(struct, dict):
        if key not in struct:
            raise KeyError(f"Missing key '{key}' in tiles structure.")
        return struct[key]

    if hasattr(struct, key):
        return getattr(struct, key)

    dtype = getattr(struct, "dtype", None)
    if dtype is not None and getattr(dtype, "names", None) and key in dtype.names:
        return struct[key]

    raise KeyError(f"Missing key '{key}' in tiles structure.")


def _to_numpy(value: Any) -> np.ndarray:
    """Convert nested scalar/object values (common in MATLAB imports) to ndarray."""
    arr = np.asarray(value)
    while arr.dtype == object and arr.size == 1:
        arr = np.asarray(arr.item())
    return arr


def _to_int_tuple_2(value: Any, name: str) -> tuple[int, int]:
    arr = _to_numpy(value).squeeze()
    if arr.size < 2:
        raise ValueError(f"{name} must contain at least 2 values.")
    return int(arr.flat[0]), int(arr.flat[1])


def _to_int(value: Any, name: str) -> int:
    arr = _to_numpy(value).squeeze()
    if arr.size != 1:
        raise ValueError(f"{name} must be scalar.")
    return int(arr.item())


def _matlab_cell_to_list(value: Any) -> list[Any]:
    """Convert MATLAB cell-like arrays to a Python list."""
    arr = _to_numpy(value)
    if arr.dtype == object:
        return [arr.flat[idx] for idx in range(arr.size)]
    return [arr]


def _get_tile_indices(tiles_index: Any, tile_num: int) -> np.ndarray:
    """Get 0-based linear indices for a 1-based tile number."""
    index_list = _matlab_cell_to_list(tiles_index)
    if tile_num < 1 or tile_num > len(index_list):
        raise ValueError(
            f"Requested tile number {tile_num} is out of range 1..{len(index_list)}"
        )

    idx = _to_numpy(index_list[tile_num - 1]).astype(np.int64).ravel()
    # ASTE MATLAB indices are 1-based linear indices in column-major order.
    return idx - 1


def tiles2field(
    field_tiles: Sequence[np.ndarray] | np.ndarray,
    indices: Sequence[np.ndarray] | np.ndarray,
    field_size: Sequence[int],
) -> np.ndarray:
    """
    Reconstruct a compact ASTE field from tile arrays.

    Parameters
    ----------
    field_tiles
        A list of 2D/3D tile arrays (or a single tile array).
    indices
        A list of linear indices (or one index array), one per tile.
        Indices are expected to be 0-based and follow MATLAB column-major
        flattening semantics.
    field_size
        Final compact field size ``(nx, ny[, nz])``.

    Returns
    -------
    np.ndarray
        Compact ASTE field of shape ``(nx, ny, nz)`` with zeros in locations
        not covered by provided tiles.
    """
    field_size_arr = np.asarray(field_size, dtype=int).ravel()
    if field_size_arr.size < 2:
        raise ValueError("field_size must contain at least nx and ny.")

    dfieldx = int(field_size_arr[0])
    dfieldy = int(field_size_arr[1])

    tile_list = list(field_tiles) if isinstance(field_tiles, (list, tuple)) else [field_tiles]
    index_list = list(indices) if isinstance(indices, (list, tuple)) else [indices]

    if len(tile_list) != len(index_list):
        raise ValueError("There must be one set of indices for each tile.")

    first_tile = np.asarray(tile_list[0])
    if first_tile.ndim == 2:
        dfieldz = 1
        tile_shape_xy = first_tile.shape
    elif first_tile.ndim == 3:
        dfieldz = first_tile.shape[2]
        tile_shape_xy = first_tile.shape[:2]
    else:
        raise ValueError("Tiles must be 2D or 3D arrays.")

    field = np.zeros((dfieldx, dfieldy, dfieldz), dtype=float)

    for tile, idx in zip(tile_list, index_list):
        tile_arr = np.asarray(tile)
        idx_arr = np.asarray(idx, dtype=np.int64).ravel()

        if tile_arr.ndim == 2:
            tile_arr = tile_arr[:, :, np.newaxis]

        if tile_arr.shape[:2] != tile_shape_xy:
            raise ValueError("All tiles must have the same x/y dimensions.")

        if tile_arr.shape[0] * tile_arr.shape[1] != idx_arr.size:
            raise ValueError(
                "Number of points in indices does not match tile x*y size."
            )

        for iz in range(dfieldz):
            tmpf = field[:, :, iz]
            flat_tmpf = tmpf.reshape(-1, order="F")
            tile_flat = tile_arr[:, :, iz].reshape(-1, order="F")
            flat_tmpf[idx_arr] = tile_flat
            field[:, :, iz] = flat_tmpf.reshape((dfieldx, dfieldy), order="F")

    return field


def nctiles2aste_v2(
    var_name: str,
    root_dir: str | Path,
    tile_list: Iterable[int],
    levels_list: Iterable[int],
    times_list: Iterable[int],
    tiles: Any,
    flag_aste: bool,
) -> tuple[np.ndarray | list[np.ndarray], np.ndarray]:
    """
    Read tile NetCDF files and return either per-tile arrays or compact ASTE fields.

    Notes
    -----
    ``tile_list``, ``levels_list``, and ``times_list`` are interpreted as
    1-based indices to match MATLAB behavior.
    """
    root_dir = Path(root_dir)
    data_dir = root_dir / var_name
    if not data_dir.is_dir():
        raise FileNotFoundError(f"Data directory not found: {data_dir}")

    file_list = sorted(data_dir.glob("*.nc"))
    if not file_list:
        raise FileNotFoundError(f"No NetCDF files found in: {data_dir}")

    tile_list = [int(tile) for tile in tile_list]
    levels_list = [int(level) for level in levels_list]
    times_list = [int(time) for time in times_list]

    nx, ny = _to_int_tuple_2(_get_struct_value(tiles, "fsize"), "tiles.fsize")
    dtilex, dtiley = _to_int_tuple_2(_get_struct_value(tiles, "tsize"), "tiles.tsize")
    wet_tiles = _to_int(_get_struct_value(tiles, "wet"), "tiles.wet")
    tiles_index = _get_struct_value(tiles, "index")

    sample_ds = xr.open_dataset(file_list[0])
    try:
        if var_name not in sample_ds:
            raise KeyError(f"Variable '{var_name}' not found in {file_list[0].name}")
        if "timstep" not in sample_ds:
            raise KeyError(f"Variable 'timstep' not found in {file_list[0].name}")

        sample_var = sample_ds[var_name]
        if sample_var.ndim == 3:
            nz = 1
            ntimes_max = int(sample_var.shape[2])
            file_ndims = 4
        elif sample_var.ndim == 4:
            nz = int(sample_var.shape[2])
            ntimes_max = int(sample_var.shape[3])
            file_ndims = 5
        else:
            raise ValueError(
                f"Unsupported variable rank for '{var_name}': {sample_var.ndim}. "
                "Expected 3D (x, y, time) or 4D (x, y, z, time)."
            )

        if dtilex != int(sample_var.shape[0]) or dtiley != int(sample_var.shape[1]):
            raise ValueError(
                "nctiles2aste ERROR: Dimensions from tiles.tsize "
                f"{dtilex}x{dtiley} do not match NetCDF variable dimensions "
                f"{sample_var.shape[0]}x{sample_var.shape[1]}."
            )
    finally:
        sample_ds.close()

    ntiles_max = len(file_list)
    ntiles = len(tile_list)
    nlevels = len(levels_list)
    ntimes = len(times_list)

    if ntiles_max != wet_tiles:
        print(
            "nctiles2aste WARNING: number of tiles in directory "
            f"{ntiles_max} is different from number of wet tiles {wet_tiles}"
        )

    if ntiles > wet_tiles:
        raise ValueError(
            "nctiles2aste ERROR: number of requested tiles exceeds wet tile count "
            f"({ntiles} > {wet_tiles})."
        )

    if max(tile_list) > wet_tiles:
        raise ValueError(
            "nctiles2aste ERROR: requested tile number exceeds wet tile count "
            f"({max(tile_list)} > {wet_tiles})."
        )

    if nlevels > nz:
        raise ValueError(
            "nctiles2aste ERROR: number of requested levels exceeds available "
            f"levels ({nlevels} > {nz})."
        )

    if max(levels_list) > nz:
        raise ValueError(
            "nctiles2aste ERROR: requested level exceeds available levels "
            f"({max(levels_list)} > {nz})."
        )

    if max(times_list) > ntimes_max:
        raise ValueError(
            "nctiles2aste ERROR: requested time step exceeds available time steps "
            f"({max(times_list)} > {ntimes_max})."
        )

    field_tiles: list[np.ndarray] = []
    indices: list[np.ndarray] = []
    time_step = np.full((ntimes,), np.nan)

    for itile, tile_num in enumerate(tile_list):
        tile_indices = _get_tile_indices(tiles_index, tile_num)
        indices.append(tile_indices)

        tile_fname = data_dir / f"{var_name}.{tile_num:04d}.nc"
        if not tile_fname.is_file():
            raise FileNotFoundError(f"Tile file not found: {tile_fname}")

        tmp = np.full((dtilex, dtiley, nlevels, ntimes), np.nan)

        ds = xr.open_dataset(tile_fname)
        try:
            var = ds[var_name]
            if itile == 0:
                for itime, time_1based in enumerate(times_list):
                    time_step[itime] = (
                        ds["timstep"].isel({ds["timstep"].dims[0]: time_1based - 1}).item()
                    )

            if file_ndims == 4:
                x_dim, y_dim, t_dim = var.dims
                for itime, time_1based in enumerate(times_list):
                    slice2d = var.isel({t_dim: time_1based - 1}).values
                    tmp[:, :, 0, itime] = np.asarray(slice2d)
            else:
                x_dim, y_dim, z_dim, t_dim = var.dims
                for itime, time_1based in enumerate(times_list):
                    for ilev, level_1based in enumerate(levels_list):
                        slice2d = var.isel(
                            {z_dim: level_1based - 1, t_dim: time_1based - 1}
                        ).values
                        tmp[:, :, ilev, itime] = np.asarray(slice2d)

            field_tiles.append(tmp)
        finally:
            ds.close()

    if flag_aste:
        print(
            "WARNING: putting fields into full ASTE compact format; "
            "this can be memory intensive"
        )
        field = np.full((nx, ny, nlevels, ntimes), np.nan)
        for itime in range(ntimes):
            time_tiles = [tile[:, :, :, itime] for tile in field_tiles]
            field[:, :, :, itime] = tiles2field(time_tiles, indices, (nx, ny, nlevels))
    else:
        field = field_tiles

    return field, time_step
