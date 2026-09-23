#!/usr/bin/env python3
"""Estimate grounding-line ice discharge from a MALI/MPAS-Landice mesh.

The script can use either observed surface velocity as depth-uniform plug flow,
or the modeled 3-D uReconstructX/Y fields.  Modeled interface velocities are
averaged into layers and integrated with layerThicknessFractions.

For every interior edge separating an ice-covered grounded cell from an
ice-covered floating cell, the script evaluates

    q_e = rho_i * H_e * (u_e dot n_grounded_to_floating) * dvEdge_e

using arithmetic cell-to-edge averages for H and u.  Positive flux is from
grounded ice into the floating shelf.

Geometry variables:
    thickness, bedTopography, cellsOnEdge, dvEdge, xCell, yCell

Velocity variables:
    observed: observedSurfaceVelocityX, observedSurfaceVelocityY
    modeled:  uReconstructX, uReconstructY, layerThicknessFractions

Optional:
    observedSurfaceVelocityUncertainty
    a regions file containing regionEdgeMasks, regionCellMasks, and regionNames

All MALI velocity fields are assumed to use Registry units of m/s.  Regional
masks are evaluated independently and may overlap, so regional values are not
assumed to sum to the global value.

Examples:
    python estimate_gl_flux.py -f landice_grid.nc
    python estimate_gl_flux.py -f output.nc --velocity-source modeled
    python estimate_gl_flux.py -f output.nc --velocity-source both
    python estimate_gl_flux.py -f output.nc -r region_masks.nc --velocity-source both
"""

from __future__ import annotations

import argparse
import math
import sys
from dataclasses import dataclass
from typing import Any, Optional

import numpy as np


SECONDS_PER_YEAR = 365.0 * 24.0 * 60.0 * 60.0
KG_PER_GT = 1.0e12
ICE_DENSITY = 910.0
WATER_DENSITY = 1028.0
SEA_LEVEL = 0.0
MIN_THICKNESS = 1.0


@dataclass
class FluxResult:
    n_grounded_cells: int
    n_floating_cells: int
    n_gl_edges: int
    signed_gt_per_year: float
    outward_gt_per_year: float
    inward_gt_per_year: float
    absolute_gt_per_year: float
    uncertainty_gt_per_year: Optional[float]


def _as_array(variable: Any, time_index: int, expected_dim: str) -> np.ndarray:
    """Read a 1-D mesh field, selecting Time if present."""
    dims = tuple(variable.dimensions)
    index = []
    for dim in dims:
        if dim.lower() == "time":
            index.append(time_index)
        elif dim == expected_dim:
            index.append(slice(None))
        else:
            # Permit singleton dimensions, but reject ambiguous multi-level data.
            size = variable.shape[len(index)]
            if size != 1:
                raise ValueError(
                    f"Variable {variable.name!r} has unsupported dimension "
                    f"{dim!r} of size {size}; expected only Time and {expected_dim}."
                )
            index.append(0)

    if expected_dim not in dims:
        raise ValueError(
            f"Variable {variable.name!r} does not have expected dimension "
            f"{expected_dim!r}; dimensions are {dims}."
        )
    data = variable[tuple(index)]
    if np.ma.isMaskedArray(data):
        data = data.filled(np.nan)
    return np.asarray(data, dtype=float).reshape(-1)


def _as_cell_interface_array(variable: Any, time_index: int) -> np.ndarray:
    """Read a cell-by-vertical velocity field, regardless of dimension order."""
    dims = tuple(variable.dimensions)
    if "nCells" not in dims or "nVertInterfaces" not in dims:
        raise ValueError(
            f"Variable {variable.name!r} must contain nCells and "
            f"nVertInterfaces; dimensions are {dims}."
        )

    index = []
    retained_dims = []
    for axis, dim in enumerate(dims):
        if dim.lower() == "time":
            index.append(time_index)
        elif dim in ("nCells", "nVertInterfaces"):
            index.append(slice(None))
            retained_dims.append(dim)
        elif variable.shape[axis] == 1:
            index.append(0)
        else:
            raise ValueError(
                f"Variable {variable.name!r} has unsupported dimension {dim!r} "
                f"of size {variable.shape[axis]}."
            )

    data = variable[tuple(index)]
    if np.ma.isMaskedArray(data):
        data = data.filled(np.nan)
    data = np.asarray(data, dtype=float)
    cell_axis = retained_dims.index("nCells")
    vertical_axis = retained_dims.index("nVertInterfaces")
    data = np.moveaxis(data, (cell_axis, vertical_axis), (0, 1))
    if data.ndim != 2:
        raise ValueError(f"Could not reduce {variable.name!r} to a 2-D field; shape is {data.shape}.")
    return data


def depth_average_velocity(
    velocity: np.ndarray,
    layer_thickness_fractions: np.ndarray,
) -> np.ndarray:
    """Average interface velocities into layers, then return the vertical mean."""
    velocity = np.asarray(velocity, dtype=float)
    fractions = np.asarray(layer_thickness_fractions, dtype=float).reshape(-1)
    if velocity.ndim != 2:
        raise ValueError("Modeled velocity must have shape (nCells, nVertical).")
    if fractions.size == 0 or np.any(~np.isfinite(fractions)) or np.any(fractions < 0.0):
        raise ValueError("layerThicknessFractions must be finite, nonnegative, and nonempty.")
    fraction_sum = float(fractions.sum())
    if fraction_sum <= 0.0:
        raise ValueError("layerThicknessFractions must have a positive sum.")
    if not np.isclose(fraction_sum, 1.0, rtol=1.0e-8, atol=1.0e-12):
        raise ValueError(
            f"layerThicknessFractions sum to {fraction_sum:.16g}, not 1.0."
        )

    if velocity.shape[1] != fractions.size + 1:
        raise ValueError(
            f"Modeled velocity has {velocity.shape[1]} interfaces but "
            f"layerThicknessFractions has {fractions.size} layers."
        )
    layer_velocity = 0.5 * (velocity[:, :-1] + velocity[:, 1:])

    return layer_velocity @ fractions


def _read_cells_on_edge(variable: Any) -> np.ndarray:
    """Read cellsOnEdge and return shape (nEdges, 2), retaining 1-based IDs."""
    dims = tuple(variable.dimensions)
    data = variable[:]
    if np.ma.isMaskedArray(data):
        data = data.filled(0)
    data = np.asarray(data, dtype=np.int64)

    if data.ndim != 2:
        raise ValueError(f"{variable.name!r} must be 2-D; got shape {data.shape}.")
    if data.shape[1] == 2:
        return data
    if data.shape[0] == 2:
        return data.T
    raise ValueError(
        f"Could not identify TWO dimension in {variable.name!r}; shape is "
        f"{data.shape}, dimensions are {dims}."
    )


def _as_region_masks(variable: Any, mesh_dim: str) -> np.ndarray:
    """Read a mesh-by-region integer mask, regardless of dimension order."""
    dims = tuple(variable.dimensions)
    if mesh_dim not in dims or "nRegions" not in dims:
        raise ValueError(
            f"Variable {variable.name!r} must contain {mesh_dim} and nRegions; "
            f"dimensions are {dims}."
        )
    data = variable[:]
    if np.ma.isMaskedArray(data):
        data = data.filled(0)
    data = np.asarray(data)
    mesh_axis = dims.index(mesh_dim)
    region_axis = dims.index("nRegions")
    data = np.moveaxis(data, (mesh_axis, region_axis), (0, 1))
    if data.ndim != 2:
        raise ValueError(f"Variable {variable.name!r} must be 2-D; shape is {data.shape}.")
    return data != 0


def _read_region_names(variable: Any, n_regions: int) -> list[str]:
    """Decode the Registry-style char regionNames(nRegions, StrLen) array."""
    dims = tuple(variable.dimensions)
    if "nRegions" not in dims:
        raise ValueError(
            f"Variable {variable.name!r} must contain nRegions; dimensions are {dims}."
        )
    data = variable[:]
    if np.ma.isMaskedArray(data):
        fill = b" " if data.dtype.kind == "S" else " "
        data = data.filled(fill)
    data = np.asarray(data)
    data = np.moveaxis(data, dims.index("nRegions"), 0)
    if data.shape[0] != n_regions:
        raise ValueError(
            f"regionNames has {data.shape[0]} entries but masks have {n_regions} regions."
        )

    names = []
    for index, row in enumerate(data):
        chars = np.asarray(row).reshape(-1)
        if chars.dtype.kind == "S":
            name = b"".join(chars.tolist()).decode("utf-8", errors="replace")
        else:
            name = "".join(str(char) for char in chars.tolist())
        name = name.split("\x00", 1)[0].strip()
        names.append(name or f"Region {index + 1}")
    return names


def compute_flux(
    thickness: np.ndarray,
    bed: np.ndarray,
    cells_on_edge_one_based: np.ndarray,
    edge_length: np.ndarray,
    x_cell: np.ndarray,
    y_cell: np.ndarray,
    velocity_x_m_per_year: np.ndarray,
    velocity_y_m_per_year: np.ndarray,
    velocity_uncertainty_m_per_year: Optional[np.ndarray] = None,
    edge_selector: Optional[np.ndarray] = None,
    cell_selector: Optional[np.ndarray] = None,
) -> FluxResult:
    """Compute grounding-line discharge for a velocity field.

    Edge-normal directions are calculated from the vector joining the two
    cell centers, which is normal to their shared edge on a planar MPAS Voronoi
    mesh.  The uncertainty calculation treats the supplied value as an
    independent, isotropic 1-sigma uncertainty for each cell's x and y
    velocity components.  Shared cells are accounted for.
    """
    thickness = np.asarray(thickness, dtype=float)
    bed = np.asarray(bed, dtype=float)
    ux = np.asarray(velocity_x_m_per_year, dtype=float)
    uy = np.asarray(velocity_y_m_per_year, dtype=float)
    x_cell = np.asarray(x_cell, dtype=float)
    y_cell = np.asarray(y_cell, dtype=float)
    cells = np.asarray(cells_on_edge_one_based, dtype=np.int64)
    edge_length = np.asarray(edge_length, dtype=float)

    n_cells = thickness.size
    if any(a.shape != (n_cells,) for a in (bed, ux, uy, x_cell, y_cell)):
        raise ValueError(
            "thickness, bed, xCell, yCell, velocity X, and velocity Y must "
            "have the same 1-D shape."
        )
    n_edges = cells.shape[0]
    if cells.shape != (n_edges, 2):
        raise ValueError("cellsOnEdge must have shape (nEdges, 2).")
    if edge_length.shape != (n_edges,):
        raise ValueError("dvEdge must have length nEdges.")
    if edge_selector is None:
        edge_selector = np.ones(n_edges, dtype=bool)
    else:
        edge_selector = np.asarray(edge_selector, dtype=bool)
        if edge_selector.shape != (n_edges,):
            raise ValueError("Regional edge mask must have length nEdges.")
    if cell_selector is None:
        cell_selector = np.ones(n_cells, dtype=bool)
    else:
        cell_selector = np.asarray(cell_selector, dtype=bool)
        if cell_selector.shape != (n_cells,):
            raise ValueError("Regional cell mask must have length nCells.")

    finite_geometry = np.isfinite(thickness) & np.isfinite(bed)
    ice = finite_geometry & (thickness > MIN_THICKNESS)
    water_depth = np.maximum(SEA_LEVEL - bed, 0.0)
    # Positive flotation residual means ice overburden exceeds ocean pressure.
    flotation_residual = ICE_DENSITY * thickness - WATER_DENSITY * water_depth
    grounded = ice & (flotation_residual > 0.0)
    floating = ice & ~grounded

    # MPAS connectivity is 1-based; zero denotes no neighboring cell.
    interior = (cells[:, 0] > 0) & (cells[:, 1] > 0)
    c0_raw = cells[:, 0] - 1
    c1_raw = cells[:, 1] - 1
    valid_ids = (
        (c0_raw >= 0) & (c0_raw < n_cells)
        & (c1_raw >= 0) & (c1_raw < n_cells)
    )
    # Clipping makes indexing safe before the interior/valid-ID mask is applied.
    c0 = np.clip(c0_raw, 0, max(n_cells - 1, 0))
    c1 = np.clip(c1_raw, 0, max(n_cells - 1, 0))

    gf = interior & valid_ids & grounded[c0] & floating[c1]
    fg = interior & valid_ids & floating[c0] & grounded[c1]
    gl = (gf | fg) & edge_selector
    edge_ids = np.flatnonzero(gl)

    if edge_ids.size == 0:
        return FluxResult(
            int((grounded & cell_selector).sum()),
            int((floating & cell_selector).sum()), 0,
            0.0, 0.0, 0.0, 0.0,
            0.0 if velocity_uncertainty_m_per_year is not None else None,
        )

    i = c0[edge_ids]
    j = c1[edge_ids]
    # The vector from c0 to c1 is normal to their shared Voronoi edge. Reverse
    # it where c1 is grounded so every normal points grounded -> floating.
    orientation = np.where(gf[edge_ids], 1.0, -1.0)
    dx = x_cell[j] - x_cell[i]
    dy = y_cell[j] - y_cell[i]
    center_distance = np.hypot(dx, dy)
    if np.any(~np.isfinite(center_distance)) or np.any(center_distance <= 0.0):
        bad = int((~np.isfinite(center_distance) | (center_distance <= 0.0)).sum())
        raise ValueError(f"Found {bad} grounding-line edges with invalid cell-center geometry.")
    nx = orientation * dx / center_distance
    ny = orientation * dy / center_distance

    h_edge = 0.5 * (thickness[i] + thickness[j])
    ux_edge = 0.5 * (ux[i] + ux[j])
    uy_edge = 0.5 * (uy[i] + uy[j])
    normal_velocity = ux_edge * nx + uy_edge * ny
    mass_flux = ICE_DENSITY * h_edge * normal_velocity * edge_length[edge_ids]  # kg/yr

    finite_edge = (
        np.isfinite(mass_flux)
        & np.isfinite(edge_length[edge_ids])
        & (edge_length[edge_ids] > 0.0)
    )
    if not np.all(finite_edge):
        bad = int((~finite_edge).sum())
        raise ValueError(f"Found {bad} grounding-line edges with invalid flux or nonpositive dvEdge.")

    signed = float(mass_flux.sum())
    outward = float(mass_flux[mass_flux > 0.0].sum())
    inward = float(mass_flux[mass_flux < 0.0].sum())
    absolute = float(np.abs(mass_flux).sum())

    uncertainty_gt_per_year = None
    if velocity_uncertainty_m_per_year is not None:
        sigma = np.asarray(velocity_uncertainty_m_per_year, dtype=float)
        if sigma.shape != (n_cells,):
            raise ValueError("Velocity uncertainty must have the same shape as thickness.")
        gl_cells = np.unique(np.concatenate((i, j)))
        if np.any(~np.isfinite(sigma[gl_cells])) or np.any(sigma[gl_cells] < 0.0):
            raise ValueError("Velocity uncertainty is missing, nonfinite, or negative on grounding-line cells.")

        # Q is linear in cell velocity. Aggregate coefficients first so that
        # correlations caused by one cell participating in several edges are
        # represented exactly under the independent-cell assumption.
        k = ICE_DENSITY * h_edge * edge_length[edge_ids]
        coeff_x = np.zeros(n_cells, dtype=float)
        coeff_y = np.zeros(n_cells, dtype=float)
        np.add.at(coeff_x, i, 0.5 * k * nx)
        np.add.at(coeff_x, j, 0.5 * k * nx)
        np.add.at(coeff_y, i, 0.5 * k * ny)
        np.add.at(coeff_y, j, 0.5 * k * ny)
        variance = np.sum((sigma * coeff_x) ** 2 + (sigma * coeff_y) ** 2)
        uncertainty_gt_per_year = math.sqrt(float(variance)) / KG_PER_GT

    return FluxResult(
        n_grounded_cells=int((grounded & cell_selector).sum()),
        n_floating_cells=int((floating & cell_selector).sum()),
        n_gl_edges=int(edge_ids.size),
        signed_gt_per_year=signed / KG_PER_GT,
        outward_gt_per_year=outward / KG_PER_GT,
        inward_gt_per_year=inward / KG_PER_GT,
        absolute_gt_per_year=absolute / KG_PER_GT,
        uncertainty_gt_per_year=uncertainty_gt_per_year,
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Estimate observed and/or modeled grounding-line discharge on a MALI mesh.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-f", "--file", required=True, help="MALI NetCDF input file")
    parser.add_argument(
        "-r", "--regions-file",
        help="Optional NetCDF file containing regionEdgeMasks, regionCellMasks, and regionNames",
    )
    parser.add_argument("--time-index", type=int, default=-1, help="Time record (Python indexing)")
    parser.add_argument(
        "--velocity-source", choices=("observed", "modeled", "both"), default="observed",
        help="Velocity field(s) used for the flux calculation",
    )
    return parser


def _require_variable(dataset: Any, name: str) -> Any:
    if name not in dataset.variables:
        available = ", ".join(sorted(dataset.variables))
        raise KeyError(f"Required variable {name!r} not found. Available variables: {available}")
    return dataset.variables[name]


def _print_result(label: str, result: FluxResult) -> None:
    print(f"\n{label} velocity")
    print("  Input velocity units:         m/s (converted to m/yr using 365 days/yr)")
    print(f"  Grounded ice cells:           {result.n_grounded_cells}")
    print(f"  Floating ice cells:           {result.n_floating_cells}")
    print(f"  Grounding-line edges:         {result.n_gl_edges}")
    print(f"  Signed GL discharge:          {result.signed_gt_per_year:.6g} Gt/yr")
    print(f"  Outward-only contribution:    {result.outward_gt_per_year:.6g} Gt/yr")
    print(f"  Inward-only contribution:     {result.inward_gt_per_year:.6g} Gt/yr")
    print(f"  Sum of absolute contributions:{result.absolute_gt_per_year: .6g} Gt/yr")
    if result.uncertainty_gt_per_year is not None:
        print(f"  Approx. 1-sigma uncertainty:  {result.uncertainty_gt_per_year:.6g} Gt/yr")
        print("    (Assumes independent cells and isotropic component uncertainties.)")


def _print_regional_results(
    label: str,
    names: list[str],
    results: list[FluxResult],
) -> None:
    """Print compact regional grounding-line metrics."""
    include_uncertainty = any(result.uncertainty_gt_per_year is not None for result in results)
    print(f"\nRegional metrics: {label} velocity")
    header = (
        f"  {'Region':<28} {'Grounded':>10} {'Floating':>10} {'GL edges':>9} "
        f"{'Signed':>12} {'Outward':>12} {'Inward':>12} {'Absolute':>12}"
    )
    if include_uncertainty:
        header += f" {'1-sigma':>12}"
    print(header)
    print(f"  {'':<28} {'cells':>10} {'cells':>10} {'':>9} "
          f"{'Gt/yr':>12} {'Gt/yr':>12} {'Gt/yr':>12} {'Gt/yr':>12}"
          + (f" {'Gt/yr':>12}" if include_uncertainty else ""))
    for name, result in zip(names, results):
        row = (
            f"  {name[:28]:<28} {result.n_grounded_cells:10d} "
            f"{result.n_floating_cells:10d} {result.n_gl_edges:9d} "
            f"{result.signed_gt_per_year:12.5g} {result.outward_gt_per_year:12.5g} "
            f"{result.inward_gt_per_year:12.5g} {result.absolute_gt_per_year:12.5g}"
        )
        if include_uncertainty:
            sigma = result.uncertainty_gt_per_year
            row += f" {sigma:12.5g}" if sigma is not None else f" {'--':>12}"
        print(row)


def main(argv: Optional[list[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        from netCDF4 import Dataset
    except ImportError:
        print(
            "ERROR: This script requires netCDF4 (for example, `conda install netcdf4` "
            "or `python -m pip install netCDF4`).",
            file=sys.stderr,
        )
        return 2

    try:
        with Dataset(args.file, "r") as ds:
            thickness = _as_array(_require_variable(ds, "thickness"), args.time_index, "nCells")
            bed = _as_array(_require_variable(ds, "bedTopography"), args.time_index, "nCells")
            cells_on_edge = _read_cells_on_edge(_require_variable(ds, "cellsOnEdge"))
            edge_length = _as_array(
                _require_variable(ds, "dvEdge"), args.time_index, "nEdges"
            )
            x_cell = _as_array(_require_variable(ds, "xCell"), args.time_index, "nCells")
            y_cell = _as_array(_require_variable(ds, "yCell"), args.time_index, "nCells")

            results: dict[str, FluxResult] = {}
            regional_results: dict[str, list[FluxResult]] = {}

            region_names: list[str] = []
            region_edge_masks: Optional[np.ndarray] = None
            region_cell_masks: Optional[np.ndarray] = None
            if args.regions_file is not None:
                with Dataset(args.regions_file, "r") as region_ds:
                    region_edge_masks = _as_region_masks(
                        _require_variable(region_ds, "regionEdgeMasks"), "nEdges"
                    )
                    region_cell_masks = _as_region_masks(
                        _require_variable(region_ds, "regionCellMasks"), "nCells"
                    )
                    if region_edge_masks.shape[0] != edge_length.size:
                        raise ValueError(
                            f"Region file has {region_edge_masks.shape[0]} edges, but "
                            f"the MALI file has {edge_length.size}."
                        )
                    if region_cell_masks.shape[0] != thickness.size:
                        raise ValueError(
                            f"Region file has {region_cell_masks.shape[0]} cells, but "
                            f"the MALI file has {thickness.size}."
                        )
                    if region_edge_masks.shape[1] != region_cell_masks.shape[1]:
                        raise ValueError(
                            "regionEdgeMasks and regionCellMasks contain different "
                            "numbers of regions."
                        )
                    region_names = _read_region_names(
                        _require_variable(region_ds, "regionNames"),
                        region_edge_masks.shape[1],
                    )

            def calculate_metrics(
                ux: np.ndarray,
                uy: np.ndarray,
                uncertainty: Optional[np.ndarray] = None,
            ) -> tuple[FluxResult, list[FluxResult]]:
                common = dict(
                    thickness=thickness,
                    bed=bed,
                    cells_on_edge_one_based=cells_on_edge,
                    edge_length=edge_length,
                    x_cell=x_cell,
                    y_cell=y_cell,
                    velocity_x_m_per_year=ux,
                    velocity_y_m_per_year=uy,
                    velocity_uncertainty_m_per_year=uncertainty,
                )
                global_result = compute_flux(**common)
                region_values = []
                if region_edge_masks is not None and region_cell_masks is not None:
                    for region_index in range(region_edge_masks.shape[1]):
                        region_values.append(
                            compute_flux(
                                **common,
                                edge_selector=region_edge_masks[:, region_index],
                                cell_selector=region_cell_masks[:, region_index],
                            )
                        )
                return global_result, region_values

            if args.velocity_source in ("observed", "both"):
                observed_ux_var = _require_variable(ds, "observedSurfaceVelocityX")
                observed_uy_var = _require_variable(ds, "observedSurfaceVelocityY")
                observed_ux = _as_array(observed_ux_var, args.time_index, "nCells")
                observed_uy = _as_array(observed_uy_var, args.time_index, "nCells")
                observed_ux *= SECONDS_PER_YEAR
                observed_uy *= SECONDS_PER_YEAR

                uncertainty = None
                if "observedSurfaceVelocityUncertainty" in ds.variables:
                    uncertainty = _as_array(
                        ds.variables["observedSurfaceVelocityUncertainty"],
                        args.time_index,
                        "nCells",
                    )
                    uncertainty *= SECONDS_PER_YEAR
                else:
                    print(
                        "WARNING: 'observedSurfaceVelocityUncertainty' is absent; "
                        "continuing without observed-velocity uncertainty.",
                        file=sys.stderr,
                    )

                results["observed"], regional_results["observed"] = calculate_metrics(
                    observed_ux, observed_uy, uncertainty
                )

            if args.velocity_source in ("modeled", "both"):
                modeled_ux_var = _require_variable(ds, "uReconstructX")
                modeled_uy_var = _require_variable(ds, "uReconstructY")
                modeled_ux_interfaces = _as_cell_interface_array(
                    modeled_ux_var, args.time_index
                )
                modeled_uy_interfaces = _as_cell_interface_array(
                    modeled_uy_var, args.time_index
                )
                if modeled_ux_interfaces.shape != modeled_uy_interfaces.shape:
                    raise ValueError(
                        "Modeled X and Y velocity fields have different shapes: "
                        f"{modeled_ux_interfaces.shape} and {modeled_uy_interfaces.shape}."
                    )
                layer_fractions = _as_array(
                    _require_variable(ds, "layerThicknessFractions"),
                    args.time_index,
                    "nVertLevels",
                )
                modeled_ux = depth_average_velocity(modeled_ux_interfaces, layer_fractions)
                modeled_uy = depth_average_velocity(modeled_uy_interfaces, layer_fractions)
                modeled_ux *= SECONDS_PER_YEAR
                modeled_uy *= SECONDS_PER_YEAR
                results["modeled"], regional_results["modeled"] = calculate_metrics(
                    modeled_ux, modeled_uy
                )
    except (KeyError, ValueError, IndexError, OSError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    print(f"Input file: {args.file}")
    if args.regions_file is not None:
        print(f"Regions file: {args.regions_file}")
        print("Regional masks are evaluated independently and may overlap.")
    print(f"Time index: {args.time_index}")
    if "observed" in results:
        _print_result("Observed surface (plug-flow)", results["observed"])
        if region_names:
            _print_regional_results(
                "Observed surface (plug-flow)", region_names, regional_results["observed"]
            )
    if "modeled" in results:
        _print_result("Modeled depth-averaged", results["modeled"])
        if region_names:
            _print_regional_results(
                "Modeled depth-averaged", region_names, regional_results["modeled"]
            )
    if "observed" in results and "modeled" in results:
        observed_flux = results["observed"].signed_gt_per_year
        modeled_flux = results["modeled"].signed_gt_per_year
        difference = modeled_flux - observed_flux
        print("\nModeled minus observed")
        print(f"  Signed-flux difference:      {difference:.6g} Gt/yr")
        if observed_flux != 0.0:
            print(f"  Relative difference:         {100.0 * difference / observed_flux:.6g}%")
        if region_names:
            print("\nRegional modeled-minus-observed signed flux")
            print(f"  {'Region':<28} {'Difference (Gt/yr)':>20} {'Relative':>12}")
            for name, observed, modeled in zip(
                region_names, regional_results["observed"], regional_results["modeled"]
            ):
                regional_difference = modeled.signed_gt_per_year - observed.signed_gt_per_year
                relative = (
                    f"{100.0 * regional_difference / observed.signed_gt_per_year:.5g}%"
                    if observed.signed_gt_per_year != 0.0 else "--"
                )
                print(f"  {name[:28]:<28} {regional_difference:20.6g} {relative:>12}")
    print("Positive signed discharge is from grounded ice toward floating ice.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
