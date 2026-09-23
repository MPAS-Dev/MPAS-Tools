from argparse import ArgumentParser, RawDescriptionHelpFormatter
from mpas_tools.io import write_netcdf
import xarray as xr
import numpy as np
from scipy.spatial import cKDTree


parser = ArgumentParser(description=__doc__,
                       formatter_class=RawDescriptionHelpFormatter)
parser.add_argument('-i', '--input_file', dest='input_file', required=True,
                    metavar='FILENAME',
                    help='MALI file (NetCDF format) containing muFriction')
parser.add_argument('-o', '--output_file', dest='output_file', required=True,
                    metavar='FILENAME',
                    help='Destination file for converted muFriction')
parser.add_argument('-d', '--direction', dest='direction', required=True,
                    choices=['weertman_to_budd', 'w2b', 'budd_to_weertman', 'b2w'],
                    help='Conversion direction: weertman_to_budd (w2b) or budd_to_weertman (b2w)')
args = parser.parse_args()

input_filename = args.input_file
output_filename = args.output_file
direction = args.direction
assert input_filename != output_filename, \
    "Input file and output file must have different names!"

mu_friction_max = 100  # Replace with your desired threshold
fill_method = "idw"      # Either "idw" or "nearest"
idw_neighbors = 8
idw_power = 2.0

def fill_large_values(
    field,
    x_cell,
    y_cell,
    maximum,
    method="idw",
    neighbors=8,
    power=2.0,
):
    """
    Replace values greater than `maximum` using spatial interpolation.

    Parameters
    ----------
    field : xr.DataArray
        Field containing an ``nCells`` dimension.
    x_cell, y_cell : xr.DataArray
        MPAS cell-center coordinates.
    maximum : float
        Values greater than this threshold will be replaced.
    method : {"idw", "nearest"}
        Spatial filling method.
    neighbors : int
        Number of neighbors used for IDW.
    power : float
        Distance exponent used for IDW.
    """
    if "nCells" not in field.dims:
        raise ValueError("field must contain an 'nCells' dimension")

    if method not in {"idw", "nearest"}:
        raise ValueError("method must be either 'idw' or 'nearest'")

    # Put nCells last, allowing the function to handle optional Time or
    # other leading dimensions.
    original_dims = field.dims
    field_work = field.transpose(
        *[dim for dim in field.dims if dim != "nCells"],
        "nCells",
    )

    values = np.asarray(field_work.values).copy()
    original_shape = values.shape
    values_2d = values.reshape(-1, original_shape[-1])

    coordinates = np.column_stack(
        [np.asarray(x_cell.values), np.asarray(y_cell.values)]
    )

    for row in values_2d:
        target_mask = np.isfinite(row) & (row > maximum)
        donor_mask = np.isfinite(row) & (row <= maximum)

        if not np.any(target_mask):
            continue

        if not np.any(donor_mask):
            raise ValueError(
                "No valid muFriction cells are available for interpolation"
            )

        donor_values = row[donor_mask]
        donor_coordinates = coordinates[donor_mask]
        target_coordinates = coordinates[target_mask]

        tree = cKDTree(donor_coordinates)

        if method == "nearest":
            _, indices = tree.query(target_coordinates, k=1)
            row[target_mask] = donor_values[indices]

        else:
            neighbor_count = min(neighbors, donor_values.size)
            distances, indices = tree.query(
                target_coordinates,
                k=neighbor_count,
            )

            # Ensure two-dimensional arrays when k=1.
            distances = np.atleast_2d(distances)
            indices = np.atleast_2d(indices)

            if distances.shape[0] != target_coordinates.shape[0]:
                distances = distances.T
                indices = indices.T

            # A tiny distance floor prevents division by zero.
            weights = 1.0 / np.maximum(distances, 1.0e-12) ** power
            neighbor_values = donor_values[indices]

            row[target_mask] = np.sum(
                weights * neighbor_values, axis=1
            ) / np.sum(weights, axis=1)

    filled = xr.DataArray(
        values_2d.reshape(original_shape),
        dims=field_work.dims,
        coords=field_work.coords,
        attrs=field.attrs,
        name=field.name,
    )

    return filled.transpose(*original_dims)


# ----------------------------------------------------------------------
# Load and process the data
# ----------------------------------------------------------------------

ds = xr.open_dataset(input_filename)

thickness = ds["thickness"]
bed_topography = ds["bedTopography"]
mu_friction = ds["muFriction"]

# Compute max(0, -(rho_ocean / rho_ice) * bedTopography).
rho_ocean = 1028.0
rho_ice = 910.0
gravity = 9.81

inner_term = -(rho_ocean / rho_ice) * bed_topography
inner_max = xr.where(inner_term > 0.0, inner_term, 0.0)

denominator = 1.e-3 * (  # account for kPa units in Albany
    rho_ice * gravity * thickness
    - rho_ice * gravity * inner_max
)

denominator = xr.where(
    denominator > 1.0e-12,
    denominator,
    1.0e-12,
)

if direction in ('weertman_to_budd', 'w2b'):
    # Convert from Weertman to Budd by dividing by effective pressure
    mu_friction /= denominator
    # First fill anomalously large muFriction values.
    mu_friction_filled = fill_large_values(
        field=mu_friction,
        x_cell=ds["xCell"],
        y_cell=ds["yCell"],
        maximum=mu_friction_max,
        method=fill_method,
        neighbors=idw_neighbors,
        power=idw_power,
    )
else:  # budd_to_weertman or b2w
    # Convert from Budd to Weertman by multiplying by effective pressure
    mu_friction *= denominator
    # No need to fill large values when converting to Weertman
    mu_friction_filled = mu_friction

# Update the dataset.
ds["muFriction"] = mu_friction_filled
ds["muFriction"].attrs.update(mu_friction.attrs)

#ds.to_netcdf(
#    output_filename,
#    engine="netcdf4",
#    format="NETCDF3_64BIT_OFFSET"
#)

write_netcdf(ds, output_filename, format='NETCDF3_64BIT_DATA')
ds.close()

print(f"Processing complete. Output saved to {output_filename!r}.")
