.. _mpas_to_xdmf:

======================
MPAS to XDMF Converter
======================

The MPAS to XDMF Converter is a tool designed to convert MPAS output files into
XDMF + HDF5 format for visualization in tools like ParaView. This converter
simplifies the process of working with MPAS data by providing a format that is
more compatible with common visualization software.

Overview
========
The MPAS to XDMF Converter allows users to transform MPAS NetCDF files into
XDMF + HDF5 format, which is widely supported by visualization tools like
ParaView. The tool supports slicing along extra dimensions, combining mesh and
time series data, and selecting specific variables for conversion.

Usage
=====
The converter can be used via its command-line interface ``mpas_to_xdmf`` or as
a Python library.

Command-Line Arguments
----------------------
The following arguments are supported:

- ``-m, --mesh``: Path to the MPAS mesh file (required).
- ``-t, --time-series``: Wildcard or list of time series files (optional).
- ``-v, --variables``: List of variables to convert. Special keys include:
  - ``allOnCells``: Variables with dimension ``nCells``.
  - ``allOnEdges``: Variables with dimension ``nEdges``.
  - ``allOnVertices``: Variables with dimension ``nVertices``.
- ``-o, --output-dir``: Directory to save XDMF files (required).
- ``-x, --xtime``: Name of the variable containing time information (optional).
- ``-d, --dim-list``: List of dimensions and indices to slice (e.g.,
  ``nVertLevels=0:10:2``).

Examples
--------
Command-Line Examples
~~~~~~~~~~~~~~~~~~~~~
Convert a mesh file and time series files, selecting specific variables:

.. code-block:: bash

    mpas_to_xdmf -m mesh.nc -t time_series.*.nc -v temperature,salinity -o output_dir

Convert all variables on cells, slicing along ``nVertLevels``:

.. code-block:: bash

    mpas_to_xdmf -m mesh.nc -t time_series.*.nc -v allOnCells -o output_dir -d nVertLevels=0:5

Python Examples
~~~~~~~~~~~~~~~
Use the package programmatically to load data and convert it to XDMF:

.. code-block:: python

    from mpas_tools.viz.mpas_to_xdmf.mpas_to_xdmf import MpasToXdmf

    # Initialize the converter
    converter = MpasToXdmf()

    # Load the mesh and time series files
    converter.load(
        mesh_filename="mesh.nc",
        time_series_filenames=["time_series.0001.nc", "time_series.0002.nc"],
        variables=["temperature", "salinity"],
        xtime_var="xtime"
    )

    # Convert to XDMF format
    converter.convert_to_xdmf(out_dir="output_dir")

Slice along extra dimensions and convert:

.. code-block:: python

    from mpas_tools.viz.mpas_to_xdmf.mpas_to_xdmf import MpasToXdmf

    # Initialize the converter
    converter = MpasToXdmf()

    # Load the mesh and time series files
    converter.load(
        mesh_filename="mesh.nc",
        time_series_filenames="time_series.*.nc",
        variables=["allOnCells"],
        xtime_var="xtime"
    )

    # Define extra dimensions to slice
    extra_dims = {"nVertLevels": [0, 1, 2, 3, 4]}

    # Convert to XDMF format with slicing
    converter.convert_to_xdmf(out_dir="output_dir", extra_dims=extra_dims)

Input and Output
================
The tool supports the following input and output formats:

Input Files
-----------
- **Mesh File**: A NetCDF file containing the MPAS mesh (e.g., ``mesh.nc``).
  This file is also used for the variables to extract if alternate data files
  (e.g. a time series) are not provided.
- **Time Series Files**: Optional NetCDF file(s) containing the data to be
  extracted, often a time series (e.g., ``time_series.*.nc``).

Output Files
------------
- **XDMF Files**: Metadata files describing the structure of the data.
- **HDF5 Files**: Binary files containing the actual data.

The output files are saved in the specified directory, with separate files
for cell-centered, edge-centered, and vertex-centered data.

Features
========
The MPAS to XDMF Converter includes several basic features:

- **Slicing Extra Dimensions**: Users can slice along extra dimensions (e.g.,
  ``nVertLevels``) by specifying indices or ranges.
- **Combining Mesh and Time Series Data**: The tool will merge mesh and time
  series files into a single dataset for conversion.
- **Selective Variable Conversion**: Users can choose specific variables or
  groups of variables (e.g., ``allOnCells``) for conversion.

Opening Files in ParaView
=========================
Once the conversion is complete, you can open the generated XDMF files in
ParaView for visualization. Follow these steps:

1. **Open the XDMF File**: In ParaView, open the ``.xdmf`` file, not the
   ``.h5`` file. The ``.xdmf`` file contains the metadata that links to the
   data stored in the ``.h5`` file.

2. **Choose the Correct Reader**: When opening the ``.xdmf`` file, ParaView
   will prompt you to select one of three readers:
   - **Xdmf3 Reader S**: This reader is optimized for static datasets without
     time information. It is not suitable for time-varying MPAS data.
   - **Xdmf3 Reader T**: This reader is designed for time-varying datasets and
     is the preferred choice for MPAS data converted with this tool.
   - **XDMF Reader**: This is an older reader that may not fully support modern
     XDMF features and should generally be avoided.

   **Recommendation**: Always select the **Xdmf3 Reader T** when prompted.
   This ensures that time-varying data is handled correctly, allowing you to
   explore the temporal evolution of your dataset in ParaView.

3. **Select Fields to Import**: After choosing the reader, ParaView will
   display a list of fields available for import. Uncheck any fields you do
   not wish to view, then click the **Apply** button to load the selected
   fields.

4. **Visualize the Data**: After selecting the correct reader and fields, the
   dataset will load into ParaView. You can then use ParaView's tools to
   visualize and analyze the data.

By following these steps, you can ensure that your MPAS data is correctly
interpreted and visualized in ParaView.

References
==========
- `ParaView Documentation <https://www.paraview.org/documentation/>`_
- `XDMF Format Specification <https://xdmf.org/index.php/Main_Page>`_
- `xarray Documentation <https://docs.xarray.dev/>`_
