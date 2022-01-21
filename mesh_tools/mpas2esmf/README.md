# mpas2esmf

This tool generates ESMF and SCRIP files from an MPAS grid file. To avoid
complications with CESM/CAM infrastructure tools, this tool should only be ran
on grids that have a `sphere_radius = 1`. 

To build, ensure `nc-config` is in your $PATH and call `make`.

By default, the ESMF and SCRIP NetCDF files created are 64BIT offset format.
