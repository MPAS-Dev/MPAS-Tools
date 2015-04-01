Sample grid generation workflow.

There are many different tools and choices for making an MPAS LI grid.  This is
a sample workflow for making a 20km resolution, uniform density, planar, GIS mesh.

NOTE: All (?) of the python scripts below have up-to-date usage instructions
if you invoke with: python <script.py> --help

Note most of the example commands below are using absolute paths to things on Matt's computer...




1. Make MPAS grid
Here using periodic_hex tool for uniform density planar mesh.
Location: in MPAS-Tools repo at grid_gen/periodic_hex
Run 'make' to compile it.

Set namelist.input for the mesh dimensions required for a 20km GIS mesh:
&periodic_grid
   nx = 90,
   ny = 156,
! ny is 135 in square grid cells, so divide by sqrt(3)/2 = 156
   dc = 20000.,
   nVertLevels = 1,
   nTracers = 1,
   nproc = 2
/
Run:
./periodic_grid


2. Strip out periodic cells - ESMF can't handle that! (This step only needed if using ESMF regridding)
In same directory.
> python mark_periodic_boundaries_for_culling.py grid.nc
Then use the cell culler:
~/documents/mpas-git/Tools/grid_gen/mesh_conversion_tools/MpasCellCuller.x grid.nc
(this tool has to be compiled - just run make)


3. Make an empty MPAS LI mesh with 10 vertical levels
python ~/documents/mpas-git/Tools/grid_gen/landice_grid_tools/create_landice_grid_from_generic_MPAS_grid.py --beta -l 10 -i culled_mesh.nc -o landice_grid.nc


4. Translate the empty mesh to the correct location
python /Users/mhoffman/documents/mpas-git/Tools/grid_gen/planar_grid_transformations/translate_planar_grid.py -f landice_grid.nc -d /Users/mhoffman/documents/greenland_geometry/10-5-2-GIS/gis_10km.ice2sea.050213.nc


5. Add lat/lon values for the appropriate projection (only necessary if using ESMF regridding, but a good idea to do anyway)
python ~/documents/mpas-git/Tools/grid_gen/planar_grid_transformations/set_lat_lon_fields_in_planar_grid.py -f landice_grid.nc -p gis-bamber
Invoke with -h to see available coordinate systems.  More can be added as needed.


6. interpolate data from CISM grid to MPAS grid
python ~/documents/mpas-git/Tools/grid_gen/landice_grid_tools/interpolate_cism_grid_to_mpas_grid.py -c /Users/mhoffman/documents/greenland_geometry/10-5-2-GIS/gis_10km.ice2sea.050213.nc -m landice_grid.nc
Note: if you leave off the -w option (and the corresponding filename) as done here, this script will use a bilinear interpolation scheme written in python.  The bilinear interpolation can have issues if the CISM mesh is much finer than the MPAS mesh because it will only use a few of the CISM cells that within the MPAS cells to determine the value at the MPAS cell.  This can be unrepresentative by sampling high-frequency topographic variations.  However, if you are going to crop the mesh to the extent of Greenland (next steps), you can get away with the potentially less accurate bilinear interpolation during this step.  You might, however, prefer to use the ESMF method when you re-populate the IC in later steps.
The interpolation script will try to interpolate a bunch of different CISM fields if they exist.  For example, it will try to interpolate temp and beta.  If you want it to try tempstag instead, you have to comment/uncomment some lines in the dictionary that matches up the CISM to MPAS field names.  It does do vertical interpolation as well, so the CISM and MPAS vertical levels do NOT need to match!
Note that for now you only really need thickness, so if this step is slow, you can modify the script to only interpolate thickness.


6. setup mask for culling
python ~/documents/mpas-git/Tools/grid_gen/landice_grid_tools/define_cullMask.py
This script has some hard-coded options you might want to look at in the code.  But by default it will mark for culling all cells that are more than 5 cells outside of the current ice extent.  (It will keep the extent of the ice sheet plus 5 more cells as a buffer for ice to grow into.)  This script should be refined as we start to use it more (e.g., keeping cells within a specified distance rather than a number of cells).


7. cull it
~/documents/mpas-git/Tools/grid_gen/mesh_conversion_tools/MpasCellCuller.x landice_grid.nc culled_landice_grid.nc
(this tool has to be compiled - just run make)


8. convert to LI grid - the culling tool spits out a bare-minimum MPAS mesh, so we now need to repeat a number of steps...
python ~/documents/mpas-git/Tools/grid_gen/landice_grid_tools/create_landice_grid_from_generic_MPAS_grid.py --beta -l 10 -i culled_landice_grid.nc -o gis20km.nc
NOTE: Don't need re-add lat/lon values!  Those are copied over by the cell culler.


9. Optional: If wanting to use ESMF regridder
  A. Create SCRIP file for the MPAS mesh
python ~/documents/mpas-git/Tools/grid_gen/create_SCRIP_files/create_SCRIP_file_from_MPAS_mesh.py -m gis20km.nc -s landice_grid.SCRIP

  B. Create SCRIP file from CISM mesh
python ~/documents/mpas-git/Tools/grid_gen/create_SCRIP_files/create_SCRIP_file_from_CISM_mesh.py -c /Users/mhoffman/documents/greenland_geometry/10-5-2-GIS/gis_10km.ice2sea.050213.nc -s cism.scrip.nc -p gis-bamber

  C. Run ESMF to create weight file
~/software/esmf/esmf_6.3.0rp1/DEFAULTINSTALLDIR/bin/binO/Darwin.gfortran.64.mpiuni.default/ESMF_RegridWeightGen --source cism.scrip.nc --destination landice_grid.SCRIP  --weight weight-c.nc --method conserve
See:
https://www.earthsystemcog.org/projects/esmf/regridding_esmf_6_3_0rp1
http://www.earthsystemmodeling.org/esmf_releases/public/ESMF_6_3_0rp1/ESMF_refdoc/node3.html#SECTION03020000000000000000
Note: ESMF is available on Mustang if you don't want to deal with installing it somewhere.


10. interpolate data from CISM grid to the culled MPAS grid
Using ESMF weight file:
python ~/documents/mpas-git/Tools/grid_gen/landice_grid_tools/interpolate_cism_grid_to_mpas_grid.py -c /Users/mhoffman/documents/greenland_geometry/10-5-2-GIS/gis_10km.ice2sea.050213.nc -m gis20km.nc -w weight-c-gis20km.nc

Using built-in bilinear interp:
python ~/documents/mpas-git/Tools/grid_gen/landice_grid_tools/interpolate_cism_grid_to_mpas_grid.py -c /Users/mhoffman/documents/greenland_geometry/10-5-2-GIS/gis_10km.ice2sea.050213.nc -m gis20km.nc


11. make decomposition graph files so you can run on multiple processors
mv culled_graph.info graph.info
gpmetis graph.info 4   (etc.)
(gpmetis should be available on most clusters.  It can also be installed through MacPorts.)



Whew!  That should be it!  If you don't have data for fields like temperature and beta, you'll have to next figure out how to initialize those... but that's another topic.



