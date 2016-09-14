.SUFFIXES: .F .f .o

OBJS = STRIPACK.o module_grid_params.o module_grid_constants.o module_data_types.o module_sphere_utilities.o module_voronoi_utils.o module_grid_gen_utils.o module_scvt.o module_write_netcdf.o module_grid_meta.o grid_gen.o

all: $(OBJS)
	$(FC) $(PROMOTION) $(LDFLAGS) -o grid_gen $(OBJS) $(shell nc-config --fflags)


grid_gen.o: module_grid_params.o module_grid_constants.o module_data_types.o module_grid_gen_utils.o module_voronoi_utils.o STRIPACK.o module_scvt.o module_grid_meta.o

module_grid_gen_utils.o: module_sphere_utilities.o

module_scvt.o: module_data_types.o module_sphere_utilities.o module_voronoi_utils.o module_grid_constants.o module_grid_params.o

module_write_netcdf.o: module_grid_params.o

module_data_types.o:

module_grid_meta.o: module_data_types.o module_grid_constants.o module_sphere_utilities.o module_write_netcdf.o

module_sphere_utilities.o: module_data_types.o

module_grid_constants.o:

module_grid_params.o:

module_voronoi_utils.o: module_grid_constants.o STRIPACK.o


.F.o:
	cpp -C -P -traditional $(CPPFLAGS) $< > $*.f90
	$(FC) $(FFLAGS) $(PROMOTION) -c $*.f90
	rm -f $*.f90

.f.o:
	cpp -C -P -traditional $(CPPFLAGS) $< > $*.for
	$(FC) $(F77FLAGS) $(PROMOTION) -c $*.for
	rm -f $.*.for

clean:
	rm -f *.o *.mod grid_gen
