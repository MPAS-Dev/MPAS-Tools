FC = $(shell nf-config --fc)
FFLAGS = -O3
FCINCLUDES = $(shell nf-config --fflags)
RPATH_FLAGS = $(shell nf-config --flibs | grep -o -e '-L\S\+\( \|$$\)' | sed 's/^-L/-Wl,-rpath,/' | tr -d '\n')
RPATH_FLAGS += $(shell nc-config --libs | grep -o -e '-L\S\+\( \|$$\)' | sed 's/^-L/-Wl,-rpath,/' | tr -d '\n')
FCLIBS = -L$(shell nc-config --libdir) $(shell nf-config --flibs) $(RPATH_FLAGS)

OBJS = project_hexes.o \
       mpas_geometry_utilities.o \
       ph_utils.o \
       ph_write_mesh.o \
       projection_setup.o \
       projections.o \
       precision.o

EXE = project_hexes

$(EXE) : $(OBJS)
	$(FC) -o $@ $(OBJS) $(FCLIBS)

project_hexes.o : mpas_geometry_utilities.o ph_utils.o ph_write_mesh.o projection_setup.o projections.o precision.o

mpas_geometry_utilities.o : precision.o

ph_utils.o: mpas_geometry_utilities.o precision.o

ph_write_mesh.o : precision.o

projection_setup.o : precision.o

projections.o : projection_setup.o precision.o

clean:
	$(RM) $(OBJS) *.mod $(EXE)


%.o : %.mod

%.o : %.f90
	$(RM) $@ $*.mod
	$(FC) $(FFLAGS) -c $< $(FCINCLUDES)
