FC = $(shell nc-config --fc)
FCINCLUDES = $(shell nc-config --fflags)
FCLIBS = $(shell nc-config --flibs)

all: mpas2esmf.f90
	$(FC) -o mpas2esmf mpas2esmf.f90 ${FCINCLUDES} ${FCLIBS}

clean:
	rm -f mpas2esmf read_mesh.mod write_desc.mod
