FC = $(shell nf-config --fc)
FFLAGS = -O3
FCINCLUDES = $(shell nf-config --fflags)
FCLIBS = -L$(shell nc-config --libdir) $(shell nf-config --flibs)

all: grid_rotate

grid_rotate: grid_rotate.f90
	$(FC) grid_rotate.f90 -o grid_rotate $(FFLAGS) $(FCINCLUDES) $(FCLIBS)

clean:
	rm grid_rotate
