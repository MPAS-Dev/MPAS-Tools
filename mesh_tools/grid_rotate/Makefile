FC = $(shell nf-config --fc)
FFLAGS = -O3
FCINCLUDES = $(shell nf-config --fflags)
RPATH_FLAGS = $(shell nf-config --flibs | grep -o -e '-L\S\+\( \|$$\)' | sed 's/^-L/-Wl,-rpath,/' | tr -d '\n')
RPATH_FLAGS += $(shell nc-config --libs | grep -o -e '-L\S\+\( \|$$\)' | sed 's/^-L/-Wl,-rpath,/' | tr -d '\n')
FCLIBS = -L$(shell nc-config --libdir) $(shell nf-config --flibs) $(RPATH_FLAGS)

all: grid_rotate

grid_rotate: grid_rotate.f90
	$(FC) grid_rotate.f90 -o grid_rotate $(FFLAGS) $(FCINCLUDES) $(FCLIBS)

clean:
	rm grid_rotate
