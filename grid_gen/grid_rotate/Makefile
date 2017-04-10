FC = gfortran
FFLAGS = -m64 -O2 -ffree-line-length-none -ffree-form -Wall

ifneq ($(NETCDF), )
ifneq ($(shell which ${NETCDF}/bin/nc-config 2> /dev/null), )
	LIBS = $(shell ${NETCDF}/bin/nc-config --libs) -lnetcdf_c++ -lnetcdff
	INCS = $(shell ${NETCDF}/bin/nc-config --cflags)
else
	LIBS= -L${NETCDF}/lib
	LIBS += -lnetcdf_c++ -lnetcdf -lstdc++ -lnetcdff
	INCS = -I${NETCDF}/include
endif
else ifneq ($(shell which nc-config 2> /dev/null), )
	LIBS = $(shell nc-config --libs) -lnetcdf_c++
	INCS = $(shell nc-config --cflags)
else
	LIBS= -L${NETCDF}/lib
	LIBS += -lnetcdf_c++ -lnetcdf -lstdc++ -lnetcdff
	INCS = -I${NETCDF}/include
endif

all: grid_rotate

grid_rotate: grid_rotate.f90
	$(FC) grid_rotate.f90 -o grid_rotate $(FFLAGS) $(INCS) $(LIBS)

clean:
	rm grid_rotate
