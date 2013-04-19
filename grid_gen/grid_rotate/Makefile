all: grid_rotate

grid_rotate: grid_rotate.f90
	gfortran grid_rotate.f90 -o grid_rotate -m64 -ffree-line-length-none -I$(NETCDF)/include  -L$(NETCDF)/lib -lnetcdf

clean:
	rm grid_rotate
