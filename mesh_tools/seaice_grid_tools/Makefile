#Assumes ${NETCDF} is defined to the root of the netcdf library

F90 = gfortran

F90FLAGS = -O2

INCS = -I${NETCDF}/include
LIBS = -L${NETCDF}/lib -lnetcdf -lnetcdff

all:
	${F90} ${F90FLAGS} ${INCS} ${LIBS} fix_regrid_output.F90 -o fix_regrid_output.exe

clean:
	rm fix_regrid_output.exe
