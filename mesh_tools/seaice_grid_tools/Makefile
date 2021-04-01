#Assumes ${NETCDF} is defined to the root of the netcdf library

F90 = gfortran

F90FLAGS = -O2

ifeq ($(NETCDF_FORTRAN_PATH),)
	INCS = -I${NETCDF}/include
	LIBS = -L${NETCDF}/lib -lnetcdf -lnetcdff
else
	INCS = -I${NETCDF_FORTRAN_PATH}/include -I${NETCDF_C_PATH}/include
	LIBS = -L${NETCDF_FORTRAN_PATH}/lib -L${NETCDF_C_PATH}/lib -lnetcdf -lnetcdff
endif

all:
	${F90} ${F90FLAGS} ${INCS} fix_regrid_output.F90 -o fix_regrid_output.exe ${LIBS}

clean:
	rm fix_regrid_output.exe
