#Assumes ${NETCDF} is defined to the root of the netcdf library

CC=g++
CFLAGS= -O3
EXECUTABLE= PointsToMpas.x

NCINCDIR=${NETCDF}/include
NCLIBDIR=${NETCDF}/lib

all:
	${CC} points-mpas.cpp ${CFLAGS} -o ${EXECUTABLE} -I${NCINCDIR} -L${NCLIBDIR} -lnetcdf_c++ -lnetcdf

debug:
	${CC} points-mpas.cpp  -g -o ${EXECUTABLE} -I${NCINCDIR} -L${NCLIBDIR} -lnetcdf_c++ -lnetcdf

clean:
	rm -f grid.nc
	rm -f graph.info
	rm -f ${EXECUTABLE}

