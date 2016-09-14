#Assumes ${NETCDF} is defined to the root of the netcdf library

CC=g++
CFLAGS= -O3
EXECUTABLE= PointsToMpas.x

NCINC=$(shell nc-config --cflags)
NCLIB=$(shell nc-config --flibs)

all:
	${CC} points-mpas.cpp ${CFLAGS} -o ${EXECUTABLE} ${NCINC} ${NCLIB} -lnetcdf_c++ -lnetcdf

debug:
	${CC} points-mpas.cpp  -g -o ${EXECUTABLE} ${NCINC} ${NCLIB} -lnetcdf_c++ -lnetcdf

clean:
	rm -f grid.nc
	rm -f graph.info
	rm -f ${EXECUTABLE}

