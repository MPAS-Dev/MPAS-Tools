#Assumes ${NETCDF} is defined to the root of the netcdf library

CC=g++
CFLAGS= -O3
EXECUTABLE= PointsToMpas.x

ifneq ($(NETCDF), )
ifneq ($(shell which ${NETCDF}/bin/nc-config 2> /dev/null), )
	LIBS = $(shell ${NETCDF}/bin/nc-config --libs) -lnetcdf_c++
	INCS = $(shell ${NETCDF}/bin/nc-config --cflags)
else
	LIBS = -L${NETCDF}/lib
	LIBS += -lnetcdf_c++ -lnetcdf -lstdc++
	INCS = -I${NETCDF}/include
endif
else ifneq ($(shell which nc-config 2> /dev/null), )
	LIBS = $(shell nc-config --libs) -lnetcdf_c++
	INCS = $(shell nc-config --cflags)
else
	LIBS = -L${NETCDF}/lib
	LIBS += -lnetcdf_c++ -lnetcdf -lstdc++
	INCS = -I${NETCDF}/include
endif

all:
	${CC} points-mpas.cpp ${CFLAGS} -o ${EXECUTABLE} ${INCS} ${LIBS} -lnetcdf_c++ -lnetcdf

debug:
	${CC} points-mpas.cpp  -g -o ${EXECUTABLE} ${INCS} ${LIBS} -lnetcdf_c++ -lnetcdf

clean:
	rm -f grid.nc
	rm -f graph.info
	rm -f ${EXECUTABLE}

