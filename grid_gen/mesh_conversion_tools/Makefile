#Assumes ${NETCDF} is defined to the root of the netcdf library

# gfortran
CC=g++
CFLAGS= -O3 -std=c++0x -fopenmp -lstdc++
DFLAGS= -g -std=c++0x -D_DEBUG -fopenmp -lstdc++

# intel
# CC=icpc
# CFLAGS= -O3 -std=c++0x -openmp -lstdc++
# DFLAGS= -g -std=c++0x -D_DEBUG -openmp -lstdc++

CONV_EXECUTABLE= MpasMeshConverter.x
CULL_EXECUTABLE= MpasCellCuller.x
MASK_EXECUTABLE= MpasMaskCreator.x

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

all:
	${CC} mpas_mask_creator.cpp netcdf_utils.cpp jsoncpp.cpp ${CFLAGS} -o ${MASK_EXECUTABLE} -I. ${INCS} ${LIBS}
	${CC} mpas_mesh_converter.cpp netcdf_utils.cpp ${CFLAGS} -o ${CONV_EXECUTABLE} ${INCS} ${LIBS}
	${CC} mpas_cell_culler.cpp netcdf_utils.cpp ${CFLAGS} -o ${CULL_EXECUTABLE}  ${INCS} ${LIBS}

debug:
	${CC} mpas_mesh_converter.cpp netcdf_utils.cpp ${DFLAGS} -o ${CONV_EXECUTABLE} ${INCS} ${LIBS}
	${CC} mpas_cell_culler.cpp netcdf_utils.cpp ${DFLAGS} -o ${CULL_EXECUTABLE} ${INCS} ${LIBS}
	${CC} mpas_mask_creator.cpp netcdf_utils.cpp jsoncpp.cpp ${DFLAGS} -o ${MASK_EXECUTABLE} -I. ${INCS} ${LIBS}

clean:
	rm -f grid.nc
	rm -f graph.info
	rm -f ${CONV_EXECUTABLE} ${CULL_EXECUTABLE} ${MASK_EXECUTABLE}

