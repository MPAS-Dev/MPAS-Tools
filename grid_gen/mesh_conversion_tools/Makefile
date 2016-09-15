#Assumes ${NETCDF} is defined to the root of the netcdf library

# gfortran
CC=g++
CFLAGS= -O3 -std=c++0x -fopenmp
DFLAGS= -g -std=c++0x -D_DEBUG -fopenmp

# intel
# CC=icpc
# CFLAGS= -O3 -std=c++0x -openmp
# DFLAGS= -g -std=c++0x -D_DEBUG -openmp

CONV_EXECUTABLE= MpasMeshConverter.x
CULL_EXECUTABLE= MpasCellCuller.x
MASK_EXECUTABLE= MpasMaskCreator.x

NETCDFLIB = $(shell nc-config --libs) -lnetcdf_c++
NETCDFINC = $(shell nc-config --cflags)

all:
	${CC} mpas_mask_creator.cpp netcdf_utils.cpp jsoncpp.cpp ${CFLAGS} -o ${MASK_EXECUTABLE} -I. ${NETCDFINC} ${NETCDFLIB} -lstdc++
	${CC} mpas_mesh_converter.cpp netcdf_utils.cpp ${CFLAGS} -o ${CONV_EXECUTABLE} ${NETCDFINC} ${NETCDFLIB} -lstdc++
	${CC} mpas_cell_culler.cpp netcdf_utils.cpp ${CFLAGS} -o ${CULL_EXECUTABLE}  ${NETCDFINC} ${NETCDFLIB} -lstdc++

debug:
	${CC} mpas_mesh_converter.cpp netcdf_utils.cpp ${DFLAGS} -o ${CONV_EXECUTABLE} ${NETCDFINC} ${NETCDFLIB} -lstdc++
	${CC} mpas_cell_culler.cpp netcdf_utils.cpp ${DFLAGS} -o ${CULL_EXECUTABLE} ${NETCDFINC} ${NETCDFLIB} -lstdc++
	${CC} mpas_mask_creator.cpp netcdf_utils.cpp jsoncpp.cpp ${DFLAGS} -o ${MASK_EXECUTABLE} -I. ${NETCDFINC} ${NETCDFLIB} -lstdc++

clean:
	rm -f grid.nc
	rm -f graph.info
	rm -f ${CONV_EXECUTABLE} ${CULL_EXECUTABLE} ${MASK_EXECUTABLE}

