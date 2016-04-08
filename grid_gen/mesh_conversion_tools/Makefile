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

NCINCDIR=${NETCDF}/include
NCLIBDIR=${NETCDF}/lib

all:
	${CC} mpas_mesh_converter.cpp netcdf_utils.cpp ${CFLAGS} -o ${CONV_EXECUTABLE} -I${NCINCDIR} -L${NCLIBDIR} -lnetcdf_c++ -lnetcdf -lstdc++
	${CC} mpas_cell_culler.cpp netcdf_utils.cpp ${CFLAGS} -o ${CULL_EXECUTABLE} -I${NCINCDIR} -L${NCLIBDIR} -lnetcdf_c++ -lnetcdf -lstdc++
	${CC} mpas_mask_creator.cpp netcdf_utils.cpp jsoncpp.cpp ${CFLAGS} -o ${MASK_EXECUTABLE} -I. -I${NCINCDIR} -L${NCLIBDIR} -lnetcdf_c++ -lnetcdf -lstdc++

debug:
	${CC} mpas_mesh_converter.cpp netcdf_utils.cpp ${DFLAGS} -o ${CONV_EXECUTABLE} -I${NCINCDIR} -L${NCLIBDIR} -lnetcdf_c++ -lnetcdf -lstdc++
	${CC} mpas_cell_culler.cpp netcdf_utils.cpp ${DFLAGS} -o ${CULL_EXECUTABLE} -I${NCINCDIR} -L${NCLIBDIR} -lnetcdf_c++ -lnetcdf -lstdc++
	${CC} mpas_mask_creator.cpp netcdf_utils.cpp jsoncpp.cpp ${DFLAGS} -o ${MASK_EXECUTABLE} -I. -I${NCINCDIR} -L${NCLIBDIR} -lnetcdf_c++ -lnetcdf -lstdc++

clean:
	rm -f grid.nc
	rm -f graph.info
	rm -f ${CONV_EXECUTABLE} ${CULL_EXECUTABLE} ${MASK_EXECUTABLE}

