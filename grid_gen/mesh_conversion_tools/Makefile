#Assumes ${NETCDF} is defined to the root of the netcdf library

CC=g++
CFLAGS= -O3 -std=c++0x
DFLAGS= -g -std=c++0x -D_DEBUG
CONV_EXECUTABLE= MpasMeshConverter.x
CULL_EXECUTABLE= MpasCellCuller.x

NCINCDIR=${NETCDF}/include
NCLIBDIR=${NETCDF}/lib

all:
	${CC} mpas_mesh_converter.cpp netcdf_utils.cpp ${CFLAGS} -o ${CONV_EXECUTABLE} -I${NCINCDIR} -L${NCLIBDIR} -lnetcdf_c++ -lnetcdf -lstdc++
	${CC} mpas_cell_culler.cpp netcdf_utils.cpp ${CFLAGS} -o ${CULL_EXECUTABLE} -I${NCINCDIR} -L${NCLIBDIR} -lnetcdf_c++ -lnetcdf -lstdc++

debug:
	${CC} mpas_mesh_converter.cpp netcdf_utils.cpp ${DFLAGS} -o ${CONV_EXECUTABLE} -I${NCINCDIR} -L${NCLIBDIR} -lnetcdf_c++ -lnetcdf -lstdc++
	${CC} mpas_cell_culler.cpp netcdf_utils.cpp ${DFLAGS} -o ${CULL_EXECUTABLE} -I${NCINCDIR} -L${NCLIBDIR} -lnetcdf_c++ -lnetcdf -lstdc++

clean:
	rm -f grid.nc
	rm -f graph.info
	rm -f ${CONV_EXECUTABLE} ${CULL_EXECUTABLE}

