#Assumes ${NETCDF} is defined to the root of the netcdf library

#### titan instructions ####
# module purge
# module load gcc/4.9.3
# module load cray-netcdf/3.6.2
# printenv NETCDF_DIR
# export NETCDF=$NETCDF_DIR
# make

#### edison instructions ####
# In this file:
# comment gnu, uncomment intel flags
# change to:
#   CXX=CC
# may need to unload parallel NetCDF and HDF5 libraries
# to avoid g++ conflicts:
#   module unload python
# then:
#   module load cray-netcdf
#   module load netcdf-cxx
#   export NETCDF=/usr/common/graphics/netcdf-cxx/4.2
#   make

# gnu
CXX ?= g++
CFLAGS ?= -O3 -std=c++0x -fopenmp -lstdc++
DFLAGS ?= -g -std=c++0x -D_DEBUG -fopenmp -lstdc++

# intel
# CXX=icpc
# CFLAGS= -O3 -std=c++0x -qopenmp -lstdc++
# DFLAGS= -g -std=c++0x -D_DEBUG -qopenmp -lstdc++

CONV_EXECUTABLE= MpasMeshConverter.x
CULL_EXECUTABLE= MpasCellCuller.x
MASK_EXECUTABLE= MpasMaskCreator.x

ifneq (${NETCDF}, )
	ifneq ($(shell which ${NETCDF}/bin/nc-config 2> /dev/null), )
		LIBS = $(shell ${NETCDF}/bin/nc-config --libs)
		INCS = $(shell ${NETCDF}/bin/nc-config --cflags)
	else
		LIBS= -L${NETCDF}/lib
		LIBS += -lnetcdf
		INCS = -I${NETCDF}/include
	endif
else ifneq ($(shell which nc-config 2> /dev/null), )
	LIBS = $(shell nc-config --libs)
	INCS = $(shell nc-config --cflags)
else
	LIBS= -L${NETCDF}/lib
	LIBS += -lnetcdf
	INCS = -I${NETCDF}/include
endif

INCS += -Inetcdf-cxx-4.2/
SRC = netcdf_utils.cpp netcdf-cxx-4.2/ncvalues.cpp netcdf-cxx-4.2/netcdf.cpp

all:
	${CXX} mpas_mesh_converter.cpp ${SRC} ${CFLAGS} -o ${CONV_EXECUTABLE} -I. ${INCS} ${LIBS}
	${CXX} mpas_cell_culler.cpp ${SRC} ${CFLAGS} -o ${CULL_EXECUTABLE} ${INCS} ${LIBS}
	${CXX} mpas_mask_creator.cpp ${SRC} jsoncpp.cpp ${CFLAGS} -o ${MASK_EXECUTABLE} -I. ${INCS} ${LIBS}

debug:
	${CXX} mpas_mesh_converter.cpp ${SRC} ${DFLAGS} -o ${CONV_EXECUTABLE} ${INCS} ${LIBS}
	${CXX} mpas_cell_culler.cpp ${SRC} ${DFLAGS} -o ${CULL_EXECUTABLE} ${INCS} ${LIBS}
	${CXX} mpas_mask_creator.cpp ${SRC} jsoncpp.cpp ${DFLAGS} -o ${MASK_EXECUTABLE} -I. ${INCS} ${LIBS}

clean:
	rm -f grid.nc
	rm -f graph.info
	rm -f ${CONV_EXECUTABLE} ${CULL_EXECUTABLE} ${MASK_EXECUTABLE}

