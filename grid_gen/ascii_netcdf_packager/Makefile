#If used, assumes ${NETCDF} is defined to the root of the netcdf library

CXX = g++
CXXFLAGS = -O3

EXE = AsciiNetCDFPackager.x


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
	$(CXX) $(INCS) ascii_to_netcdf_packager.cpp $(CXXFLAGS) $(LIBS) -o $(EXE)

clean:
	rm -f $(EXE)
