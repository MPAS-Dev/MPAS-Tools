CXX=g++
#CXXFLAGS= -O3
CXXFLAGS= -O3 -m64
#NETCDF=/opt/netcdf/3.6.3/64

DBGFLAGS= -g -m64

#OFFSETFLAG=-D_64BITOFFSET

platform=$(shell uname)
ifeq ($(platform),Darwin)
	PLATFORM=_MACOS
else
	PLATFORM=_LINUX
endif

ifeq ($(PLATFORM),_LINUX)
	CXXLIBS = -I$(NETCDF)/include -L$(NETCDF)/lib -lGL -lglut -lnetcdf_c++ -lnetcdf -lGLU -lstdc++
endif

ifeq ($(PLATFORM),_MACOS)
	CXXLIBS = -L$(NETCDF)/lib -lnetcdf -lnetcdf_c++ -I$(NETCDF)/include -lstdc++ -framework OpenGL -framework GLUT
endif

all:
	$(CXX) mpas_draw.cpp vec_utils.cpp netcdf_utils.cpp $(CXXLIBS) $(CXXFLAGS) $(OFFSETFLAG) -D$(PLATFORM) -o MpasDraw.x

debug:
	$(CXX) mpas_draw.cpp vec_utils.cpp netcdf_utils.cpp $(CXXLIBS) $(DBGFLAGS) $(OFFSETFLAG) -D$(PLATFORM) -o MpasDraw.x

clean:
	rm MpasDraw.x

