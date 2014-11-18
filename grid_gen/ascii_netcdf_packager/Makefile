CXX = g++
CPPFLAGS = -I${NETCDF}/include
CXXFLAGS = -O3
LIBS = -L${NETCDF}/lib -lnetcdf -lnetcdf_c++
EXE = AsciiNetCDFPackager.x

all:
	$(CXX) $(CPPFLAGS) ascii_to_netcdf_packager.cpp $(CXXFLAGS) $(LIBS) -o $(EXE)

clean:
	rm -f $(EXE)
