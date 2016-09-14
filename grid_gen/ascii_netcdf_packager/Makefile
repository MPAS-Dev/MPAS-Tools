CXX = g++
CPPFLAGS = $(shell nc-config --cflags)
CXXFLAGS = -O3
LIBS = $(shell nc-config --libs) -lnetcdf_c++
EXE = AsciiNetCDFPackager.x

all:
	$(CXX) $(CPPFLAGS) ascii_to_netcdf_packager.cpp $(CXXFLAGS) $(LIBS) -o $(EXE)

clean:
	rm -f $(EXE)
