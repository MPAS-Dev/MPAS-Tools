.SUFFIXES: .cxx .o

OMP =
CXX = g++
CXXFLAGS = -O3 $(OMP) -I$(NETCDF)/include

OBJS = main.o Point.o Triangle.o DensityFunction.o PointSet.o

all: libfortune.a
	( $(MAKE) periodic_general )

libfortune.a:
	( cd fortune; $(MAKE) all "CC=$(CXX)" "CFLAGS=$(CXXFLAGS)")
	( ln -s fortune/libfortune.a . )

periodic_general: $(OBJS)
	$(CXX) -O3 $(OMP) -o periodic_general $(OBJS) -L$(NETCDF)/lib -L. -lnetcdf -lfortune

clean:
	( cd fortune; $(MAKE) clean )
	rm -f *.o periodic_general libfortune.a

.cxx.o:
	$(CXX) $(CXXFLAGS) -c $< -I./fortune
