# IBM with Xlf compilers
#FC = xlf90
#CC = xlc
#FFLAGS = -qrealsize=8 -g -C
#CFLAGS = -g
#LDFLAGS = -g -C

# pgf90
#FC = pgf90
#CC = pgcc
#FFLAGS = -r8 -O3
#CFLAGS = -O3
#LDFLAGS = -O3

# gfortran
FC = gfortran
CC = gcc
FFLAGS = -O3 -m64 -ffree-line-length-none -fdefault-real-8 -fconvert=big-endian
CFLAGS = -O3 -m64
LDFLAGS = -O3 -m64

# ifort
#FC = ifort
#CC = icc
#FFLAGS = -real-size 64 -O3
#CFLAGS = -O3
#LDFLAGS = -O3

# absoft
#FC = f90
#CC = gcc
#FFLAGS = -dp -O3
#CFLAGS = -O3
#LDFLAGS = -O3


CPP = cpp -P -traditional
CPPFLAGS = 
CPPINCLUDES = 
INCLUDES = -I$(NETCDF)/include

# Specify NetCDF libraries, checking if netcdff is required (it will be present in v4 of netCDF)
LIBS = -L$(NETCDF)/lib
NCLIB = -lnetcdf
NCLIBF = -lnetcdff
ifneq ($(wildcard $(NETCDF)/lib/libnetcdff.*), ) # CHECK FOR NETCDF4
        LIBS += $(NCLIBF)
endif # CHECK FOR NETCDF4
LIBS += $(NCLIB)



RM = rm -f

##########################

.SUFFIXES: .F .o


OBJS = periodic_grid.o \
       module_cell_indexing.o \
       module_write_netcdf.o

all: periodic_grid

periodic_grid.o: module_cell_indexing.o module_write_netcdf.o 

periodic_grid: $(OBJS)
	$(FC) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	$(RM) *.o *.mod periodic_grid

.F.o:
	$(RM) $@ $*.mod
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(INCLUDES)
	$(RM) $*.f90
