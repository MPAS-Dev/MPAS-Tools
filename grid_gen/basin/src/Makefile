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

## ifort
 FC = ifort
 CC = icc
 FFLAGS = -real-size 64 #-g -traceback -check all
 CFLAGS = #-g
 LDFLAGS = #-g -traceback -check all

# gfortran
#FC = gfortran
#CC = gcc
#FFLAGS = -O3 -m64 -ffree-line-length-none -fdefault-real-8 -fconvert=big-endian -ffree-form
#CFLAGS =
#LDFLAGS =

# absoft
#FC = f90
#CC = gcc
#FFLAGS = -dp -O3
#CFLAGS = -O3
#LDFLAGS = -O3
#NETCDF = /Users/maltrud/local


CPP = cpp -C -P -traditional
CPPFLAGS = 
CPPINCLUDES = 
INCLUDES = -I$(NETCDF)/include

LIBS =  -L$(NETCDF)/lib
NCLIB = -lnetcdf
NCLIBF = -lnetcdff
ifneq ($(wildcard $(NETCDF)/lib/libnetcdff.*), ) # CHECK FOR NETCDF4
	LIBS += $(NCLIBF)
endif # CHECK FOR NETCDF4
LIBS += $(NCLIB)

RM = rm -f

##########################

.SUFFIXES: .F .o


OBJS = basin.o \
       utilities.o \
       module_read_netcdf.o \
       module_read_topo.o \
       module_read_TS.o \
       module_read_monthly.o \
       module_cullLoops.o \
       module_write_netcdf.o

all: map

basin.o: utilities.o module_write_netcdf.o module_read_netcdf.o module_read_topo.o module_read_TS.o module_read_monthly.o module_cullLoops.o

map: $(OBJS)
	$(FC) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	$(RM) *.o *.mod pop *.f90

.F.o:
	$(RM) $@ $*.mod
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(INCLUDES)
	#$(RM) $*.f90
