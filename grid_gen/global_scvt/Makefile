#FC = ifort
#FFLAGS = -FR -m64 -O3 -fast -ipo -openmp
#F77FLAGS = -FI -m64 -O3 -fast -ipo -openmp
#CPPFLAGS = -DRKIND=8
#PROMOTION = -r8
#LDFLAGS = -m64 -O3 -fast -ipo -openmp

FC = $(shell nc-config --fc)
FFLAGS = -ffree-form -O3 -fopenmp -ffree-line-length-none $(shell nc-config --fflags)
F77FLAGS = -ffixed-form -O3 -fopenmp -fsecond-underscore
CPPFLAGS = -DRKIND=8 -D_RKIND=_8
PROMOTION = -fdefault-real-8
LDFLAGS = -O3 -fopenmp $(shell nc-config --flibs)

#FC = pgf90
#FFLAGS = -Mfree -O3 -mp -byteswapio
#F77FLAGS = -O3 -byteswapio
#CPPFLAGS = -DRKIND=8
#PROMOTION = -r8
#LDFLAGS = -O3 -mp -byteswapio

all: grid_gen grid_ref

grid_gen:
	( cd src; make FC="$(FC)" FFLAGS="$(FFLAGS)" F77FLAGS="$(F77FLAGS)" CPPFLAGS="$(CPPFLAGS)" PROMOTION="$(PROMOTION)" LDFLAGS="$(LDFLAGS)" )
	if [ ! -e grid_gen ]; then ln -s src/grid_gen .; fi

grid_ref:
	$(FC) refine/grid_ref.f refine/svtgen.f -o refine/grid_ref
	if [ ! -e grid_ref ]; then ln -s refine/grid_ref .; fi

clean:
	( cd src; make clean )
	rm -f grid_gen grid_ref refine/grid_ref

