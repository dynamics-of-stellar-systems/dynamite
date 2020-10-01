#
#   ICF directory
#

FFLAGS = -O

L_ARCH = $(ARCH)

LIB_NAME = d-$(L_ARCH).a

BLAS = src/blas/$(LIB_NAME)

MINPACK2 = src/icf/$(LIB_NAME)

libs: blas icf 

blas: 
	cd src/blas; make

icf: 
	cd src/icf; make

#   Files for the MINPACK-2 incomplete Cholesky factorization.

LIBS = $(MINPACK2) $(BLAS) 

SRC = dmain.o dtimer.o dnrm2.o mmio.o

icfx : $(SRC) 
	$(FC) $(FFLAGS) $(SRC) $(LIBS) -o icf

