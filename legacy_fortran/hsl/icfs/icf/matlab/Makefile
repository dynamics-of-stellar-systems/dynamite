#
#   ICF directory
#

L_ARCH = $(ARCH)

LIB_NAME = d-$(L_ARCH).a

BLAS = ../src/blas/$(LIB_NAME)

MINPACK2 = ../src/icf/$(LIB_NAME)

LIBS = $(MINPACK2) $(BLAS) 

libs: blas icf 

blas: 
	cd ../src/blas; make

icf: 
	cd ../src/icf; make

CMEX = mex
CMEXFLAGS = -DARCH=$(ARCH)

# Path for mex.h

MINCLUDE = /software/solaris-2/com/packages/matlab-5.2/extern/include
# MINCLUDE = /software/aix-4/com/packages/matlab-5.2/extern/include

matlab: icfmex cgmex

icfmex : icfmex.c 
	$(CMEX) $(CMEXFLAGS) -I$(MINCLUDE) icfmex.c $(LIBS) 

cgmex : cgmex.c 
	$(CMEX) $(CMEXFLAGS) -I$(MINCLUDE) cgmex.c $(LIBS) 


