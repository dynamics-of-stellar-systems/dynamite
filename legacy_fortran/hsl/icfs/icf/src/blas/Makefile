#
# BLAS directory.
#

L_ARCH   = $(ARCH)
LIB_NAME = d-$(L_ARCH).a

OPTFLAGS = -O

AR     = ar rcv
RANLIB = ranlib 

FILES = daxpy.o  dcopy.o  ddot.o   dnrm2.o  dscal.o

CFLAGS = $(OPTFLAGS) 
FFLAGS = $(OPTFLAGS)

blas: $(FILES) 
	$(AR) $(LIB_NAME) $(FILES)  
	$(RANLIB) $(LIB_NAME)
	rm -f *.o

clean:
	- rm -f *.o
	- rm -f *.a

.c.o:
	$(CC) $(CFLAGS) -c $*.c
.f.o:
	$(FC) $(FFLAGS) -c $*.f


