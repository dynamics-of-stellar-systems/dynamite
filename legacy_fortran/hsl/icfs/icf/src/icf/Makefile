#
# ICF directory.
#

L_ARCH   = $(ARCH)
LIB_NAME = d-$(L_ARCH).a

OPTFLAGS = -O 

AR     = ar rcv
RANLIB = ranlib 

FILES = dicf.o dicfs.o dpcg.o dsel2.o dssyax.o dstrsol.o \
        ihsort.o insort.o srtdat2.o

CFLAGS = $(OPTFLAGS) 
FFLAGS = $(OPTFLAGS)

icf: $(FILES) 
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


