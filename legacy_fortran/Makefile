
##############################################################
#
# This file can be used to compile the executables using any
# standard Fortran 90 compiler (e.g. the Linux INTEL Fortran 95)
# This makefile was written by Remco van den Bosch, using the
# Leiden, 2009
#
##############################################################
FAST=TRUE
#PROFILE=TRUE
compiler=gfortran

#PROFILEDIR=~/profile/

ifeq ($(compiler),ifort)
ifdef FAST
  #flags = -fast
  flags += -O3 -IPO -xT -prefetch -pad
  flags += -fno-fnalias
else
  flags += -C -g -inline-debug-info
endif



ifdef PROFILEDIR
ifdef PROFILE
flags  += -prof-dir $(PROFILEDIR)
flags  +=  -prof-gen
else
  ifdef FAST
    flags  += -prof-dir $(PROFILEDIR)
    flags  += -prof-use
  endif
endif
endif


fortran90 = $(compiler) $(flags)     -c
fortran77 = $(compiler) $(flags) -W0 -c
link      = $(compiler) $(flags)     -o
endif

ifeq ($(compiler),F)
  fortran90 = F -ieee=stop -C -O4   -c
  fortran77 = g77.old -O2   -c
  link = F -ieee=stop -C -O4   -o
endif

ifeq ($(compiler),g95)
  fortran90 = ~/local/bin/g95 -c -O3 -m64  -c
  fortran77 = ~/local/bin/g95 -c -O3 -m64  -c
  link = ~/local/bin/g95 -o -O3 -m64 -o
endif

ifeq ($(compiler),gfortran)
  ifdef FAST
   flags =    -ffast-math -O3 -march=native -fomit-frame-pointer -m64
   flags +=    -funroll-loops -ftree-loop-linear
#   flags += -fno-protect-parens
#   flags +=  -fwhole-file #-fwhole-program # -flto
   flags += -std=legacy
  else
   flags =  -m64 -g  -ggdb -fbounds-check -fcheck-array-temporaries
   #flags = -m64   -g -ffast-math -O3 -march=native -fomit-frame-pointer -m64
  endif

  ifdef PROFILEDIR
  ifdef PROFILE
    flags  += -fprofile-generate=$(PROFILEDIR)
  else
    ifdef FAST
        flags  += -fprofile-use=$(PROFILEDIR)
    endif
  endif
  endif

   # on MAC OS X we can use the built-in BLAS/LAPACK
  #flags += -fexternal-blas -framework Accelerate
  flags += -fexternal-blas

  fortran90 = gfortran $(flags) -c
  fortran77 = gfortran $(flags) -c
  link = gfortran $(flags) -o
endif

# define local Galahad directory --> might need to be changed for your system
GALAHADDIR =$(CURDIR)/galahad-2.3

# Define galahad compiled libary choice (platform dependant)

# MAC, gfortran
GALAHADTYPE= mac.osx.gfo/double/

# linux, ifort
# GALAHADTYPE= pc.lnx.ifc/double
#GALAHADTYPE= pc.lnx.gfo/double
##########################################################

all : orbitstart orbitstart_bar orblib_new_mirror orblib_bar triaxmass triaxmass_bar triaxmassbin triaxmassbin_bar triaxnnls_CRcut triaxnnls_noCRcut triaxnnls_bar
nogal : orbitstart orbitstart_bar orblib_new_mirror orblib_bar triaxmass triaxmassbin triaxmass_bar triaxmassbin_bar

#########################################################
#
# Store a copy of all compiled version in RCS as backup
#
#########################################################

RCS : Makefile iniparam_f.f90 orblib_f_new_mirror.f90 triaxpotent.f90 interpolpotent.f90 orblibprogram.f90 triaxmassbin.f90 triaxmass.f90 orbitstart.f90 triaxnnls_CRcut.f90 triaxnnls_noCRcut.f90 dmpotent.f90
	ci -q -l -mcompile Makefile iniparam_f.f90 orblib_f_new_mirror.f90 triaxpotent.f90 interpolpotent.f90 orblibprogram.f90 triaxmassbin.f90 triaxmass.f90 orbitstart.f90 triaxnnls_CRcut.f90 triaxnnls_noCRcut.f90 dmpotent.f90

##########################################################
#
# Choice between two integrators
#
# DOP853 is the more accurate and is less prone to fail
# integrating the orbit trajectory, but it is slower too.
#
##########################################################

#DOPCDE = sub/dopri5.f
#DOPOBJ = dopri5.o

DOPCDE = sub/dop853.f
DOPOBJ = dop853.o


######################################################
#
# MGE Code.
#
######################################################

MGECDE = ran1_nr.f sub/numeric_kinds_f.f90 sub/dqxgs.f sub/ellipint.f90 iniparam_f.f90 triaxpotent.f90 dmpotent.f90 interpolpotent.f90 sub/specfunc_beta.f90
MGEOBJ = ran1_nr.o numeric_kinds_f.o dqxgs.o iniparam_f.o triaxpotent.o dmpotent.o interpolpotent.o ellipint.o specfunc_beta.o

iniparam_f.o: iniparam_f.f90 numeric_kinds_f.o
	$(fortran90) iniparam_f.f90

triaxpotent.o : triaxpotent.f90 ellipint.o numeric_kinds_f.o iniparam_f.o
	$(fortran90) triaxpotent.f90

triaxmass_f.o : nag.o triaxpotent.f90
	$(fortran90) triaxmass_f.f90

triaxmassbin_f.o : nag.o triaxpotent.f90
	$(fortran90) triaxmassbin_f.f90

dmpotent.o : dmpotent.f90 triaxpotent.o numeric_kinds_f.o iniparam_f.o
	$(fortran90) dmpotent.f90

interpolpotent.o : interpolpotent.f90 triaxpotent.o numeric_kinds_f.o iniparam_f.o
	$(fortran90) interpolpotent.f90

#####################################################
#
# AUX routines
#
###################################################

ran1_nr.o : ./ran1_nr.f
	$(fortran77) ./ran1_nr.f

numeric_kinds_f.o : ./sub/numeric_kinds_f.f90
	$(fortran90) ./sub/numeric_kinds_f.f90

ellipint.o : ./sub/ellipint.f90
	$(fortran90) ./sub/ellipint.f90

numrep.o: ./sub/numrep.f90 numeric_kinds_f.o
	$(fortran77) ./sub/numrep.f90

nag.o : sub/nag.f
	$(fortran77) sub/nag.f

numrec_arloc.o : sub/numrec_arloc.f
	$(fortran77) sub/numrec_arloc.f

dop853.o: ./sub/dop853.f
	$(fortran77) ./sub/dop853.f

dopri5.o: ./sub/dopri5.f
	$(fortran77) ./sub/dopri5.f

nnls95.o : ./sub/nnls95.f
	$(fortran77) ./sub/nnls95.f

gausherm.o : ./sub/gausherm.f
	$(fortran77) ./sub/gausherm.f

dqxgs.o : ./sub/dqxgs.f
	$(fortran77)  ./sub/dqxgs.f

specfunc_beta.o : ./sub/specfunc_beta.f90
	$(fortran90)  ./sub/specfunc_beta.f90

####################################################
#
# triaxmass
#
###################################################

triaxmass: $(MGEOBJ) nag.o triaxmass_f.o triaxmass.f90
	$(link) triaxmass $(MGEOBJ) triaxmass_f.o nag.o triaxmass.f90

triaxmass_bar: $(MGEOBJ) nag.o  triaxmass_f.o triaxmass_bar.f90
	$(link) triaxmass_bar $(MGEOBJ) triaxmass_f.o nag.o triaxmass_bar.f90

####################################################
#
# triaxmassbin
#
###################################################

triaxmassbin: $(MGEOBJ) triaxmassbin_f.o triaxmassbin.f90
	$(link) triaxmassbin $(MGEOBJ) triaxmassbin_f.o triaxmassbin.f90

triaxmassbin_bar: $(MGEOBJ) triaxmassbin_f.o triaxmassbin_bar.f90
	$(link) triaxmassbin_bar $(MGEOBJ) triaxmassbin_f.o triaxmassbin_bar.f90

######################################################
#
# orblib_f, The orbit library program.
#
#######################################################

ORBLIBNMCDE= $(MGECDE) sub/numrec_arloc.f $(DOPCDE) orblib_f_new_mirror.f90
ORBLIBNMOBJ= $(MGEOBJ) numrec_arloc.o $(DOPOBJ) orblib_f_new_mirror.o

ifeq ($(FAST),TRUE)
  ORBLIBNMOBJ=$(ORBLIBNMCDE)
endif

orblib_f_new_mirror.o : orblib_f_new_mirror.f90 $(MGEOBJ) numrec_arloc.o $(DOPOBJ)
	$(fortran90) orblib_f_new_mirror.f90

orblib_new_mirror : $(ORBLIBNMOBJ) orblibprogram.f90
	$(link) orblib_new_mirror  $(ORBLIBNMOBJ) orblibprogram.f90

orblib_bar : $(ORBLIBNMOBJ) orblibprogram_bar.f90
	$(link) orblib_bar  $(ORBLIBNMOBJ) orblibprogram_bar.f90

orbitstart_f.o : $(ORBLIBNMOBJ) $(MGEOBJ) orbitstart_f.f90
	$(fortran90) orbitstart_f.f90

orbitstart : $(ORBLIBNMOBJ) orbitstart_f.o orbitstart.f90
	$(link) orbitstart  $(ORBLIBNMOBJ) orbitstart_f.o orbitstart.f90

orbitstart_bar : $(ORBLIBNMOBJ) orbitstart_f.o orbitstart_bar.f90
	$(link) orbitstart_bar  $(ORBLIBNMOBJ) orbitstart_f.o orbitstart_bar.f90

##############################################################
#
# triaxnnls_CRcut
# triaxnnls_noCRcut
##############################################################

OBJECTS5_out = nnls95.o iniparam_f.o numeric_kinds_f.o gausherm.o triaxnnls_CRcut.o


triaxnnls_CRcut : $(OBJECTS5_out)
	$(link) triaxnnls_CRcut $(OBJECTS5_out)   -I$(GALAHADDIR)/modules/$(GALAHADTYPE)  -L$(GALAHADDIR)/objects/$(GALAHADTYPE)  -lgalahad  -lgalahad_hsl -lgalahad_metis -lgalahad_lapack -lgalahad_blas

triaxnnls_CRcut.o : triaxnnls_CRcut.f90 numeric_kinds_f.o iniparam_f.o
	$(fortran90) triaxnnls_CRcut.f90  -I$(GALAHADDIR)/modules/$(GALAHADTYPE)  -L$(GALAHADDIR)/objects/$(GALAHADTYPE)

##############################################################
#
##############################################################

OBJECTS5_out_nl = nnls95.o iniparam_f.o numeric_kinds_f.o gausherm.o triaxnnls_bar.o

triaxnnls_bar : $(OBJECTS5_out_nl)
	$(link) triaxnnls_bar $(OBJECTS5_out_nl)   -I$(GALAHADDIR)/modules/$(GALAHADTYPE)  -L$(GALAHADDIR)/objects/$(GALAHADTYPE)  -lgalahad  -lgalahad_hsl -lgalahad_metis -lgalahad_lapack -lgalahad_blas

triaxnnls_bar.o : triaxnnls_bar.f90 numeric_kinds_f.o iniparam_f.o
	$(fortran90) triaxnnls_bar.f90  -I$(GALAHADDIR)/modules/$(GALAHADTYPE)  -L$(GALAHADDIR)/objects/$(GALAHADTYPE)


##############################################################
#
##############################################################

OBJECTS7 = nnls95.o iniparam_f.o numeric_kinds_f.o gausherm.o triaxnnls_noCRcut.o

triaxnnls_noCRcut : $(OBJECTS7)
	$(link) triaxnnls_noCRcut $(OBJECTS7)   -I$(GALAHADDIR)/modules/$(GALAHADTYPE)  -L$(GALAHADDIR)/objects/$(GALAHADTYPE)  -lgalahad  -lgalahad_hsl -lgalahad_metis -lgalahad_lapack -lgalahad_blas

triaxnnls_noCRcut.o : triaxnnls_noCRcut.f90 numeric_kinds_f.o iniparam_f.o
	$(fortran90) triaxnnls_noCRcut.f90  -I$(GALAHADDIR)/modules/$(GALAHADTYPE)  -L$(GALAHADDIR)/objects/$(GALAHADTYPE)


##########################################
#
#  Utilities
#
##########################################

distclean:
	rm -f *.o *.il *.mod *.dyn work.* *.d *.exe ifc* *.dpi *.spi gmon.out orbitstart orbitstart_bar orblib_new_mirror orblib_bar triaxmass triaxmass_bar triaxmassbin triaxmassbin_bar triaxnnls_CRcut triaxnnls_noCRcut triaxnnls_bar

clean:
	rm -f *.o *.il *.mod work.* *.d ifc*

profclean:
	rm -f *.o
