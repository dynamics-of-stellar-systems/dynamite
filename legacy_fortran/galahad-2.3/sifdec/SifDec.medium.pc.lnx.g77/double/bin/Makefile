CPPCMD = /usr/bin/cpp
CPPINCLUDEDIR = 
CPPOTHERFLAGS = -P -C -nostdinc -nostdinc++
CPP = $(CPP) -I$(CPPINCLUDEDIR) $(CPPOTHERFLAGS)

UMAKE = ${SIFDEC}/build/scripts/umake

CMP = GNU g77 
MACHINE = Intel-like-PC 
MCH = pc 
COMPILE = g77 -c 
LOAD = g77 
ISF9095 = no 
COMPILE9095 = 
LOAD9095 = 
FFLAGS = -O -fno-second-underscore 
SPECIALFLAGS = 
NBYTES = 8 
CCOMPILE = 
CLOAD = CLoadCmd
CFLAGS = CFlags
R1 = 1.1920930E-07
R2 = 5.9604646E-08
R3 = 1.1754945E-38
R4 = 1.1754945E-38
R5 = 3.4028234E+38
D1 = 2.2204460492503132D-16
D2 = 1.1102230246251566D-16
D3 = 2.225073858507202D-308
D4 = 2.225073858507202D-308
D5 = 1.797693134862314D+308

PRECISION = double
PRC = CD
PRC90 = \!D

SIZE = medium
SIZ = CMED
SIZ90 = \!MED

PATHSEP = /
SHELL = /bin/sh

TOP = ../..
CURRENT_DIR = /double/bin

UMAKE_DEFINES = 
UMAKE_CONFIGS = -C$(MYSIFDEC)/config 

UMAKE_CPP_FLAGS = -Ptraditional -Pw

IRULESRC = $(MYSIFDEC)/config 
UMAKE_CMD = $(UMAKE) -I. -I$(IRULESRC) $(UMAKE_DEFINES) $(UMAKE_CONFIGS) $(UMAKE_CPP_FLAGS)
ICONFIGFILES = $(IRULESRC)/Umake.tmpl $(IRULESRC)/site.def $(IRULESRC)/$(MACROFILE) SiteIConfigFiles $(EXTRA_ICONFIGFILES)

TMP = /tmp

SYSTEM = Linux
KEYSYS = CUNIX

MKDIR = /bin/mkdir
CP = /bin/cp
RM = /bin/rm -f
MV = /bin/mv
CAT = /bin/cat
CHMOD = /bin/chmod
SED = /bin/sed
LN = /bin/ln
LS = /bin/ls
AR = /usr/bin/ar
ARFLAGS= -cr
RMDIR = /bin/rmdir
GREP = /bin/grep
AWK = /usr/bin/gawk
HEAD = /usr/bin/head
TAIL = /usr/bin/tail
WC = /usr/bin/wc
MAKE = /usr/bin/make

MEXFORTRAN = /usr/bin/mex 
MEXFFLAGS = -fortran

INSTALL = install

COMMANDLIST = CMP MACHINE MCH COMPILE LOAD COMPILE9095 LOAD9095 FFLAGS NBYTES CCOMPILE CLOAD CFLAGS R1 R2 R3 R4 R5 D1 D2 D3 D4 D5 PRECISION PRC PRC90 SIZE SIZ SIZ90 SYSTEM KEYSYS MKDIR CP RM MV CAT CHMOD SED LN LS AR ARFLAGS RMDIR GREP AWK HEAD TAIL WC MAKE MEXFORTRAN MEXFFLAGS

.SUFFIXES: .pro .f .o .f90

.f.o:
	true

.f90.o:
	true

DEPLIBS =
LDPOSTLIB =

DESTDIR =	
SIFDECSRCDIR = $(SIFDEC)/common/src/sifdec
SIFDECDESTDIR = $(MYSIFDEC)/$(PRECISION)/bin
INTCONFIG = $(MYSIFDEC)/$(PRECISION)/config

all::

baseconfig:
	cd $(MYSIFDEC)/$(PRECISION)/config; 	$(MAKE) $(MFLAGS) all;

all:: baseconfig

SIFDECOBJ = sifdec.o

SD_OTHER_OBJS = decode.o 	gps.o 		inlanc.o 		makefn.o 		mafnad.o 		makegr.o 		magrad.o 		printp.o 		rest.o 		trans.o 		utils.o

SD_OTHER_SRCS = $(SIFDECSRCDIR)/decode.f 	$(SIFDECSRCDIR)/gps.f 		$(SIFDECSRCDIR)/inlanc.f 		$(SIFDECSRCDIR)/makefn.f 		$(SIFDECSRCDIR)/mafnad.f 		$(SIFDECSRCDIR)/makegr.f 		$(SIFDECSRCDIR)/magrad.f 		$(SIFDECSRCDIR)/printp.f 		$(SIFDECSRCDIR)/rest.f 		$(SIFDECSRCDIR)/trans.f 		$(SIFDECSRCDIR)/utils.f

$(SIFDECOBJ) : $(SIFDECSRCDIR)/sifdec.f $(INTCONFIG)/sifdec.siz $(INTCONFIG)/cast.sed 
	$(RM) $@ ${@:.o=.f} 
	@printf ' %-12s %-15s\t\t' "Building" $@ 
	$(SED) "/#{sizing}/r $(INTCONFIG)/sifdec.siz" $(SIFDECSRCDIR)/${@:.o=.f} | $(SED) -f $(INTCONFIG)/cast.sed > ${@:.o=.f} 
	$(COMPILE) $(FFLAGS) -o $@ ${@:.o=.f} 
	@printf '[Ok]\n'

clean:: 
	$(RM) ${SIFDECOBJ:.o=.f}

${SD_OTHER_OBJS}: $(SD_OTHER_SRCS) 
	@printf ' %-12s %-15s\t\t' "Building" $@ 
	$(COMPILE) $(FFLAGS) -o $@ $(SIFDECSRCDIR)/${@:.o=.f} 
	@printf '[Ok]\n'

purge:: 
	$(RM) $(SIFDECOBJ) $(SD_OTHER_OBJS)

LOCAL_OBJ = local.o

OBJECTS = $(SIFDECOBJ) $(SD_OTHER_OBJS) $(INTCONFIG)/$(LOCAL_OBJ)

sifdec: $(OBJECTS) 
	$(RM) sifdec 
	@printf ' %-12s %-15s\t\t' "Linking" sifdec 
	$(LOAD) $(FFLAGS) -o $@ $(F_LDOPTIONS) $(OBJECTS) $(LDLIBS) $(EXTRA_LOAD_FLAGS) 
	@printf '[Ok]\n' 
all:: sifdec 

purge:: 
	$(RM) sifdec

emptyrule::

clean:: 
	$(RM) FilesToClean ExtraFilesToClean "#"*

Makefile:: 
	-@if [ -f Makefile ]; then set -x; \
	$(RM) Makefile.bak; $(MV) Makefile Makefile.bak; \
	else exit 0; fi 
	$(UMAKE_CMD) -DTOPDIR=$(TOP) -DCURDIR=$(CURRENT_DIR) -Dlinux -DIsg77 -DMediumSize -DDoublePrecision -UNoCcomp -Ptraditional -Pw

Makefiles::

includes::

