CPPCMD		 = /usr/bin/cpp
CPPINCLUDEDIR	 = 
CPPOTHERFLAGS	 = -P -C -nostdinc -nostdinc++
CPP		 = $(CPP) -I$(CPPINCLUDEDIR) $(CPPOTHERFLAGS)

UMAKE		 = ${CUTER}/build/scripts/umake

CMP		 = gfo	
MACHINE		 = Intel-like-PC	
MCH		 = pc	
COMPILE		 = gfortran -c -ffixed-form -fno-second-underscore	
LOAD		 = gfortran -ffixed-form -fno-second-underscore		
ISF9095		 = yes	
COMPILE9095	 = gfortran -c -fno-second-underscore	
LOAD9095	 = gfortran -fno-second-underscore	
FFLAGS		 = -O	
SPECIALFLAGS	 = 	
NBYTES		 = 8	
CCOMPILE	 = 
CLOAD		 = CLoadCmd
CFLAGS		 = CFlags
SPECIALLIBS = SpecialLibs
R1		 = 1.1920930E-07
R2		 = 5.9604646E-08
R3		 = 1.1754945E-38
R4		 = 1.1754945E-38
R5		 = 3.4028234E+38
D1		 = 2.2204460492503132D-16
D2		 = 1.1102230246251566D-16
D3		 = 2.225073858507202D-308
D4		 = 2.225073858507202D-308
D5		 = 1.797693134862314D+308

PRECISION = double
PRC = CD
PRC90 = \!D

SIZE = large
SIZ = CBIG
SIZ90 = \!BIG

PATHSEP = /
SHELL = /bin/sh

TOP = ../..
CURRENT_DIR = /double/bin

UMAKE_DEFINES = 
UMAKE_CONFIGS = -C$(MYCUTER)/config 

UMAKE_CPP_FLAGS = -Ptraditional -Pw

IRULESRC = $(MYCUTER)/config	
UMAKE_CMD = $(UMAKE) -I. -I$(IRULESRC) $(UMAKE_DEFINES) $(UMAKE_CONFIGS) $(UMAKE_CPP_FLAGS)

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
LD = /usr/bin/ld

MEXFORTRAN = /usr/bin/mex 
MEXFFLAGS = -fortran

INSTALL = install

COMMANDLIST = CMP MACHINE MCH COMPILE LOAD COMPILE9095 LOAD9095 FFLAGS NBYTES CCOMPILE CLOAD CFLAGS SPECIALLIBS R1 R2 R3 R4 R5 D1 D2 D3 D4 D5 PRECISION PRC PRC90 SIZE SIZ SIZ90 SYSTEM KEYSYS MKDIR CP RM MV CAT CHMOD SED LN LS AR ARFLAGS RMDIR GREP AWK HEAD TAIL WC MAKE LD MEXFORTRAN MEXFFLAGS

.SUFFIXES: .pro .f .o .f90

.f.o:
	true

.f90.o:
	true

DEPLIBS =
LDPOSTLIB =

DESTDIR =	
SIFDECSRCDIR = $(CUTER)/common/src/sifdec
SIFDECDESTDIR = $(MYCUTER)/$(PRECISION)/bin
INTCONFIG = $(MYCUTER)/$(PRECISION)/config
INTERFACESSRC = $(CUTER)/build/prototypes
TOOLSSRCDIR = $(CUTER)/common/src/tools
PKGSRCDIR = $(CUTER)/common/src/pkg
TOOLSBINDIR = $(MYCUTER)/$(PRECISION)/bin

all::

baseconfig:
	cd $(MYCUTER)/$(PRECISION)/config; 	$(MAKE) $(MFLAGS) all;

all:: baseconfig

LOCAL_OBJ = local.o
LINPACOBJ = linpac.o

${LINPACOBJ}: $(TOOLSSRCDIR)/linpac.f $(INTCONFIG)/cast.sed 
	$(RM) $@ ${@:.o=.f} 
	@printf ' %-12s %-15s\t\t' "Building" $@ 
	$(SED) -f $(INTCONFIG)/cast.sed $(TOOLSSRCDIR)/${@:.o=.f} > ${@:.o=.f} 
	$(COMPILE) $(FFLAGS) -o $@ ${@:.o=.f} 
	@printf '[Ok]\n'

clean:: 
	$(RM) ${LINPACOBJ:.o=.f}

purge:: 
	$(RM) $(LINPACOBJ)

MAINPROGSOBJ = cobma.o		filma.o		genma.o	 ipoptma.o 	hrbma.o	 la04ma.o 	lmbma.o		mnsma.o		npsma.o		oslma.o		pdsma.o		prxma.o		snpma.o		stnma.o		tenma.o		uncma.o		cgpma.o		cgdma.o		lbsma.o		lbbma.o		trnma.o		va15ma.o		ve09ma.o		ve14ma.o		vf13ma.o 	nitsma.o

F_KNITRO =

C_KNITRO = 
C_LOQO =

CMAINPROGSOBJ = $(C_KNITRO) $(C_LOQO)

cobma.o: $(TOOLSSRCDIR)/cobma.f $(INTCONFIG)/tools.siz $(INTCONFIG)/cast.sed 
	$(RM) $@ ${@:.o=.f} 
	@printf ' %-12s %-15s\t\t' "Building" $@ 
	$(SED) "/#{sizing}/r $(INTCONFIG)/tools.siz" $(TOOLSSRCDIR)/cobma.f | $(SED) -f $(INTCONFIG)/cast.sed > ${@:.o=.f} 
	$(COMPILE) $(FFLAGS) -o $@ ${@:.o=.f} 
	@printf '[Ok]\n'

filma.o: $(TOOLSSRCDIR)/filma.f $(INTCONFIG)/tools.siz $(INTCONFIG)/cast.sed 
	$(RM) $@ ${@:.o=.f} 
	@printf ' %-12s %-15s\t\t' "Building" $@ 
	$(SED) "/#{sizing}/r $(INTCONFIG)/tools.siz" $(TOOLSSRCDIR)/filma.f | $(SED) -f $(INTCONFIG)/cast.sed > ${@:.o=.f} 
	$(COMPILE) $(FFLAGS) -o $@ ${@:.o=.f} 
	@printf '[Ok]\n'

genma.o: $(TOOLSSRCDIR)/genma.f $(INTCONFIG)/tools.siz $(INTCONFIG)/cast.sed 
	$(RM) $@ ${@:.o=.f} 
	@printf ' %-12s %-15s\t\t' "Building" $@ 
	$(SED) "/#{sizing}/r $(INTCONFIG)/tools.siz" $(TOOLSSRCDIR)/genma.f | $(SED) -f $(INTCONFIG)/cast.sed > ${@:.o=.f} 
	$(COMPILE) $(FFLAGS) -o $@ ${@:.o=.f} 
	@printf '[Ok]\n'

hrbma.o: $(TOOLSSRCDIR)/hrbma.f $(INTCONFIG)/tools.siz $(INTCONFIG)/cast.sed 
	$(RM) $@ ${@:.o=.f} 
	@printf ' %-12s %-15s\t\t' "Building" $@ 
	$(SED) "/#{sizing}/r $(INTCONFIG)/tools.siz" $(TOOLSSRCDIR)/hrbma.f | $(SED) -f $(INTCONFIG)/cast.sed > ${@:.o=.f} 
	$(COMPILE) $(FFLAGS) -o $@ ${@:.o=.f} 
	@printf '[Ok]\n'

ipoptma.o: $(TOOLSSRCDIR)/ipoptma.f $(INTCONFIG)/tools.siz $(INTCONFIG)/cast.sed $(INTCONFIG)/ipoptma.inc $(INTCONFIG)/ipoptma.wrk 
	$(RM) $@ ${@:.o=.f} 
	@printf ' %-12s %-15s\t\t' "Building" $@ 
	$(SED) "/#{sizing}/r $(INTCONFIG)/tools.siz" $(TOOLSSRCDIR)/ipoptma.f | $(SED) -f $(INTCONFIG)/cast.sed > ${@:.o=.f} 
	$(COMPILE) $(FFLAGS) -o $@ ${@:.o=.f} 
	@printf '[Ok]\n'

la04ma.o: $(TOOLSSRCDIR)/la04ma.f $(INTCONFIG)/tools.siz $(INTCONFIG)/cast.sed 
	$(RM) $@ ${@:.o=.f} 
	@printf ' %-12s %-15s\t\t' "Building" $@ 
	$(SED) "/#{sizing}/r $(INTCONFIG)/tools.siz" $(TOOLSSRCDIR)/la04ma.f | $(SED) -f $(INTCONFIG)/cast.sed > ${@:.o=.f} 
	$(COMPILE) $(FFLAGS) -o $@ ${@:.o=.f} 
	@printf '[Ok]\n'

lmbma.o: $(TOOLSSRCDIR)/lmbma.f $(INTCONFIG)/tools.siz $(INTCONFIG)/cast.sed 
	$(RM) $@ ${@:.o=.f} 
	@printf ' %-12s %-15s\t\t' "Building" $@ 
	$(SED) "/#{sizing}/r $(INTCONFIG)/tools.siz" $(TOOLSSRCDIR)/lmbma.f | $(SED) -f $(INTCONFIG)/cast.sed > ${@:.o=.f} 
	$(COMPILE) $(FFLAGS) -o $@ ${@:.o=.f} 
	@printf '[Ok]\n'

mnsma.o: $(TOOLSSRCDIR)/mnsma.f $(INTCONFIG)/tools.siz $(INTCONFIG)/cast.sed 
	$(RM) $@ ${@:.o=.f} 
	@printf ' %-12s %-15s\t\t' "Building" $@ 
	$(SED) "/#{sizing}/r $(INTCONFIG)/tools.siz" $(TOOLSSRCDIR)/mnsma.f | $(SED) -f $(INTCONFIG)/cast.sed > ${@:.o=.f} 
	$(COMPILE) $(FFLAGS) -o $@ ${@:.o=.f} 
	@printf '[Ok]\n'

npsma.o: $(TOOLSSRCDIR)/npsma.f $(INTCONFIG)/tools.siz $(INTCONFIG)/cast.sed 
	$(RM) $@ ${@:.o=.f} 
	@printf ' %-12s %-15s\t\t' "Building" $@ 
	$(SED) "/#{sizing}/r $(INTCONFIG)/tools.siz" $(TOOLSSRCDIR)/npsma.f | $(SED) -f $(INTCONFIG)/cast.sed > ${@:.o=.f} 
	$(COMPILE) $(FFLAGS) -o $@ ${@:.o=.f} 
	@printf '[Ok]\n'

oslma.o: $(TOOLSSRCDIR)/oslma.f $(INTCONFIG)/tools.siz $(INTCONFIG)/cast.sed 
	$(RM) $@ ${@:.o=.f} 
	@printf ' %-12s %-15s\t\t' "Building" $@ 
	$(SED) "/#{sizing}/r $(INTCONFIG)/tools.siz" $(TOOLSSRCDIR)/oslma.f | $(SED) -f $(INTCONFIG)/cast.sed > ${@:.o=.f} 
	$(COMPILE) $(FFLAGS) -o $@ ${@:.o=.f} 
	@printf '[Ok]\n'

pdsma.o: $(TOOLSSRCDIR)/pdsma.f $(INTCONFIG)/tools.siz $(INTCONFIG)/cast.sed 
	$(RM) $@ ${@:.o=.f} 
	@printf ' %-12s %-15s\t\t' "Building" $@ 
	$(SED) "/#{sizing}/r $(INTCONFIG)/tools.siz" $(TOOLSSRCDIR)/pdsma.f | $(SED) -f $(INTCONFIG)/cast.sed > ${@:.o=.f} 
	$(COMPILE) $(FFLAGS) -o $@ ${@:.o=.f} 
	@printf '[Ok]\n'

prxma.o: $(TOOLSSRCDIR)/prxma.f $(INTCONFIG)/tools.siz $(INTCONFIG)/cast.sed 
	$(RM) $@ ${@:.o=.f} 
	@printf ' %-12s %-15s\t\t' "Building" $@ 
	$(SED) "/#{sizing}/r $(INTCONFIG)/tools.siz" $(TOOLSSRCDIR)/prxma.f | $(SED) -f $(INTCONFIG)/cast.sed > ${@:.o=.f} 
	$(COMPILE) $(FFLAGS) -o $@ ${@:.o=.f} 
	@printf '[Ok]\n'

snpma.o: $(TOOLSSRCDIR)/snpma.f $(INTCONFIG)/tools.siz $(INTCONFIG)/cast.sed 
	$(RM) $@ ${@:.o=.f} 
	@printf ' %-12s %-15s\t\t' "Building" $@ 
	$(SED) "/#{sizing}/r $(INTCONFIG)/tools.siz" $(TOOLSSRCDIR)/snpma.f | $(SED) -f $(INTCONFIG)/cast.sed > ${@:.o=.f} 
	$(COMPILE) $(FFLAGS) -o $@ ${@:.o=.f} 
	@printf '[Ok]\n'

stnma.o: $(TOOLSSRCDIR)/stnma.f $(INTCONFIG)/tools.siz $(INTCONFIG)/cast.sed 
	$(RM) $@ ${@:.o=.f} 
	@printf ' %-12s %-15s\t\t' "Building" $@ 
	$(SED) "/#{sizing}/r $(INTCONFIG)/tools.siz" $(TOOLSSRCDIR)/stnma.f | $(SED) -f $(INTCONFIG)/cast.sed > ${@:.o=.f} 
	$(COMPILE) $(FFLAGS) -o $@ ${@:.o=.f} 
	@printf '[Ok]\n'

tenma.o: $(TOOLSSRCDIR)/tenma.f $(INTCONFIG)/tools.siz $(INTCONFIG)/cast.sed 
	$(RM) $@ ${@:.o=.f} 
	@printf ' %-12s %-15s\t\t' "Building" $@ 
	$(SED) "/#{sizing}/r $(INTCONFIG)/tools.siz" $(TOOLSSRCDIR)/tenma.f | $(SED) -f $(INTCONFIG)/cast.sed > ${@:.o=.f} 
	$(COMPILE) $(FFLAGS) -o $@ ${@:.o=.f} 
	@printf '[Ok]\n'

cgpma.o: $(TOOLSSRCDIR)/cgpma.f $(INTCONFIG)/tools.siz $(INTCONFIG)/cast.sed 
	$(RM) $@ ${@:.o=.f} 
	@printf ' %-12s %-15s\t\t' "Building" $@ 
	$(SED) "/#{sizing}/r $(INTCONFIG)/tools.siz" $(TOOLSSRCDIR)/cgpma.f | $(SED) -f $(INTCONFIG)/cast.sed > ${@:.o=.f} 
	$(COMPILE) $(FFLAGS) -o $@ ${@:.o=.f} 
	@printf '[Ok]\n'

cgdma.o: $(TOOLSSRCDIR)/cgdma.f $(INTCONFIG)/tools.siz $(INTCONFIG)/cast.sed 
	$(RM) $@ ${@:.o=.f} 
	@printf ' %-12s %-15s\t\t' "Building" $@ 
	$(SED) "/#{sizing}/r $(INTCONFIG)/tools.siz" $(TOOLSSRCDIR)/cgdma.f | $(SED) -f $(INTCONFIG)/cast.sed > ${@:.o=.f} 
	$(COMPILE) $(FFLAGS) -o $@ ${@:.o=.f} 
	@printf '[Ok]\n'

lbsma.o: $(TOOLSSRCDIR)/lbsma.f $(INTCONFIG)/tools.siz $(INTCONFIG)/cast.sed 
	$(RM) $@ ${@:.o=.f} 
	@printf ' %-12s %-15s\t\t' "Building" $@ 
	$(SED) "/#{sizing}/r $(INTCONFIG)/tools.siz" $(TOOLSSRCDIR)/lbsma.f | $(SED) -f $(INTCONFIG)/cast.sed > ${@:.o=.f} 
	$(COMPILE) $(FFLAGS) -o $@ ${@:.o=.f} 
	@printf '[Ok]\n'

lbbma.o: $(TOOLSSRCDIR)/lbbma.f $(INTCONFIG)/tools.siz $(INTCONFIG)/cast.sed 
	$(RM) $@ ${@:.o=.f} 
	@printf ' %-12s %-15s\t\t' "Building" $@ 
	$(SED) "/#{sizing}/r $(INTCONFIG)/tools.siz" $(TOOLSSRCDIR)/lbbma.f | $(SED) -f $(INTCONFIG)/cast.sed > ${@:.o=.f} 
	$(COMPILE) $(FFLAGS) -o $@ ${@:.o=.f} 
	@printf '[Ok]\n'

trnma.o: $(TOOLSSRCDIR)/trnma.f $(INTCONFIG)/tools.siz $(INTCONFIG)/cast.sed 
	$(RM) $@ ${@:.o=.f} 
	@printf ' %-12s %-15s\t\t' "Building" $@ 
	$(SED) "/#{sizing}/r $(INTCONFIG)/tools.siz" $(TOOLSSRCDIR)/trnma.f | $(SED) -f $(INTCONFIG)/cast.sed > ${@:.o=.f} 
	$(COMPILE) $(FFLAGS) -o $@ ${@:.o=.f} 
	@printf '[Ok]\n'

uncma.o: $(TOOLSSRCDIR)/uncma.f $(INTCONFIG)/tools.siz $(INTCONFIG)/cast.sed 
	$(RM) $@ ${@:.o=.f} 
	@printf ' %-12s %-15s\t\t' "Building" $@ 
	$(SED) "/#{sizing}/r $(INTCONFIG)/tools.siz" $(TOOLSSRCDIR)/uncma.f | $(SED) -f $(INTCONFIG)/cast.sed > ${@:.o=.f} 
	$(COMPILE) $(FFLAGS) -o $@ ${@:.o=.f} 
	@printf '[Ok]\n'

va15ma.o: $(TOOLSSRCDIR)/va15ma.f $(INTCONFIG)/tools.siz $(INTCONFIG)/cast.sed 
	$(RM) $@ ${@:.o=.f} 
	@printf ' %-12s %-15s\t\t' "Building" $@ 
	$(SED) "/#{sizing}/r $(INTCONFIG)/tools.siz" $(TOOLSSRCDIR)/va15ma.f | $(SED) -f $(INTCONFIG)/cast.sed > ${@:.o=.f} 
	$(COMPILE) $(FFLAGS) -o $@ ${@:.o=.f} 
	@printf '[Ok]\n'

ve09ma.o: $(TOOLSSRCDIR)/ve09ma.f $(INTCONFIG)/tools.siz $(INTCONFIG)/cast.sed 
	$(RM) $@ ${@:.o=.f} 
	@printf ' %-12s %-15s\t\t' "Building" $@ 
	$(SED) "/#{sizing}/r $(INTCONFIG)/tools.siz" $(TOOLSSRCDIR)/ve09ma.f | $(SED) -f $(INTCONFIG)/cast.sed > ${@:.o=.f} 
	$(COMPILE) $(FFLAGS) -o $@ ${@:.o=.f} 
	@printf '[Ok]\n'

ve14ma.o: $(TOOLSSRCDIR)/ve14ma.f $(INTCONFIG)/tools.siz $(INTCONFIG)/cast.sed 
	$(RM) $@ ${@:.o=.f} 
	@printf ' %-12s %-15s\t\t' "Building" $@ 
	$(SED) "/#{sizing}/r $(INTCONFIG)/tools.siz" $(TOOLSSRCDIR)/ve14ma.f | $(SED) -f $(INTCONFIG)/cast.sed > ${@:.o=.f} 
	$(COMPILE) $(FFLAGS) -o $@ ${@:.o=.f} 
	@printf '[Ok]\n'

vf13ma.o: $(TOOLSSRCDIR)/vf13ma.f $(INTCONFIG)/tools.siz $(INTCONFIG)/cast.sed 
	$(RM) $@ ${@:.o=.f} 
	@printf ' %-12s %-15s\t\t' "Building" $@ 
	$(SED) "/#{sizing}/r $(INTCONFIG)/tools.siz" $(TOOLSSRCDIR)/vf13ma.f | $(SED) -f $(INTCONFIG)/cast.sed > ${@:.o=.f} 
	$(COMPILE) $(FFLAGS) -o $@ ${@:.o=.f} 
	@printf '[Ok]\n'

nitsma.o: $(TOOLSSRCDIR)/nitsma.f $(INTCONFIG)/tools.siz $(INTCONFIG)/cast.sed 
	$(RM) $@ ${@:.o=.f} 
	@printf ' %-12s %-15s\t\t' "Building" $@ 
	$(SED) "/#{sizing}/r $(INTCONFIG)/tools.siz" $(TOOLSSRCDIR)/nitsma.f | $(SED) -f $(INTCONFIG)/cast.sed > ${@:.o=.f} 
	$(COMPILE) $(FFLAGS) -o $@ ${@:.o=.f} 
	@printf '[Ok]\n'

clean:: 
	$(RM) ${MAINPROGSOBJ:.o=.f}

purge:: 
	$(RM) $(MAINPROGSOBJ)

F9095TOOLSSRC = $(TOOLSSRCDIR)/cuterinter.f90 	$(TOOLSSRCDIR)/readin.f90

F9095TOOLSOBJ = cuterinter.o readin.o

${F9095TOOLSOBJ}: $(F9095TOOLSSRC) $(INTCONFIG)/cast90.sed 
	$(RM) $@ ${@:.o=.f90} 
	@printf ' %-12s %-15s\t\t' "Building" $@ 
	$(SED) -f $(INTCONFIG)/cast90.sed $(TOOLSSRCDIR)/${@:.o=.f90} > ${@:.o=.f90} 
	$(COMPILE9095) $(FFLAGS) -o $@ ${@:.o=.f90} 
	@printf '[Ok]\n'

clean:: 
	$(RM) ${F9095TOOLSOBJ:.o=.f90}

purge:: 
	$(RM) $(F9095TOOLSOBJ)

all:: $(LINPACOBJ) $(MAINPROGSOBJ) $(F_KNITRO) $(CMAINPROGSOBJ) $(CMAINHELPERS) $(F9095TOOLSOBJ)

GENSRCDIR = $(CUTER)/common/src/pkg/gen
GENSRC = $(GENSRCDIR)/gen.f
GENOBJ = gen.o

$(GENOBJ): $(GENSRCDIR)/gen.f $(INTCONFIG)/tools.siz $(INTCONFIG)/cast.sed 
	$(RM) $@ ${@:.o=.f} 
	@printf ' %-12s %-15s\t\t' "Building" $@ 
	$(SED) "/#{sizing}/r $(INTCONFIG)/tools.siz" $(GENSRCDIR)/gen.f | $(SED) -f $(INTCONFIG)/cast.sed > ${@:.o=.f} 
	$(COMPILE) $(FFLAGS) -o $@ ${@:.o=.f} 
	@printf '[Ok]\n'

purge:: 
	$(RM) $(GENOBJ)
all:: $(GENOBJ)

GEN90SRC = $(GENSRCDIR)/gen90.f90
GEN90OBJ = gen90.o

$(GEN90OBJ): $(GEN90SRC) $(INTCONFIG)/cast90.sed 
	$(RM) $@ ${@:.o=.f90} 
	@printf ' %-12s %-15s\t\t' "Building" $@ 
	$(SED) -f $(INTCONFIG)/cast90.sed $(GENSRCDIR)/${@:.o=.f90} > ${@:.o=.f90} 
	$(COMPILE9095) $(FFLAGS) -o $@ ${@:.o=.f90} 
	@printf '[Ok]\n'

GEN90MA = $(TOOLSSRCDIR)/gen90ma.f90
GEN90MO = gen90ma.o

${GEN90MO}: $(GEN90MA) $(INTCONFIG)/cast90.sed 
	$(RM) $@ ${@:.o=.f90} 
	@printf ' %-12s %-15s\t\t' "Building" $@ 
	$(SED) -f $(INTCONFIG)/cast90.sed $(TOOLSSRCDIR)/${@:.o=.f90} > ${@:.o=.f90} 
	$(COMPILE9095) $(FFLAGS) -o $@ ${@:.o=.f90} 
	@printf '[Ok]\n'

STATSMA = $(TOOLSSRCDIR)/statsma.f90
STATSMO = statsma.o

${STATSMO}: $(STATSMA) $(INTCONFIG)/cast90.sed 
	$(RM) $@ ${@:.o=.f90} 
	@printf ' %-12s %-15s\t\t' "Building" $@ 
	$(SED) -f $(INTCONFIG)/cast90.sed $(TOOLSSRCDIR)/${@:.o=.f90} > ${@:.o=.f90} 
	$(COMPILE9095) $(FFLAGS) -o $@ ${@:.o=.f90} 
	@printf '[Ok]\n'

all:: ${GEN90OBJ}
all:: ${GEN90MO}

clean:: 
	$(RM) ${GEN90OBJ:.o=.f90}

purge:: 
	$(RM) ${GEN90OBJ}

purge:: 
	$(RM) ${GEN90MO}

all:: ${STATSMO}

purge:: 
	$(RM) ${STATSMO}

MATLABSRCDIR = $(CUTER)/common/src/tools

LOCALMATLABSRCS = utools.f ctools.f gtools.f

utools.f: $(MATLABSRCDIR)/utools.ff 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/cast.sed $(MATLABSRCDIR)/utools.ff > $@ 
	@printf '[Ok]\n'

ctools.f: $(MATLABSRCDIR)/ctools.ff 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/cast.sed $(MATLABSRCDIR)/ctools.ff > $@ 
	@printf '[Ok]\n'

gtools.f: $(MATLABSRCDIR)/gtools.ff 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/cast.sed $(MATLABSRCDIR)/gtools.ff > $@ 
	@printf '[Ok]\n'

purge:: 
	$(RM) $(LOCALMATLABSRCS)
all:: $(LOCALMATLABSRCS)

MFLAGS =

emptyrule::

clean:: 
	$(RM) *.CKP *.ln *.BAK *.bak core errs ,* *~ tags TAGS make.log MakeOut "#"*

Makefile:: 
	-@if [ -f Makefile ]; then set -x; \
	$(RM) Makefile.bak; $(MV) Makefile Makefile.bak; \
	else exit 0; fi 
	$(UMAKE_CMD) -DTOPDIR=$(TOP) -DCURDIR=$(CURRENT_DIR) -Dlinux -DIsgfo -DLargeSize -DDoublePrecision -DIsa90Compiler -DNoCcomp -Ptraditional -Pw

Makefiles::

includes::

