CPPCMD = /usr/bin/cpp
CPPINCLUDEDIR =
CPPOTHERFLAGS = -P -C -nostdinc -nostdinc++
CPP = $(CPPCMD) -I$(CPPINCLUDEDIR) $(CPPOTHERFLAGS)
UMAKE = ${CUTER}/bin/umake.sh
CMP = GNU g77
MACHINE = Intel-like-PC
MCH = pc
COMPILE = g77 -c
LOAD = g77
ISF9095 = no
COMPILE9095 =
LOAD9095 =
FFLAGS = -O -fno-second-underscore -fPIC
SPECIALFLAGS =
FDEBUGFLAGS = FDebugFlags
FPROFILEFLAGS = FProfileFlags
NBYTES = 8
CCOMPILE = gcc -c -std=c99
CLOAD = gcc
CFLAGS = -O -fPIC
CXXFLAGS = -O -fPIC
SPECIALLIBS = -ldl -m64
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
UMAKE_CONFIGS = -C$(MYCUTER)/config
UMAKE_CPP_FLAGS =
IRULESRC = $(MYCUTER)/config
UMAKE_CMD = $(UMAKE) -I. -I$(IRULESRC) $(UMAKE_DEFINES) $(UMAKE_CONFIGS) $(UMAKE_CPP_FLAGS)
TMP = /tmp
SYSTEM = Linux
KEYSYS = CUNIX
MKDIR = mkdir
CP = cp
RM = rm -f
MV = mv
CAT = cat
CHMOD = chmod
SED = sed
LN = ln
LS = ls
AR = ar
ARFLAGS= -cr
RMDIR = rmdir
GREP = grep
AWK = gawk
HEAD = head
TAIL = tail
WC = wc
MAKE = make
LD = ld
MEXFORTRAN = mex
MEXFFLAGS = -fortran
INSTALL = install
COMMANDLIST = CMP MACHINE MCH COMPILE LOAD COMPILE9095 LOAD9095 FFLAGS FDEBUGFLAGS FPROFILEFLAGS NBYTES CCOMPILE CLOAD CFLAGS SPECIALLIBS R1 R2 R3 R4 R5 D1 D2 D3 D4 D5 PRECISION PRC PRC90 SIZE SIZ SIZ90 SYSTEM KEYSYS MKDIR CP RM MV CAT CHMOD SED LN LS AR ARFLAGS RMDIR GREP AWK HEAD TAIL WC MAKE LD MEXFORTRAN MEXFFLAGS
.SUFFIXES: .pro .f .o .f90
.f.o:
	true
.f90.o:
	true
DEPLIBS =
LDPOSTLIB =
DESTDIR =
CUTERSRCDIR = $(CUTER)/common/src/tools
CUTERDESTDIR = $(MYCUTER)/$(PRECISION)/bin
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
MAINPROGSOBJ = cobma.o dfoma.o 	filma.o 	genma.o ipoptma.o 	hrbma.o la04ma.o 	lmbma.o 	mnsma.o 	npsma.o 	oslma.o 	pdsma.o 	prxma.o 	snpma.o 	stnma.o 	tenma.o 	uncma.o 	cgpma.o 	cgdma.o 	lbsma.o 	lbbma.o 	trnma.o 	va15ma.o 	ve09ma.o 	ve14ma.o 	vf13ma.o 	nitsma.o
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

clean:: 
	$(RM) ${CMAINPROGSOBJ:.o=.c}
all:: $(LINPACOBJ) $(MAINPROGSOBJ) $(CMAINPROGSOBJ) $(CMAINHELPERS)
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
GENCSRC = $(GENSRCDIR)/genc.c
GENCOBJ = genc.o
GENCMA = $(TOOLSSRCDIR)/gencma.c
GENCMO = gencma.o

$(GENCOBJ): $(GENSRCDIR)/genc.c $(CUTER)/common/include/cuter.h 
	$(RM) $@ ${@:.o=.c} 
	@printf ' %-12s %-15s\t\t' "Building" $@ 
	$(CCOMPILE) $(CFLAGS) -I$(CUTER)/common/include -o $@ $(GENSRCDIR)/genc.c 
	@printf '[Ok]\n'

$(GENCMO): $(TOOLSSRCDIR)/gencma.c $(CUTER)/common/include/cuter.h 
	$(RM) $@ ${@:.o=.c} 
	@printf ' %-12s %-15s\t\t' "Building" $@ 
	$(CCOMPILE) $(CFLAGS) -I$(CUTER)/common/include -DIsg77 -o $@ $(TOOLSSRCDIR)/gencma.c 
	@printf '[Ok]\n'
all:: ${GENCOBJ}
all:: ${GENCMO}

purge:: 
	$(RM) ${GENCOBJ}

purge:: 
	$(RM) ${GENCMO}
DERCHKMA = $(TOOLSSRCDIR)/derchkma.c
DERCHKMO = derchkma.o

$(DERCHKMO): $(TOOLSSRCDIR)/derchkma.c $(CUTER)/common/include/cuter.h 
	$(RM) $@ ${@:.o=.c} 
	@printf ' %-12s %-15s\t\t' "Building" $@ 
	$(CCOMPILE) $(CFLAGS) -I$(CUTER)/common/include -DIsg77 -o $@ $(TOOLSSRCDIR)/derchkma.c 
	@printf '[Ok]\n'
all:: ${DERCHKMO}

purge:: 
	$(RM) ${DERCHKMO}

MFLAGS =
emptyrule::

clean:: 
	$(RM) *.CKP *.ln *.BAK *.bak core errs ,* *~ tags TAGS make.log MakeOut "#"*

Makefile:: 
	-@if [ -f Makefile ]; then set -x; \
	$(RM) Makefile.bak; $(MV) Makefile Makefile.bak; \
	else exit 0; fi 
	$(UMAKE_CMD) -DTOPDIR=$(TOP) -DCURDIR=$(CURRENT_DIR) -Dlinux -DIsg77 -DMediumSize -DDoublePrecision -DIsgcc -UNoCcomp
Makefiles::
includes::
