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
TOP = ..
CURRENT_DIR = /bin
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
.SUFFIXES: .pro
INTERFACESSRC = $(CUTER)/build/prototypes
INTERFACESBIN = $(MYCUTER)/bin
INTCONFIG = $(MYCUTER)/$(PRECISION)/config
DESTDIR =
TARGETSHELL=sh
INTERFACES_AT_LARGE = cob cob90 derchk dfo 	envcheck 	fil 	gen 	gen90 	helpmsg 	hrb ipopt la04 	mns 	mx 	mx_pre 	nps 	osl 	pds 	prx 	nits runcuter 	snp 	stats 	stn 	ten 	unc 	cgp 	cgd 	lbs 	lbb 	trn 	va15 	ve09 	ve12 	ve14 	vf13
C_INTERFACES = derchk genc
LOQO_INTERFACES =
KNITRO_INTERFACES =
TAO_INTERFACES =
all::
baseconfig:
	cd $(MYCUTER)/$(PRECISION)/config; 	$(MAKE) $(MFLAGS) all;
all:: baseconfig
INTERFACES_DEPS = $(INTCONFIG)/cmds $(INTCONFIG)/this_build $(INTCONFIG)/script.sed

cob: $(INTERFACESSRC)/cob.$(TARGETSHELL).pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/cob.$(TARGETSHELL).pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

cob90: $(INTERFACESSRC)/cob90.$(TARGETSHELL).pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/cob90.$(TARGETSHELL).pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

dfo: $(INTERFACESSRC)/dfo.$(TARGETSHELL).pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/dfo.$(TARGETSHELL).pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

fil: $(INTERFACESSRC)/fil.$(TARGETSHELL).pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/fil.$(TARGETSHELL).pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

gen: $(INTERFACESSRC)/gen.$(TARGETSHELL).pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/gen.$(TARGETSHELL).pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

gen90: $(INTERFACESSRC)/gen90.$(TARGETSHELL).pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/gen90.$(TARGETSHELL).pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

hrb: $(INTERFACESSRC)/hrb.$(TARGETSHELL).pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/hrb.$(TARGETSHELL).pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

ipopt: $(INTERFACESSRC)/ipopt.$(TARGETSHELL).pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/ipopt.$(TARGETSHELL).pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

la04: $(INTERFACESSRC)/la04.$(TARGETSHELL).pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/la04.$(TARGETSHELL).pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

lmb: $(INTERFACESSRC)/lmb.$(TARGETSHELL).pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/lmb.$(TARGETSHELL).pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

derchk: $(INTERFACESSRC)/derchk.$(TARGETSHELL).pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/derchk.$(TARGETSHELL).pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

genc: $(INTERFACESSRC)/genc.$(TARGETSHELL).pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/genc.$(TARGETSHELL).pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

mns: $(INTERFACESSRC)/mns.$(TARGETSHELL).pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/mns.$(TARGETSHELL).pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

mx: $(INTERFACESSRC)/mx.$(TARGETSHELL).pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/mx.$(TARGETSHELL).pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

mx_pre: $(INTERFACESSRC)/mx_pre.$(TARGETSHELL).pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/mx_pre.$(TARGETSHELL).pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

nps: $(INTERFACESSRC)/nps.$(TARGETSHELL).pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/nps.$(TARGETSHELL).pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

osl: $(INTERFACESSRC)/osl.$(TARGETSHELL).pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/osl.$(TARGETSHELL).pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

pds: $(INTERFACESSRC)/pds.$(TARGETSHELL).pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/pds.$(TARGETSHELL).pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

prx: $(INTERFACESSRC)/prx.$(TARGETSHELL).pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/prx.$(TARGETSHELL).pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

nits: $(INTERFACESSRC)/nits.$(TARGETSHELL).pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/nits.$(TARGETSHELL).pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

snp: $(INTERFACESSRC)/snp.$(TARGETSHELL).pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/snp.$(TARGETSHELL).pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

stn: $(INTERFACESSRC)/stn.$(TARGETSHELL).pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/stn.$(TARGETSHELL).pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

stats: $(INTERFACESSRC)/stats.$(TARGETSHELL).pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/stats.$(TARGETSHELL).pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

ten: $(INTERFACESSRC)/ten.$(TARGETSHELL).pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/ten.$(TARGETSHELL).pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

unc: $(INTERFACESSRC)/unc.$(TARGETSHELL).pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/unc.$(TARGETSHELL).pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

cgp: $(INTERFACESSRC)/cgp.$(TARGETSHELL).pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/cgp.$(TARGETSHELL).pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

cgd: $(INTERFACESSRC)/cgd.$(TARGETSHELL).pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/cgd.$(TARGETSHELL).pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

lbs: $(INTERFACESSRC)/lbs.$(TARGETSHELL).pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/lbs.$(TARGETSHELL).pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

lbb: $(INTERFACESSRC)/lbb.$(TARGETSHELL).pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/lbb.$(TARGETSHELL).pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

trn: $(INTERFACESSRC)/trn.$(TARGETSHELL).pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/trn.$(TARGETSHELL).pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

va15: $(INTERFACESSRC)/va15.$(TARGETSHELL).pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/va15.$(TARGETSHELL).pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

ve09: $(INTERFACESSRC)/ve09.$(TARGETSHELL).pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/ve09.$(TARGETSHELL).pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

ve12: $(INTERFACESSRC)/ve12.$(TARGETSHELL).pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/ve12.$(TARGETSHELL).pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

ve14: $(INTERFACESSRC)/ve14.$(TARGETSHELL).pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/ve14.$(TARGETSHELL).pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

vf13: $(INTERFACESSRC)/vf13.$(TARGETSHELL).pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/vf13.$(TARGETSHELL).pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

envcheck: $(INTERFACESSRC)/envcheck.$(TARGETSHELL).pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/envcheck.$(TARGETSHELL).pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

helpmsg: $(INTERFACESSRC)/helpmsg.$(TARGETSHELL).pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/helpmsg.$(TARGETSHELL).pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

runpackage: $(INTERFACESSRC)/runpackage.$(TARGETSHELL).pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/runpackage.$(TARGETSHELL).pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

runcuter: $(INTERFACESSRC)/runcuter.$(TARGETSHELL).pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/runcuter.$(TARGETSHELL).pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

purge:: 
	$(RM) $(INTERFACES_AT_LARGE) $(C_INTERFACES) $(LOQO_INTERFACES) $(KNITRO_INTERFACES) $(TAO_INTERFACES)
all:: $(INTERFACES_AT_LARGE) $(C_INTERFACES) $(LOQO_INTERFACES) $(KNITRO_INTERFACES) $(TAO_INTERFACES)

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
