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

TOP = ..
CURRENT_DIR = /bin

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

.SUFFIXES: .pro

INTERFACESSRC = $(SIFDEC)/build/prototypes
INTERFACESBIN = $(MYSIFDEC)/bin
INTCONFIG = $(MYSIFDEC)/$(PRECISION)/config
DESTDIR =

INTERFACES_AT_LARGE_DOT_PRO = 		 $(INTERFACESSRC)/classall.pro 		 $(INTERFACESSRC)/classify.pro 		 $(INTERFACESSRC)/helpmsg.pro 		 $(INTERFACESSRC)/select.pro 		 $(INTERFACESSRC)/sifdecode.pro

INTERFACES_AT_LARGE = classall 		 classify 		 helpmsg 		 select 		 sifdecode

all::

INTERFACES_DEPS = $(INTCONFIG)/cmds $(INTCONFIG)/this_build $(INTCONFIG)/script.sed

classall: $(INTERFACESSRC)/classall.pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/classall.pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

classify: $(INTERFACESSRC)/classify.pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/classify.pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

helpmsg: $(INTERFACESSRC)/helpmsg.pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/helpmsg.pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

select: $(INTERFACESSRC)/select.pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/select.pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

sifdecode: $(INTERFACESSRC)/sifdecode.pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/sifdecode.pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

purge:: 
	$(RM) $(INTERFACES_AT_LARGE)

show.awk: $(INTERFACESSRC)/show.awk 
	@printf ' %-12s %-15s\t\t' "Installing" show.awk 
	$(INSTALL) $(INSTALLFLAGS) $(INTERFACESSRC)/show.awk $(INTERFACESBIN)/show.awk 
	@printf '[Ok]\n'

purge:: 
	$(RM) show.awk

param.awk: $(INTERFACESSRC)/param.awk 
	@printf ' %-12s %-15s\t\t' "Installing" param.awk 
	$(INSTALL) $(INSTALLFLAGS) $(INTERFACESSRC)/param.awk $(INTERFACESBIN)/param.awk 
	@printf '[Ok]\n'

purge:: 
	$(RM) param.awk

SLCTSRCDIR = $(SIFDEC)/common/src/select
SLCTSRCS = $(SLCTSRCDIR)/slct.f $(SLCTSRCDIR)/clsf.f

slct: $(SLCTSRCDIR)/slct.f 
	$(RM) slct 
	@printf ' %-12s %-15s\t\t' "Linking" slct 
	$(COMPILE) $(FFLAGS) -o $@.o $(SLCTSRCDIR)/$@.f 
	$(LOAD) $(FFLAGS) -o $@ $(F_LDOPTIONS) $@.o $(LDLIBS) $(EXTRA_LOAD_FLAGS) 
	@printf '[Ok]\n' 
all:: slct 

purge:: 
	$(RM) slct

clsf: $(SLCTSRCDIR)/clsf.f 
	$(RM) clsf 
	@printf ' %-12s %-15s\t\t' "Linking" clsf 
	$(COMPILE) $(FFLAGS) -o $@.o $(SLCTSRCDIR)/$@.f 
	$(LOAD) $(FFLAGS) -o $@ $(F_LDOPTIONS) $@.o $(LDLIBS) $(EXTRA_LOAD_FLAGS) 
	@printf '[Ok]\n' 
all:: clsf 

purge:: 
	$(RM) clsf

clean:: 
	$(RM) slct.o clsf.o
all:: slct clsf

all:: $(INTERFACES_AT_LARGE) show.awk param.awk

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

