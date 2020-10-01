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
CCOMPILE	 = gcc -c
CLOAD		 = gcc
CFLAGS		 = -O
SPECIALLIBS = 
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

SIZE = medium
SIZ = CMED
SIZ90 = \!MED

PATHSEP = /
SHELL = /bin/sh

TOP = ../..
CURRENT_DIR = /double/config

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

INTCONFIG = $(MYCUTER)/$(PRECISION)/config

INTERFACESSRC = $(CUTER)/build/prototypes

CMDS = cmds

all::

slash.sed:
	$(LN) -s $(CUTER)/build/scripts/slash.sed ./

${CMDS} : $(INTERFACESSRC)/cmds.pro slash.sed
	$(RM) cmds
	@printf ' %-27s\t\t' "Getting UNIX commands right"
	echo 's/\$${SYSTEM}/'`echo $(SYSTEM) | $(SED) -f slash.sed`'/g' > castscript.sed
	echo 's/\$${KEYSYS}/'`echo $(KEYSYS) | $(SED) -f slash.sed`'/g' >> castscript.sed
	echo 's/\$${TMP}/'`echo $(TMP) | $(SED) -f slash.sed`'/g' >> castscript.sed
	echo 's/\$${MKDIR}/'`echo $(MKDIR) | $(SED) -f slash.sed`'/g' >> castscript.sed
	echo 's/\$${CP}/'`echo $(CP) | $(SED) -f slash.sed`'/g' >> castscript.sed
	echo 's/\$${RM}/'`echo $(RM) | $(SED) -f slash.sed`'/g' >> castscript.sed
	echo 's/\$${MV}/'`echo $(MV) | $(SED) -f slash.sed`'/g' >> castscript.sed
	echo 's/\$${CAT}/'`echo $(CAT) | $(SED) -f slash.sed`'/g' >> castscript.sed
	echo 's/\$${CHMOD}/'`echo $(CHMOD) | $(SED) -f slash.sed`'/g' >> castscript.sed
	echo 's/\$${SED}/'`echo $(SED) | $(SED) -f slash.sed`'/g' >> castscript.sed
	echo 's/\$${LN}/'`echo $(LN) | $(SED) -f slash.sed`'/g' >> castscript.sed
	echo 's/\$${LS}/'`echo $(LS) | $(SED) -f slash.sed`'/g' >> castscript.sed
	echo 's/\$${AR}/'`echo $(AR) | $(SED) -f slash.sed`'/g' >> castscript.sed
	echo 's/\$${RMDIR}/'`echo $(RMDIR) | $(SED) -f slash.sed`'/g' >> castscript.sed
	echo 's/\$${GREP}/'`echo $(GREP) | $(SED) -f slash.sed`'/g' >> castscript.sed
	echo 's/\$${AWK}/'`echo $(AWK) | $(SED) -f slash.sed`'/g' >> castscript.sed
	echo 's/\$${HEAD}/'`echo $(HEAD) | $(SED) -f slash.sed`'/g' >> castscript.sed
	echo 's/\$${TAIL}/'`echo $(TAIL) | $(SED) -f slash.sed`'/g' >> castscript.sed
	echo 's/\$${WC}/'`echo $(WC) | $(SED) -f slash.sed`'/g' >> castscript.sed
	echo 's/\$${MAKE}/'`echo $(MAKE) | $(SED) -f slash.sed`'/g' >> castscript.sed
	echo 's/\$${COMPILE}/'`echo $(COMPILE) | $(SED) -f slash.sed`'/g' >> castscript.sed
	echo 's/\$${LOAD}/'`echo $(LOAD) | $(SED) -f slash.sed`'/g' >> castscript.sed
	echo 's/\$${ISF9095}/'`echo $(ISF9095) | $(SED) -f slash.sed`'/g' >> castscript.sed
	echo 's/\$${COMPILE9095}/'`echo $(COMPILE9095) | $(SED) -f slash.sed`'/g' >> castscript.sed
	echo 's/\$${LOAD9095}/'`echo $(LOAD9095) | $(SED) -f slash.sed`'/g' >> castscript.sed
	echo 's/\$${FFLAGS}/'`echo $(FFLAGS) | $(SED) -f slash.sed`'/g' >> castscript.sed
	echo 's/\$${SPECIALFLAGS}/'`echo $(SPECIALFLAGS) | $(SED) -f slash.sed`'/g' >> castscript.sed
	echo 's/\$${CCOMPILE}/'`echo $(CCOMPILE) | $(SED) -f slash.sed`'/g' >> castscript.sed
	echo 's/\$${CLOAD}/'`echo $(CLOAD) | $(SED) -f slash.sed`'/g' >> castscript.sed
	echo 's/\$${CFLAGS}/'`echo $(CFLAGS) | $(SED) -f slash.sed`'/g' >> castscript.sed
	echo 's/\$${SPECIALLIBS}/'`echo $(SPECIALLIBS) | $(SED) -f slash.sed`'/g' >> castscript.sed
	echo 's/\$${MEXFORTRAN}/'`echo $(MEXFORTRAN) | $(SED) -f slash.sed`'/g' >> castscript.sed
	echo 's/\$${MEXFFLAGS}/'`echo $(MEXFFLAGS) | $(SED) -f slash.sed`'/g' >> castscript.sed
	$(SED) -f castscript.sed $(INTERFACESSRC)/cmds.pro > $@
	$(RM) castscript.sed
	@printf '[Ok]\n'

all:: slash.sed ${CMDS}

clean:: slash.sed

purge:: ${CMDS}

script.sed: $(INTERFACESSRC)/script.sed.pro
	echo 's/\$${PRECISION}/'`echo $(PRECISION) | $(SED) -f slash.sed`'/g' > castscript.sed
	echo 's/\$${MYCUTER}/'`echo $(MYCUTER) | $(SED) -f slash.sed`'/g' >> castscript.sed
	echo 's/\$${CUTER}/'`echo $(CUTER) | $(SED) -f slash.sed`'/g' >> castscript.sed
	$(SED) -f castscript.sed $? > $@
	$(RM) castscript.sed

cast.sed: $(INTERFACESSRC)/cast.sed.pro
	echo 's/\$${SIZ}/'`echo $(SIZ) | $(SED) -f slash.sed`'/g' > castscript.sed
	echo 's/\$${KEYSYS}/'`echo $(KEYSYS) | $(SED) -f slash.sed`'/g' >> castscript.sed
	echo 's/\$${PRC}/'`echo $(PRC) | $(SED) -f slash.sed`'/g' >> castscript.sed
	$(SED) -f castscript.sed $? > $@
	$(RM) castscript.sed

cast90.sed: $(INTERFACESSRC)/cast90.sed.pro
	echo 's/\$${SIZ90}/'`echo $(SIZ90) | $(SED) -f slash.sed`'/g' > castscript.sed
	echo 's/\$${PRC90}/'`echo $(PRC90) | $(SED) -f slash.sed`'/g' >> castscript.sed
	$(SED) -f castscript.sed $? > $@
	$(RM) castscript.sed

ipoptma.inc: $(INTERFACESSRC)/ipoptma.inc.pro 
	$(RM) ipoptma.inc 
	@printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/cast.sed $(INTERFACESSRC)/ipoptma.inc.pro > ipoptma.inc 
	@printf '[Ok]\n'

ipoptma.wrk: $(INTERFACESSRC)/ipoptma.wrk.pro 
	$(RM) ipoptma.wrk 
	@printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/cast.sed $(INTERFACESSRC)/ipoptma.wrk.pro > ipoptma.wrk 
	@printf '[Ok]\n'

LOCAL = $(CUTER)/common/src/local.f

LOCAL_OBJ = local.o

${LOCAL_OBJ}: $(LOCAL) $(INTCONFIG)/timer $(INTCONFIG)/cast.sed
	$(RM) $@
	$(RM) local.ff
	@printf ' %-12s %-15s\t\t' "Casting" local.f
	echo 's/\$${NBYTES}/'`echo $(NBYTES) | $(SED) -f slash.sed`'/g' > castscript.sed
	echo 's/\$${R1}/'`echo $(R1) | $(SED) -f slash.sed`'/g' >> castscript.sed
	echo 's/\$${R2}/'`echo $(R2) | $(SED) -f slash.sed`'/g' >> castscript.sed
	echo 's/\$${R3}/'`echo $(R3) | $(SED) -f slash.sed`'/g' >> castscript.sed
	echo 's/\$${R4}/'`echo $(R4) | $(SED) -f slash.sed`'/g' >> castscript.sed
	echo 's/\$${R5}/'`echo $(R5) | $(SED) -f slash.sed`'/g' >> castscript.sed
	echo 's/\$${D1}/'`echo $(D1) | $(SED) -f slash.sed`'/g' >> castscript.sed
	echo 's/\$${D2}/'`echo $(D2) | $(SED) -f slash.sed`'/g' >> castscript.sed
	echo 's/\$${D3}/'`echo $(D3) | $(SED) -f slash.sed`'/g' >> castscript.sed
	echo 's/\$${D4}/'`echo $(D4) | $(SED) -f slash.sed`'/g' >> castscript.sed
	echo 's/\$${D5}/'`echo $(D5) | $(SED) -f slash.sed`'/g' >> castscript.sed
	$(SED) -f castscript.sed $(LOCAL) > local.ff
	$(RM) castscript.sed
	$(CHMOD) 666 local.ff
	@printf '[Ok]\n'
	@printf ' %-12s %-15s\t\t' "Adding" timer
	$(CAT) $(INTCONFIG)/timer >> $(INTCONFIG)/local.ff
	@printf '[Ok]\n'
	@printf ' %-12s %-15s\t\t' "Building" $@
	$(SED) -f $(INTCONFIG)/cast.sed $(INTCONFIG)/${@:.o=.ff} > ${@:.o=.f}
	$(COMPILE) $(FFLAGS) -o $@ ${@:.o=.f}
	@printf '[Ok]\n'

all:: script.sed

all:: cast.sed cast90.sed ipoptma.inc ipoptma.wrk

all:: $(LOCAL_OBJ)

clean:: 
	$(RM) *.CKP *.ln *.BAK *.bak core errs ,* *~ tags TAGS make.log MakeOut "#"*

purge:: 
	$(RM) $(CMDS)

purge:: 
	$(RM) $(LOCAL_OBJ)

MFLAGS =

emptyrule::

clean:: 
	$(RM) *.CKP *.ln *.BAK *.bak core errs ,* *~ tags TAGS make.log MakeOut "#"*

Makefile:: 
	-@if [ -f Makefile ]; then set -x; \
	$(RM) Makefile.bak; $(MV) Makefile Makefile.bak; \
	else exit 0; fi 
	$(UMAKE_CMD) -DTOPDIR=$(TOP) -DCURDIR=$(CURRENT_DIR) -Dlinux -DIsgfo -DMediumSize -DDoublePrecision -DIsa90Compiler -DIsgcc -UNoCcomp -Ptraditional -Pw

Makefiles::

includes::

