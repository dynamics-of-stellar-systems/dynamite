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
CURRENT_DIR = /double/config
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
INTCONFIG = $(MYCUTER)/$(PRECISION)/config
EXT =
INTERFACESSRC = $(CUTER)/build/prototypes
CMDSDOTPRO = $(INTERFACESSRC)/cmds.sh.pro
CMDS = cmds$(EXT)
CS = castscript$(EXT).sed
all::
${CMDS}: $(CMDSDOTPRO)
echo 's/\$${SYSTEM}/'`echo $(SYSTEM) | $(SED) -f $(CUTER)/build/scripts/slash.sed`'/g' > $(CS)
echo 's/\$${KEYSYS}/'`echo $(KEYSYS) | $(SED) -f $(CUTER)/build/scripts/slash.sed`'/g' >> $(CS)
echo 's/\$${TMP}/'`echo $(TMP) | $(SED) -f $(CUTER)/build/scripts/slash.sed`'/g' >> $(CS)
echo 's/\$${MKDIR}/'`echo $(MKDIR) | $(SED) -f $(CUTER)/build/scripts/slash.sed`'/g' >> $(CS)
echo 's/\$${CP}/'`echo $(CP) | $(SED) -f $(CUTER)/build/scripts/slash.sed`'/g' >> $(CS)
echo 's/\$${RM}/'`echo $(RM) | $(SED) -f $(CUTER)/build/scripts/slash.sed`'/g' >> $(CS)
echo 's/\$${MV}/'`echo $(MV) | $(SED) -f $(CUTER)/build/scripts/slash.sed`'/g' >> $(CS)
echo 's/\$${CAT}/'`echo $(CAT) | $(SED) -f $(CUTER)/build/scripts/slash.sed`'/g' >> $(CS)
echo 's/\$${CHMOD}/'`echo $(CHMOD) | $(SED) -f $(CUTER)/build/scripts/slash.sed`'/g' >> $(CS)
echo 's/\$${SED}/'`echo $(SED) | $(SED) -f $(CUTER)/build/scripts/slash.sed`'/g' >> $(CS)
echo 's/\$${LN}/'`echo $(LN) | $(SED) -f $(CUTER)/build/scripts/slash.sed`'/g' >> $(CS)
echo 's/\$${LS}/'`echo $(LS) | $(SED) -f $(CUTER)/build/scripts/slash.sed`'/g' >> $(CS)
echo 's/\$${AR}/'`echo $(AR) | $(SED) -f $(CUTER)/build/scripts/slash.sed`'/g' >> $(CS)
echo 's/\$${RMDIR}/'`echo $(RMDIR) | $(SED) -f $(CUTER)/build/scripts/slash.sed`'/g' >> $(CS)
echo 's/\$${GREP}/'`echo $(GREP) | $(SED) -f $(CUTER)/build/scripts/slash.sed`'/g' >> $(CS)
echo 's/\$${AWK}/'`echo $(AWK) | $(SED) -f $(CUTER)/build/scripts/slash.sed`'/g' >> $(CS)
echo 's/\$${HEAD}/'`echo $(HEAD) | $(SED) -f $(CUTER)/build/scripts/slash.sed`'/g' >> $(CS)
echo 's/\$${TAIL}/'`echo $(TAIL) | $(SED) -f $(CUTER)/build/scripts/slash.sed`'/g' >> $(CS)
echo 's/\$${WC}/'`echo $(WC) | $(SED) -f $(CUTER)/build/scripts/slash.sed`'/g' >> $(CS)
echo 's/\$${MAKE}/'`echo $(MAKE) | $(SED) -f $(CUTER)/build/scripts/slash.sed`'/g' >> $(CS)
echo 's/\$${COMPILE}/'`echo $(COMPILE) | $(SED) -f $(CUTER)/build/scripts/slash.sed`'/g' >> $(CS)
echo 's/\$${LOAD}/'`echo $(LOAD) | $(SED) -f $(CUTER)/build/scripts/slash.sed`'/g' >> $(CS)
echo 's/\$${ISF9095}/'`echo $(ISF9095) | $(SED) -f $(CUTER)/build/scripts/slash.sed`'/g' >> $(CS)
echo 's/\$${COMPILE9095}/'`echo $(COMPILE9095) | $(SED) -f $(CUTER)/build/scripts/slash.sed`'/g' >> $(CS)
echo 's/\$${LOAD9095}/'`echo $(LOAD9095) | $(SED) -f $(CUTER)/build/scripts/slash.sed`'/g' >> $(CS)
echo 's/\$${FFLAGS}/'`echo $(FFLAGS) | $(SED) -f $(CUTER)/build/scripts/slash.sed`'/g' >> $(CS)
echo 's/\$${SPECIALFLAGS}/'`echo $(SPECIALFLAGS) | $(SED) -f $(CUTER)/build/scripts/slash.sed`'/g' >> $(CS)
echo 's/\$${FDEBUGFLAGS}/'`echo $(FDEBUGFLAGS) | $(SED) -f $(CUTER)/build/scripts/slash.sed`'/g' >> $(CS)
echo 's/\$${FPROFILEFLAGS}/'`echo $(FPROFILEFLAGS) | $(SED) -f $(CUTER)/build/scripts/slash.sed`'/g' >> $(CS)
echo 's/\$${CPP}/'`echo $(CPP) | $(SED) -f $(CUTER)/build/scripts/slash.sed`'/g' >> $(CS)
echo 's/\$${CCOMPILE}/'`echo $(CCOMPILE) | $(SED) -f $(CUTER)/build/scripts/slash.sed`'/g' >> $(CS)
echo 's/\$${CLOAD}/'`echo $(CLOAD) | $(SED) -f $(CUTER)/build/scripts/slash.sed`'/g' >> $(CS)
echo 's/\$${CFLAGS}/'`echo $(CFLAGS) | $(SED) -f $(CUTER)/build/scripts/slash.sed`'/g' >> $(CS)
echo 's/\$${SPECIALLIBS}/'`echo $(SPECIALLIBS) | $(SED) -f $(CUTER)/build/scripts/slash.sed`'/g' >> $(CS)
echo 's/\$${MEXFORTRAN}/'`echo $(MEXFORTRAN) | $(SED) -f $(CUTER)/build/scripts/slash.sed`'/g' >> $(CS)
echo 's/\$${MEXFFLAGS}/'`echo $(MEXFFLAGS) | $(SED) -f $(CUTER)/build/scripts/slash.sed`'/g' >> $(CS)
$(SED) -f $(CS) $(CMDSDOTPRO) > $@
$(RM) $(CS)
all:: ${CMDS}
clean::
purge:: ${CMDS}
script.sed: $(INTERFACESSRC)/script.sed.pro
	echo 's/\$${PRECISION}/'`echo $(PRECISION) | $(SED) -f $(CUTER)/build/scripts/slash.sed`'/g' > castscript.sed
	echo 's/\$${MYCUTER}/'`echo $(MYCUTER) | $(SED) -f $(CUTER)/build/scripts/slash.sed`'/g' >> castscript.sed
	echo 's/\$${CUTER}/'`echo $(CUTER) | $(SED) -f $(CUTER)/build/scripts/slash.sed`'/g' >> castscript.sed
	$(SED) -f castscript.sed $? > $@
	$(RM) castscript.sed
cast.sed: $(INTERFACESSRC)/cast.sed.pro
	echo 's/\$${SIZ}/'`echo $(SIZ) | $(SED) -f $(CUTER)/build/scripts/slash.sed`'/g' > castscript.sed
	echo 's/\$${KEYSYS}/'`echo $(KEYSYS) | $(SED) -f $(CUTER)/build/scripts/slash.sed`'/g' >> castscript.sed
	echo 's/\$${PRC}/'`echo $(PRC) | $(SED) -f $(CUTER)/build/scripts/slash.sed`'/g' >> castscript.sed
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
	echo 's/\$${NBYTES}/'`echo $(NBYTES) | $(SED) -f $(CUTER)/build/scripts/slash.sed`'/g' > castscript.sed
	echo 's/\$${R1}/'`echo $(R1) | $(SED) -f $(CUTER)/build/scripts/slash.sed`'/g' >> castscript.sed
	echo 's/\$${R2}/'`echo $(R2) | $(SED) -f $(CUTER)/build/scripts/slash.sed`'/g' >> castscript.sed
	echo 's/\$${R3}/'`echo $(R3) | $(SED) -f $(CUTER)/build/scripts/slash.sed`'/g' >> castscript.sed
	echo 's/\$${R4}/'`echo $(R4) | $(SED) -f $(CUTER)/build/scripts/slash.sed`'/g' >> castscript.sed
	echo 's/\$${R5}/'`echo $(R5) | $(SED) -f $(CUTER)/build/scripts/slash.sed`'/g' >> castscript.sed
	echo 's/\$${D1}/'`echo $(D1) | $(SED) -f $(CUTER)/build/scripts/slash.sed`'/g' >> castscript.sed
	echo 's/\$${D2}/'`echo $(D2) | $(SED) -f $(CUTER)/build/scripts/slash.sed`'/g' >> castscript.sed
	echo 's/\$${D3}/'`echo $(D3) | $(SED) -f $(CUTER)/build/scripts/slash.sed`'/g' >> castscript.sed
	echo 's/\$${D4}/'`echo $(D4) | $(SED) -f $(CUTER)/build/scripts/slash.sed`'/g' >> castscript.sed
	echo 's/\$${D5}/'`echo $(D5) | $(SED) -f $(CUTER)/build/scripts/slash.sed`'/g' >> castscript.sed
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
all:: cast.sed ipoptma.inc ipoptma.wrk
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
	$(UMAKE_CMD) -DTOPDIR=$(TOP) -DCURDIR=$(CURRENT_DIR) -Dlinux -DIsg77 -DMediumSize -DDoublePrecision -DIsgcc -UNoCcomp
Makefiles::
includes::
