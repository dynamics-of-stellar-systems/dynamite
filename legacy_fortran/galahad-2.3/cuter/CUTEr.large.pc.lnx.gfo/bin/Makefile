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

TOP = ..
CURRENT_DIR = /bin

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

.SUFFIXES: .pro

INTERFACESSRC = $(CUTER)/build/prototypes
INTERFACESBIN = $(MYCUTER)/bin
INTCONFIG = $(MYCUTER)/$(PRECISION)/config
DESTDIR =

INTERFACES_AT_LARGE = cob	 		 envcheck 		 fil	 		 gen	 		 gen90 		 helpmsg 		 hrb	 ipopt la04 		 lmb	 		 mns	 		 mx 		 nps	 		 osl	 		 pds	 		 prx	 		 nits	 		 runpackage 		 sdcob 		 sdfil 		 sdgen 		 sdgen90 		 sdhrb sdipopt sdla04 		 sdlmb 		 sdmns 		 sdmx 		 sdnits 		 sdnps 		 sdosl 		 sdpds 		 sdprx 		 sdsnp 		 sdstats 		 sdstn 		 sdten 		 sdunc 		 sdcgp 		 sdcgd 		 sdlbs 		 sdlbb 		 sdtrn 		 sdva15 		 sdve09 		 sdve12 		 sdve14 		 sdvf13 		 snp	 		 stats 		 stn	 		 ten	 		 unc	 		 cgp	 		 cgd	 		 lbs	 		 lbb	 		 trn	 		 va15 		 ve09 		 ve12 		 ve14 		 vf13

TAO_INTERFACES =

all::

baseconfig:
	cd $(MYCUTER)/$(PRECISION)/config; 	$(MAKE) $(MFLAGS) all;

all:: baseconfig

INTERFACES_DEPS = $(INTCONFIG)/cmds $(INTCONFIG)/this_build $(INTCONFIG)/script.sed

cob: $(INTERFACESSRC)/cob.pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/cob.pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

fil: $(INTERFACESSRC)/fil.pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/fil.pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

gen: $(INTERFACESSRC)/gen.pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/gen.pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

gen90: $(INTERFACESSRC)/gen90.pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/gen90.pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

hrb: $(INTERFACESSRC)/hrb.pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/hrb.pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

ipopt: $(INTERFACESSRC)/ipopt.pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/ipopt.pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

la04: $(INTERFACESSRC)/la04.pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/la04.pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

lmb: $(INTERFACESSRC)/lmb.pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/lmb.pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

mns: $(INTERFACESSRC)/mns.pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/mns.pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

mx: $(INTERFACESSRC)/mx.pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/mx.pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

nps: $(INTERFACESSRC)/nps.pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/nps.pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

osl: $(INTERFACESSRC)/osl.pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/osl.pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

pds: $(INTERFACESSRC)/pds.pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/pds.pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

prx: $(INTERFACESSRC)/prx.pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/prx.pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

nits: $(INTERFACESSRC)/nits.pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/nits.pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

snp: $(INTERFACESSRC)/snp.pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/snp.pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

stn: $(INTERFACESSRC)/stn.pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/stn.pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

stats: $(INTERFACESSRC)/stats.pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/stats.pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

ten: $(INTERFACESSRC)/ten.pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/ten.pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

unc: $(INTERFACESSRC)/unc.pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/unc.pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

cgp: $(INTERFACESSRC)/cgp.pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/cgp.pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

cgd: $(INTERFACESSRC)/cgd.pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/cgd.pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

lbs: $(INTERFACESSRC)/lbs.pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/lbs.pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

lbb: $(INTERFACESSRC)/lbb.pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/lbb.pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

trn: $(INTERFACESSRC)/trn.pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/trn.pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

va15: $(INTERFACESSRC)/va15.pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/va15.pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

ve09: $(INTERFACESSRC)/ve09.pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/ve09.pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

ve12: $(INTERFACESSRC)/ve12.pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/ve12.pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

ve14: $(INTERFACESSRC)/ve14.pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/ve14.pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

vf13: $(INTERFACESSRC)/vf13.pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/vf13.pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

envcheck: $(INTERFACESSRC)/envcheck.pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/envcheck.pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

helpmsg: $(INTERFACESSRC)/helpmsg.pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/helpmsg.pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

runpackage: $(INTERFACESSRC)/runpackage.pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/runpackage.pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

sdcob: $(INTERFACESSRC)/sdcob.pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/sdcob.pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

sdfil: $(INTERFACESSRC)/sdfil.pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/sdfil.pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

sdgen: $(INTERFACESSRC)/sdgen.pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/sdgen.pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

sdgen90: $(INTERFACESSRC)/sdgen90.pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/sdgen90.pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

sdhrb: $(INTERFACESSRC)/sdhrb.pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/sdhrb.pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

sdipopt: $(INTERFACESSRC)/sdipopt.pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/sdipopt.pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

sdla04: $(INTERFACESSRC)/sdla04.pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/sdla04.pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

sdlmb: $(INTERFACESSRC)/sdlmb.pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/sdlmb.pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

sdmns: $(INTERFACESSRC)/sdmns.pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/sdmns.pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

sdmx: $(INTERFACESSRC)/sdmx.pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/sdmx.pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

sdnps: $(INTERFACESSRC)/sdnps.pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/sdnps.pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

sdosl: $(INTERFACESSRC)/sdosl.pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/sdosl.pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

sdpds: $(INTERFACESSRC)/sdpds.pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/sdpds.pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

sdprx: $(INTERFACESSRC)/sdprx.pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/sdprx.pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

sdnits: $(INTERFACESSRC)/sdnits.pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/sdnits.pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

sdsnp: $(INTERFACESSRC)/sdsnp.pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/sdsnp.pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

sdstats: $(INTERFACESSRC)/sdstats.pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/sdstats.pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

sdstn: $(INTERFACESSRC)/sdstn.pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/sdstn.pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

sdten: $(INTERFACESSRC)/sdten.pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/sdten.pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

sdunc: $(INTERFACESSRC)/sdunc.pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/sdunc.pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

sdcgp: $(INTERFACESSRC)/sdcgp.pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/sdcgp.pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

sdcgd: $(INTERFACESSRC)/sdcgd.pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/sdcgd.pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

sdlbs: $(INTERFACESSRC)/sdlbs.pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/sdlbs.pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

sdlbb: $(INTERFACESSRC)/sdlbb.pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/sdlbb.pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

sdtrn: $(INTERFACESSRC)/sdtrn.pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/sdtrn.pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

sdva15: $(INTERFACESSRC)/sdva15.pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/sdva15.pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

sdve09: $(INTERFACESSRC)/sdve09.pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/sdve09.pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

sdve12: $(INTERFACESSRC)/sdve12.pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/sdve12.pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

sdve14: $(INTERFACESSRC)/sdve14.pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/sdve14.pro > $@ 
	$(CHMOD) 755 $@ 
	@printf '[Ok]\n'

sdvf13: $(INTERFACESSRC)/sdvf13.pro $(INTERFACES_DEPS) 
	$(RM) $@ 
	printf ' %-12s %-15s\t\t' "Casting" $@ 
	$(SED) -f $(INTCONFIG)/script.sed $(INTERFACESSRC)/sdvf13.pro > $@ 
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
	$(UMAKE_CMD) -DTOPDIR=$(TOP) -DCURDIR=$(CURRENT_DIR) -Dlinux -DIsgfo -DLargeSize -DDoublePrecision -DIsa90Compiler -DNoCcomp -Ptraditional -Pw

Makefiles::

includes::

