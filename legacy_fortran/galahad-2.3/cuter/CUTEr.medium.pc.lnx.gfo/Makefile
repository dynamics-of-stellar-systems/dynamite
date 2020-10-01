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

TOP = .
CURRENT_DIR = .

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

.SUFFIXES: .pro .f .f90

.f.o:
	true

.f90.o:
	true

INTERFACESSRC = $(CUTER)/build/prototypes
INTCONFIG = $(MYCUTER)/$(PRECISION)/config

all::

baseconfig:
	cd $(MYCUTER)/$(PRECISION)/config; 	$(MAKE) $(MFLAGS) all;

all:: baseconfig

clean:: 
	$(RM) *.CKP *.ln *.BAK *.bak core errs ,* *~ tags TAGS make.log MakeOut "#"*

SUBDIRS = bin double

all:: 
	@case '${MFLAGS}' in *[ik]*) set +e;; esac; \
	for i in $(SUBDIRS) ; \
	do \
		(cd $$i ; echo "making" all "in $(CURRENT_DIR)/$$i..."; \
			$(MAKE) $(MFLAGS) all); \
	done

purge:: 
	$(RM) $(INTCONFIG)/cmds

purge:: 
	$(RM) $(INTCONFIG)/local.f

MFLAGS =

emptyrule::

clean:: 
	$(RM) *.CKP *.ln *.BAK *.bak core errs ,* *~ tags TAGS make.log MakeOut "#"*

Makefile:: 
	-@if [ -f Makefile ]; then set -x; \
	$(RM) Makefile.bak; $(MV) Makefile Makefile.bak; \
	else exit 0; fi 
	$(UMAKE_CMD) -DTOPDIR=$(TOP) -DCURDIR=$(CURRENT_DIR) -Dlinux -DIsgfo -DMediumSize -DDoublePrecision -DIsa90Compiler -DIsgcc -UNoCcomp -Ptraditional -Pw

clean:: 
	@case '${MFLAGS}' in *[ik]*) set +e;; esac; \
	for i in $(SUBDIRS) ; \
	do \
		(cd $$i ; echo "cleaning" "in $(CURRENT_DIR)/$$i..."; \
			$(MAKE) $(MFLAGS) clean); \
	done

purge:: 
	@case '${MFLAGS}' in *[ik]*) set +e;; esac; \
	for i in $(SUBDIRS) ; \
	do \
		(cd $$i ; echo "purging" "in $(CURRENT_DIR)/$$i..."; \
			$(MAKE) $(MFLAGS) purge); \
	done

Makefiles:: 
	-@case '${MFLAGS}' in *[ik]*) set +e;; esac; \
	case '${MFLAGS}' in *n*) executeit="no";; esac; \
	for i in $(SUBDIRS) ; \
	do \
		case "$(CURRENT_DIR)" in \
		.) curdir= ;; \
		*) curdir=$(CURRENT_DIR) ;; \
		esac; \
		echo "making Makefiles in $$curdir/$$i..."; \
		itmp="$$i" \
		curtmp="$(CURRENT_DIR)" \
		toptmp=""; \
		case "$$itmp" in \
		 ../?*) \
			while echo "$$itmp" | grep '^../' > /dev/null; \
			 do \
			 toptmp="/`basename $$curtmp`$$toptmp"; \
			 curtmp="`dirname $$curtmp`"; \
			 itmp="`echo $$itmp | sed 's;../;;'`"; \
			 done \
		 ;; \
		esac; \
		case "$$itmp" in \
		*/?*/?*/?*/?*)	newtop=../../../../..;; \
		*/?*/?*/?*)	newtop=../../../..;; \
		*/?*/?*)	newtop=../../..;; \
		*/?*)		newtop=../..;; \
		*)		newtop=..;; \
		esac; \
		newtop="$$newtop$$toptmp"; \
		case "$(TOP)" in \
		/?*) umaketop=$(TOP) \
		 umakeprefix= ;; \
		.) umaketop=$$newtop \
		 umakeprefix=$$newtop/ ;; \
		*) umaketop=$$newtop/$(TOP) \
		 umakeprefix=$$newtop/ ;; \
		esac; \
		cd $$i; \
		$(RM) Makefile.bak; \
		if [ -f Makefile ]; then \
			echo "	$(MV) Makefile Makefile.bak"; \
			if [ "$$executeit" != "no" ]; then \
				$(MV) Makefile Makefile.bak; \
			fi; \
		fi; \
		if [ "$$executeit" != "no" ]; then \
			$(UMAKE_CMD) -I$$umakeprefix$(IRULESRC) -Dlinux -DIsgfo -DMediumSize -DDoublePrecision -DIsa90Compiler -DIsgcc -UNoCcomp -Ptraditional -Pw -DTOPDIR=$$umaketop -DCURDIR=$$curdir/$$i; \
		fi; \
		$(MAKE) $(MFLAGS) Makefiles; \
		cd $$newtop; \
	done

includes:: 
	@case '${MFLAGS}' in *[ik]*) set +e;; esac; \
	for i in $(SUBDIRS) ; \
	do \
		(cd $$i ; echo including "in $(CURRENT_DIR)/$$i..."; \
			$(MAKE) $(MFLAGS) includes); \
	done

