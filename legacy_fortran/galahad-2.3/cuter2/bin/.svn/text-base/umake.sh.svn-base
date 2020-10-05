#!/bin/bash
# umake : a very simple imake
# D. Orban @ ECE, 2002.

# Record number of input arguments
let n=$#

# Set default config_dir as 'config' in current directory
CONFIG_DIR='config'
UMAKEFILE=''
Doptions=()
Uoptions=()
Ioptions=()
Poptions=()
let error=0

# Process arguments
let i=1

while ((  i <= n  )); do
    opt=${!i}
    Dopt=`echo $opt | egrep '^-D'`
    Uopt=`echo $opt | egrep '^-U'`
    Iopt=`echo $opt | egrep '^-I'`
    Copt=`echo $opt | egrep '^-C' | sed -e 's/^-C//'`
    Popt=`echo $opt | egrep '^-P' | sed -e 's/^-P/-/'`
    if [[  "$Dopt" != ""  ]]; then          # Symbol is defined
        Doptions=( "${Doptions[@]} $Dopt" )
    elif [[  "$Uopt" != ""  ]]; then        # Symbol is undefined
        Uoptions=( "${Uoptions[@]} $Uopt" )
    elif [[  "$Iopt" != ""  ]]; then        # Include dir given
        Ioptions=( "${Ioptions[@]} $Iopt" )
    elif [[  "$Copt" != ""  ]]; then        # Config dir given
        CONFIG_DIR=$Copt
    elif [[  "$Popt" != ""  ]]; then        # CPP option given
        Poptions=( "${Poptions[@]} $Popt" )
    else                                    # Umakefile given
        UMAKEFILE=$opt
    fi
    (( i++ ))
done

# If Umakefile is not given, try file called Umakefile

if [[  "$UMAKEFILE" == ""  ]]; then
    if [[  -f ./Umakefile  ]]; then
        UMAKEFILE='Umakefile'
    else
        echo "No valid Umakefile found or given"
        exit 1
    fi
fi

[[  "$UMAKEFILE" == "" || "$CONFIG_DIR" == ""  ]] && let error=1

# Usage message

if ((  n < 1 || error == 1  )); then
    echo 'Usage:  umake umakefile [-Dsymbol...] [-Iinclude_dir...] [-Cconfig_dir]'
    echo 'where'
    echo '  umakefile is to be substituted'
    echo '            into the Umake configuration files and'
    echo '            subsequently processed by cpp,'
    echo '            the C PreProcessor. If more than one'
    echo '            umakefile is given, only the last is'
    echo '            processed and the others are discarded'
    echo '  -Dsymbol  defines the symbol <symbol> to be passed'
    echo '            to the C PreProcessor'
    echo '  -Usymbol  undefines the symbol <symbol> to be passed'
    echo '            to the C PreProcessor'
    echo '  -Iinclude_dir defines include dir.'
    echo '            Defaults to -I. -Iconfig'
    echo '  -Cconfig_dir specifies location of Umake.tmpl'
    echo '            Defaults to -Cconfig'
    echo '  -Poption  passes the option -option to cpp.'
fi

TMPL_FILE=${CONFIG_DIR}/Umake.tmpl

# Include default -I options if none are given
(( ${#Ioptions} == 0 )) && Ioptions='-I. -Iconfig'

# If sed files are not in $CUTER/build/scripts, please update
SEDLOC=${CUTER}/build/scripts

# Debug
#echo " Umakefile : $UMAKEFILE"
#echo " CONFIG_DIR: $CONFIG_DIR"
#echo " Doptions  : $Doptions"
#echo " Ioptions  : $Ioptions"
#echo " Poptions  : $Poptions"

# Read in basic commands

. ${MYCUTER}/cmds.basic

# Perform two passes of preprocessing through the Umakefile,
# substitute backslashes and newlines and remove
# consecutive blank lines (leaving only one).

UCPP=( "${CPP} -P ${Poptions[@]} ${Ioptions[@]} ${Doptions[@]} ${Uoptions[@]}" )

# If current dir already contains a Makefile, save it
if [[  -f Makefile  ]]; then
    ${MV} Makefile Makefile.bak
fi

# Protect tabs to account for CPPs which erase them
if [[  ! -f ${CONFIG_DIR}/Umake.rules.tabsafe  ]]; then
    ${SED} -f ${SEDLOC}/protect.sed ${CONFIG_DIR}/Umake.rules > ${CONFIG_DIR}/Umake.rules.tabsafe
fi
if [[  ! -f ${UMAKEFILE}.tabsafe  ]]; then
    ${SED} -f ${SEDLOC}/protect.sed ${UMAKEFILE} > ${UMAKEFILE}.tabsafe
fi
if [[  -f makefile.cmds  ]]; then
    ${SED} -f ${SEDLOC}/protect.sed makefile.cmds > makefile.cmds.tabsafe
fi

# Protect tabs again in included files, then recover them
${SED} -e "s/INCLUDE_UMAKEFILE/<${UMAKEFILE}.tabsafe>/" ${TMPL_FILE} | ${UCPP[@]} | ${SED} -f ${SEDLOC}/protect.sed | ${UCPP[@]} | ${SED} -e 's/BKSL/\\/g' -f ${SEDLOC}/newline.sed -e 's/^ *//g' | ${SED} -f ${SEDLOC}/flush.sed | ${SED} -f ${SEDLOC}/gettabs.sed | ${SED} -e 's/ *\\ *$/ \\/g' > Makefile
#${SED} -e "s/INCLUDE_UMAKEFILE/<${UMAKEFILE}.tabsafe>/" ${TMPL_FILE} > tmp1
#${UCPP[@]} tmp1 > tmp2
#${SED} -f ${SEDLOC}/protect.sed tmp2 > tmp3
#${UCPP[@]} tmp3 > tmp4
#${SED} -e 's/BKSL/\\/g' -f ${SEDLOC}/newline.sed -e 's/^ *//g' tmp4 > tmp5
#${SED} -f ${SEDLOC}/flush.sed tmp5 > tmp6
#${SED} -f ${SEDLOC}/gettabs.sed tmp6 > tmp7
#${SED} -e 's/ *\\ *$/ \\/g' tmp7 > Makefile

if [[  -f ${CONFIG_DIR}/Umake.rules.tabsafe  ]]; then
    ${RM} ${CONFIG_DIR}/Umake.rules.tabsafe
fi
${RM} Umakefile.tabsafe
${RM} makefile.cmds.tabsafe

exit 0
