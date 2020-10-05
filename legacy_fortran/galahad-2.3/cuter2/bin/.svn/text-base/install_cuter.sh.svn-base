#!/bin/bash
#
# Build script for CUTEr
#
# syntax: install_cuter 
#
# N.Gould, D. Orban & Ph. Toint
#

#  check input arguments (if any)

if [[  $# != 0  ]]; then
   echo "Use: install_cuter"
   exit 1
fi

export CUTER=$PWD
ARCH=$CUTER/build/arch

#
# Prelude to installation phase
#

UMAKE_COMMAND="${CUTER}/bin/umake.sh -I. -Iconfig"
UMAKE_OPTIONS=( " " )
UMAKE_CPP_OPTIONS=" "

CORRECT_PLATFORM="false"
while [[  $CORRECT_PLATFORM == "false"  ]]; do 
   echo ' Select platform'
   echo
   echo '   (1) Compaq (DEC) alpha'
   echo '   (2) Cray'
   echo '   (3) HP workstation'
   echo '   (4) IBM RS/6000'
   echo '   (5) PC'
   echo '   (6) SGI workstation'
   echo '   (7) SUN workstation'
   echo '   (8) Mac OSX'
   echo '   (0) Other (user-defined)'
   
   read CHOICE
   
   case  $CHOICE  in
       "1")
            CORRECT_PLATFORM="true"
            export MCH="alp"
            export MACHINE="Compaq-Alpha"
            UMAKE_OPTIONS=( "${UMAKE_OPTIONS[@]} -D__osf__" )

            CORRECT_OS="false"
            while [[  $CORRECT_OS == "false"  ]]; do 

               echo ' Select operating system'
               echo
               echo '   (1) Digital Unix'
               echo '   (2) Tru-64 Unix'
               echo '   (3) Linux'
               
               read CHOICE

               case  $CHOICE  in
               1)
                  CORRECT_OS="true"
                  export OS="dux"
                  ;;
               2)
                  CORRECT_OS="true"
                  export OS="t64"
                  UMAKE_OPTIONS=( "${UMAKE_OPTIONS[@]} -DTru64" )
                  ;;
               3)
                  CORRECT_OS="true"
                  export OS="lnx"
                  UMAKE_OPTIONS=( "${UMAKE_OPTIONS[@]} -DCompaqLinux" )
                  ;;
               *)
                  echo ' Please give an integer between 1 and 3'
               esac
            done
            ;;
       "2")
            CORRECT_PLATFORM="true"
            export MCH="cry"
            export MACHINE="CRAY-T3E"
            UMAKE_OPTIONS=( "${UMAKE_OPTIONS[@]} -D_CRAY" )
            export OS="unc"
            ;;
       "3")
            CORRECT_PLATFORM="true"
            export MCH="hp"
            export MACHINE="HP-workstation"
            UMAKE_OPTIONS=( "${UMAKE_OPTIONS[@]} -Dhpux" )

            CORRECT_OS="false"
            while [[  $CORRECT_OS == "false"  ]]; do 

               echo ' Select operating system'
               echo
               echo '   (1) HP-UX'
               echo '   (2) Linux'
            
               read CHOICE

               case  $CHOICE  in
               1)
                  CORRECT_OS="true"
                  export OS="hpu"
                  UMAKE_OPTIONS=( "${UMAKE_OPTIONS[@]} -DHpUx" )
                  ;;
               2)
                  CORRECT_OS="true"
                  export OS="lnx"
                  UMAKE_OPTIONS=( "${UMAKE_OPTIONS[@]} -DHpLinux" )
                  ;;
               *)
                  echo ' Please give an integer between 1 and 2'
               esac
            done
            ;;
       "4")
            CORRECT_PLATFORM="true"
            export MCH="rs6"
            export MACHINE="IBM-RS6000"
            UMAKE_OPTIONS=( "${UMAKE_OPTIONS[@]} -Dibm" )

            CORRECT_OS="false"
            while [[  $CORRECT_OS == "false"  ]]; do 

               echo ' Select operating system'
               echo
               echo '   (1) AIX'
               echo '   (2) Linux'
           
               read CHOICE

               case  $CHOICE  in
               1)
                  CORRECT_OS="true"
                  export OS="aix"
                  UMAKE_OPTIONS=( "${UMAKE_OPTIONS[@]} -DIbmAix" )
                  ;;
               2)
                  CORRECT_OS="true"
                  export OS="lnx"
                  UMAKE_OPTIONS=( "${UMAKE_OPTIONS[@]} -DIbmLinux" )
                  ;;
               *)
                  echo ' Please give an integer between 1 and 2'
               esac
            done
            ;;
       "5")
            CORRECT_PLATFORM="true"
            export MCH="pc"
            export MACHINE="Intel-like-PC"

            CORRECT_OS="false"
            while [ $CORRECT_OS == "false" ]; do

               echo ' Select operating system'
               echo
               echo '   (1) Windows 2000/XP with MinGW/Msys'
               echo '   (2) Linux'

               read CHOICE

               case  $CHOICE  in
               1)
                  CORRECT_OS="true"
                  export OS="mgw"
                  export OPSYS="MGW"
                  UMAKE_OPTIONS=( "${UMAKE_OPTIONS[@]} -DMinGW" )
                  UMAKE_CPP_OPTIONS=( "${UMAKE_CPP_OPTIONS[@]} -Ptraditional" )
                  ;;
               2)
                  CORRECT_OS="true"
                  export OS="lnx"
                  export OPSYS="Linux"
                  UMAKE_OPTIONS=( "${UMAKE_OPTIONS[@]} -Dlinux" )
                  # The following line is for Red Hat systems
                  UMAKE_CPP_OPTIONS=( "${UMAKE_CPP_OPTIONS[@]} -Ptraditional -Pw" )
                  ;;
               *)
                  echo ' Please give an integer between 1 and 2'
               esac
            done
            ;;
       "6")
            CORRECT_PLATFORM="true"
            export MCH="sgi"
            export MACHINE="SGI-workstation"
            UMAKE_OPTIONS=( "${UMAKE_OPTIONS[@]} -Dsgi" )

            CORRECT_OS="false"
            while [[  $CORRECT_OS == "false"  ]]; do 

               echo ' Select operating system'
               echo
               echo '   (1) IRIX'
               echo '   (2) Linux'
            
               read CHOICE

               case  $CHOICE  in
               1)
                  CORRECT_OS="true"
                  export OS="irx"
                  UMAKE_OPTIONS=( "${UMAKE_OPTIONS[@]} -DSgiIrix" )
                  ;;
               2)
                  CORRECT_OS="true"
                  export OS="lnx"
                  UMAKE_OPTIONS=( "${UMAKE_OPTIONS[@]} -DSgiLinux" )
                  ;;
               *)
                  echo ' Please give an integer between 1 and 2'
               esac
            done
            ;;
       "7")
            CORRECT_PLATFORM="true"
            export MCH="sun"
            export MACHINE="Sun-workstation"
            UMAKE_OPTIONS=( "${UMAKE_OPTIONS[@]} -Dsunq" )

            CORRECT_OS="false"
            while [[  $CORRECT_OS == "false"  ]]; do 

               echo ' Select operating system'
               echo
               echo '   (1) Solaris'
               echo '   (2) Linux'

               read CHOICE

               case  $CHOICE  in
               1)
                  CORRECT_OS="true"
                  export OS="sol"
                  UMAKE_OPTIONS=( "${UMAKE_OPTIONS[@]} -DSunSolaris" )
                  ;;
               2)
                  CORRECT_OS="true"
                  export OS="lnx"
                  UMAKE_OPTIONS=( "${UMAKE_OPTIONS[@]} -DSunLinux" )
                  ;;
               *)
                  echo ' Please give an integer between 1 and 2'
               esac
            done
            ;;
       "8")
            CORRECT_PLATFORM="true"
            export MCH="mac"
            export MACHINE="mac"
            export OS="osx"
            UMAKE_OPTIONS=( "${UMAKE_OPTIONS[@]} -Dmac" )
            ;;
       "0")
            CORRECT_PLATFORM="true"
            export MCH="gnx"
            export MACHINE="Generic Unix"
            export OS="unknown"
            ;;
       *)
            echo ' Please give an integer between 0 and 8'
   esac
done

# Get only the necessary platform-dependent commands
# to complete the installation. A more complete set
# will be generated after the user has chosen the
# compilers and precision they wish to use.

CPP='cpp -P'
EXT=$$
sed -f ${CUTER}/build/scripts/protect.sed ${CUTER}/build/scripts/makefile.cmds | ${CPP} -I${CUTER}/config -I. -DFromInstallCUTEr ${UMAKE_OPTIONS[@]} | sed -f ${CUTER}/build/scripts/gettabs.sed > /tmp/Makefile.cmds.${EXT}
cd /tmp
#[[  -f "cmds"  ]] && rm -f cmds
make -s -f Makefile.cmds.${EXT} EXT=.${EXT} 2> /dev/null
. /tmp/cmds.${EXT}
$RM Makefile.cmds.${EXT} cmds.${EXT}
cd ${CUTER}

echo ' Select Fortran compiler'

COMP=( `$SED -e '/#.*$/d' ${ARCH}/f.arch | $SED -e '/^ *$/d' | $GREP ${MCH} | $GREP ${OS} | $AWK 'BEGIN{FS=";"}{print $5}'` )
SYMB=( `$SED -e '/#.*$/d' ${ARCH}/f.arch | $SED -e '/^ *$/d' | $GREP ${MCH} | $GREP ${OS} | $AWK 'BEGIN{FS=";"}{print $4}'` )
CMPEXT=( `$SED -e '/#.*$/d' ${ARCH}/f.arch | $SED -e '/^ *$/d' | $GREP ${MCH} | $GREP ${OS} | $AWK 'BEGIN{FS=";"}{print $3}'` )
OPSYS=( `$SED -e '/#.*$/d' ${ARCH}/f.arch | $SED -e '/^ *$/d' | $GREP ${MCH} | $GREP ${OS} | $AWK 'BEGIN{FS=";"}{print $2}'` )
SYS=( `$SED -e '/#.*$/d' ${ARCH}/f.arch | $SED -e '/^ *$/d' | $GREP ${MCH} | $GREP ${OS} | $AWK 'BEGIN{FS=";"}{print $1}'` )
ALLCOMP=( `$SED -e '/#.*$/d' ${ARCH}/f.arch | $SED -e '/^ *$/d' | $GREP 'all' | $GREP 'all' | $AWK 'BEGIN{FS=";"}{print $5}'` )
ALLSYMB=( `$SED -e '/#.*$/d' ${ARCH}/f.arch | $SED -e '/^ *$/d' | $GREP 'all' | $GREP 'all' | $AWK 'BEGIN{FS=";"}{print $4}'` )
ALLCMPEXT=( `$SED -e '/#.*$/d' ${ARCH}/f.arch | $SED -e '/^ *$/d' | $GREP 'all' | $GREP 'all' | $AWK 'BEGIN{FS=";"}{print $3}'` )
ALLOPSYS=( `$SED -e '/#.*$/d' ${ARCH}/f.arch | $SED -e '/^ *$/d' | $GREP 'all' | $GREP 'all' | $AWK 'BEGIN{FS=";"}{print $2}'` )
ALLSYS=( `$SED -e '/#.*$/d' ${ARCH}/f.arch | $SED -e '/^ *$/d' | $GREP 'all' | $GREP 'all' | $AWK 'BEGIN{FS=";"}{print $1}'` )

COMP=(${COMP[@]} ${ALLCOMP[@]})
SYMB=(${SYMB[@]} ${ALLSYMB[@]})
CMPEXT=(${CMPEXT[@]}  ${ALLCMPEXT[@]})
OPSYS=(${OPSYS[@]} ${ALLOPSYS[@]})
SYS=(${SYS[@]} ${ALLSYS[@]})

let NUMBER=${#COMP[@]}

LIST=( ${COMP[@]} )
SYMBLIST=( ${SYMB[@]} )
CMPLIST=( ${CMPEXT[@]} )
OPSYSLIST=( ${OPSYS[@]} )
SYSLIST=( ${SYS[@]} )

CORRECT_COMPILER="false"
while [[  $CORRECT_COMPILER == "false"  ]]; do 

   let count=0
   while (( count < NUMBER ))
    do
        cmpInfo=`echo ${COMP[$count]} | tr '_' ' '`
        (( count++ ))
        printf '\t[%-2d]  %s\n' $count "$cmpInfo"
    done

   read CHOICE
   (( CHOICE-- ))

   let i=0
   while [[  ($i -lt $NUMBER) &&  ($CORRECT_COMPILER == "false") ]]; do
      if ((  CHOICE == i  )); then
        CORRECT_COMPILER="true"
        CHOICE=$i
      else
        (( i++ ))
      fi
   done
   if [[  $CORRECT_COMPILER == "true"  ]]; then
     COMPILER=${LIST[$CHOICE]}
     CMP=${CMPLIST[$CHOICE]}
     OSTAG=${OPSYSLIST[$CHOICE]}
     SYSTAG=${SYSLIST[$CHOICE]}
     UMAKE_OPTIONS=( "${UMAKE_OPTIONS[@]} ${SYMBLIST[$CHOICE]}" )
   else
     echo " Please give an integer between 1 and $NUMBER"
   fi
done

echo ' Select C compiler'

C_COMP=( `$SED -e '/#.*$/d' ${ARCH}/c.arch | $SED -e '/^ *$/d' | $GREP ${MCH} | $GREP ${OS} | $AWK 'BEGIN{FS=";"}{print $5}'` )
C_SYMB=( `$SED -e '/#.*$/d' ${ARCH}/c.arch | $SED -e '/^ *$/d' | $GREP ${MCH} | $GREP ${OS} | $AWK 'BEGIN{FS=";"}{print $4}'` )
C_CMPEXT=( `$SED -e '/#.*$/d' ${ARCH}/c.arch | $SED -e '/^ *$/d' | $GREP ${MCH} | $GREP ${OS} | $AWK 'BEGIN{FS=";"}{print $3}'` )
C_OPSYS=( `$SED -e '/#.*$/d' ${ARCH}/c.arch | $SED -e '/^ *$/d' | $GREP ${MCH} | $GREP ${OS} | $AWK 'BEGIN{FS=";"}{print $2}'` )
C_SYS=( `$SED -e '/#.*$/d' ${ARCH}/c.arch | $SED -e '/^ *$/d' | $GREP ${MCH} | $GREP ${OS} | $AWK 'BEGIN{FS=";"}{print $1}'` )
C_ALLCOMP=( `$SED -e '/#.*$/d' ${ARCH}/c.arch | $SED -e '/^ *$/d' | $GREP 'all' | $GREP 'all' | $AWK 'BEGIN{FS=";"}{print $5}'` )
C_ALLSYMB=( `$SED -e '/#.*$/d' ${ARCH}/c.arch | $SED -e '/^ *$/d' | $GREP 'all' | $GREP 'all' | $AWK 'BEGIN{FS=";"}{print $4}'` )
C_ALLCMPEXT=( `$SED -e '/#.*$/d' ${ARCH}/c.arch | $SED -e '/^ *$/d' | $GREP 'all' | $GREP 'all' | $AWK 'BEGIN{FS=";"}{print $3}'` )
C_ALLOPSYS=( `$SED -e '/#.*$/d' ${ARCH}/c.arch | $SED -e '/^ *$/d' | $GREP 'all' | $GREP 'all' | $AWK 'BEGIN{FS=";"}{print $2}'` )
C_ALLSYS=( `$SED -e '/#.*$/d' ${ARCH}/c.arch | $SED -e '/^ *$/d' | $GREP 'all' | $GREP 'all' | $AWK 'BEGIN{FS=";"}{print $1}'` )

C_COMP=( ${C_COMP[@]} ${C_ALLCOMP[@]} )
C_SYMB=( ${C_SYMB[@]} ${C_ALLSYMB[@]} )
C_CMPEXT=( ${C_CMPEXT[@]}  ${C_ALLCMPEXT[@]}  )
C_OPSYS=( ${C_OPSYS[@]} ${C_ALLOPSYS[@]} )
C_SYS=( ${C_SYS[@]} ${C_ALLSYS[@]} )

let NUMBER=${#C_COMP[@]}

C_LIST=( ${C_COMP[@]} )
C_SYMBLIST=( ${C_SYMB[@]} )
C_CMPLIST=( ${C_CMPEXT[@]} )
C_OPSYSLIST=( ${C_OPSYS[@]} )
C_SYSLIST=( ${C_SYS[@]} )

useCcomp='true'
CORRECT_COMPILER="false"
while [[  $CORRECT_COMPILER == "false"  ]]; do 

   let count=0
   printf '\t[%-2d]  %s\n' $count "No C extension to CUTEr"
   while (( count < NUMBER ))
   do
     (( count++ ))
     cmpInfo=`echo ${C_COMP[$count-1]} | tr '_' ' '`
     printf '\t[%-2d]  %s\n' $count "$cmpInfo"
   done

   read CHOICE
   
   let i=0
   while [[  (i -le NUMBER) &&  ($CORRECT_COMPILER == "false") ]]; do
      if ((  CHOICE == i  )); then
        CORRECT_COMPILER="true"
        let CHOICE=$i
      else
        (( i++ ))
      fi
   done
   if [[  ($CORRECT_COMPILER == "true") && (i -ne 0)  ]]; then
     ((CHOICE --))
     C_COMPILER=${C_LIST[$CHOICE]}
     CCMP=${C_CMPLIST[$CHOICE]}
     COSTAG=${C_OPSYSLIST[$CHOICE]}
     CSYSTAG=${C_SYSLIST[$CHOICE]}
     UMAKE_OPTIONS=( "${UMAKE_OPTIONS[@]} ${C_SYMBLIST[$CHOICE]}" )
   elif (( i == 0  )); then
     useCcomp='false'
   else
     echo " Please give an integer between 0 and $NUMBER"
   fi
done

if [[ $useCcomp == 'false' ]]; then
    echo "No C compiler"
else
    echo "Using C compiler: $C_COMPILER"
fi

if [[  $useCcomp == 'false'  ]]; then
  UMAKE_OPTIONS=( "${UMAKE_OPTIONS[@]} -DNoCcomp" )
else
  UMAKE_OPTIONS=( "${UMAKE_OPTIONS[@]} -UNoCcomp" )
fi

export INSTALLPRECISION=""
export INSTALLSIZE=""

while [[  $INSTALLPRECISION != "S" && $INSTALLPRECISION != "D"  ]]; do
    echo " Set install precision (D=double, s=single): "
    read INSTALLPRECISION
    [[  $INSTALLPRECISION == ""  ]] && export INSTALLPRECISION="D"
    export INSTALLPRECISION=`echo $INSTALLPRECISION | tr '[a-z]' '[A-Z]'`
done

case $INSTALLPRECISION in
    S) export PRECISION="single"
       PRC="CS   "
       PRC90='!S '
       ;;
    D) export PRECISION="double"
       PRC="CD   "
       PRC90='!D '
       ;;
    *) echo "$0: error occured while setting precision. Aborting."
       exit 1
esac

while [[ $INSTALLSIZE != "L" && $INSTALLSIZE != "M" && \
        $INSTALLSIZE != "S" && $INSTALLSIZE != "C" ]]
do
    echo " Set install size (L=large, m=medium, s=small, c=customized): "
    read INSTALLSIZE
    [[  $INSTALLSIZE == ""  ]] && export INSTALLSIZE="L"
    export INSTALLSIZE=`echo $INSTALLSIZE | tr '[a-z]' '[A-Z]'`
done

case $INSTALLSIZE in
    C) export SIZE="custom"
       SIZ="CCUS"
       SIZ90='!CUS'
       UMAKE_OPTIONS=( "${UMAKE_OPTIONS[@]} -DCustomSize" )
       ;;
    L) export SIZE="large"
       SIZ="CBIG"
       SIZ90='!BIG'
       UMAKE_OPTIONS=( "${UMAKE_OPTIONS[@]} -DLargeSize" )
       ;;
    M) export SIZE="medium"
       SIZ="CMED"
       SIZ90='!MED'
       UMAKE_OPTIONS=( "${UMAKE_OPTIONS[@]} -DMediumSize" )
       ;;
    S) export SIZE="small"
       SIZ="CTOY"
       SIZ90='!TOY'
       UMAKE_OPTIONS=( "${UMAKE_OPTIONS[@]} -DSmallSize" )
       ;;
    *) echo "$0: error occured while setting size. Aborting."
       exit 1
esac

#
# Now select the directory to install to
#

YESNO=""
let chosen=0

while (( chosen == 0 )); do

    MYCUTER=$CUTER/CUTEr.${SIZE}.${MCH}.${OS}.${CMP}

    echo
    echo ' By default, CUTEr with your selections will be installed in'
    echo "  $MYCUTER"
    echo ' Is this OK (Y/n)?'

    while [[  $YESNO != 'Y' && $YESNO != 'N'  ]]; do
        read YESNO
        [[  $YESNO == ""  ]] && YESNO="Y"
        YESNO=`echo $YESNO | tr '[a-z]' '[A-Z]'`
    done

    case  $YESNO  in
        N) echo " Enter alternative directory for the installation:"
           read MYCUTER
           [[ "$MYCUTER" == "$CUTER" ]] && echo '   Conflictual choice!'
           ;;
        Y) ;;
    esac

    [[ ! -d $MYCUTER ]] && $MKDIR $MYCUTER 2> /dev/null || true
    if (( $? != 0 )); then
        echo " Er... Unable to create directory $MYCUTER"
        YESNO=""
    else
        export MYCUTER
        echo " Using directory $MYCUTER"
        let chosen=1
    fi

done

# The following variable is defined in usage_mycuter.xxx.pro
UMAKE_SPECIAL='${x_umake_flags[@]}'
if [[  (-d $MYCUTER/single && $PRECISION == 'double') || ($PRECISION == 'single' &&  -d $MYCUTER/double)  ]]; then
    $LN -sf $CUTER/build/scripts/important $MYCUTER/IMPORTANT
    let important=1
fi

#
# Create a "working" directory for the version being installed
#

[[ ! -d $MYCUTER  ]] && $MKDIR $MYCUTER

export MYCUTERPREC=$MYCUTER/$PRECISION
if [[  -d $MYCUTERPREC  ]]; then
    echo " Directory"
    echo "  $MYCUTERPREC"
    echo " already exists. Proceeding will overwrite it." 
    echo " Do you wish to continue (Y/n)?"
    YESNO=""
    while [[  $YESNO != 'Y' && $YESNO != 'N'  ]]; do
        read YESNO
        [[  $YESNO == ""  ]] && YESNO="Y"
        YESNO=`echo $YESNO | tr '[a-z]' '[A-Z]'`
    done

    case  $YESNO  in
        Y) ;;
        N) echo "$0 : Aborting."
           exit 1
    esac
else
    $MKDIR $MYCUTERPREC
fi

#
# Check for a similar distribution already installed
#

cd $CUTER
EXISTS="0"
CUTERVERSION="$PRECISION $SIZE $MACHINE $OS $CMP"

echo ''
printf '     %-30s\t' "Logging operations"

if [[  -d log  ]]; then
   cd log
   for i  in  "*.log"; do
       EXISTS=`$GREP "$CUTERVERSION" $i`
   done
   if [[  $EXISTS != ""  ]]; then
      echo ''
      echo ' Hmmm... there already seems to be such a version installed in'
      echo $EXISTS | $AWK -F! '{print $3}'
      echo ' Are you sure you wish to continue (Y/n)?'
      YESNO=""
      while [[  $YESNO != "Y" && $YESNO != "N"  ]]; do
        read YESNO
        [[  $YESNO == ""  ]] && YESNO="Y"
        YESNO=`echo $YESNO | tr '[a-z]' '[A-Z]'`
      done

      case  $YESNO  in
         Y)
          ;;
         N)
          echo "$0 : aborting"
          exit 0
          ;;
         *)
          ;;
      esac
   fi
   cd ../
else
    $MKDIR log           # install log files
fi

printf '[Ok]\n'

#  Add a new entry to the log file

printf '     %-30s\t' "Creating directory structure"

echo  "`date` ! $CUTERVERSION ! $MYCUTER" >> $CUTER/log/install.log

[[ ! -d $MYCUTER/bin  ]] &&        $MKDIR $MYCUTER/bin
[[ ! -d $MYCUTERPREC/config  ]] && $MKDIR $MYCUTERPREC/config
[[ ! -d $MYCUTERPREC/cuter-specs  ]] &&  $MKDIR $MYCUTERPREC/cuter-specs
[[ ! -d $MYCUTERPREC/lib  ]] &&    $MKDIR $MYCUTERPREC/lib
[[ ! -d $MYCUTERPREC/bin  ]] &&    $MKDIR $MYCUTERPREC/bin
[[ ! -d $MYCUTER/config  ]] &&     $MKDIR $MYCUTER/config

printf '[Ok]\n'

# The following 'source' puts the timer routine in $CONFIG

. $ARCH/compiler.${SYSTAG}.${OSTAG}.${CMP}
. $ARCH/size.$SIZE

#  Define some useful directories

CONFIG=$MYCUTERPREC/config
BIN=$MYCUTERPREC/bin
SRC=$CUTER/common/src

#
#  Further initializations
#

set +C
DATE=`date`
echo "#version for ${MACHINE} under ${SYSTEM} with ${CMP} ($DATE)" > $CONFIG/this_build

printf '     %-30s\t' "Copying Umakefiles"

# Put the Umake template and rules files in place
$CP $CUTER/config/* $MYCUTER/config

# Put the Umakefiles in place
$CP $CUTER/build/scripts/umakefile.mycuterbin        $MYCUTER/bin/Umakefile
$CP $CUTER/build/scripts/umakefile.mycuter           $MYCUTER/Umakefile
$CP $CUTER/build/scripts/umakefile.mycuterprec       $MYCUTERPREC/Umakefile
$CP $CUTER/build/scripts/umakefile.mycuterprecbin    $MYCUTERPREC/bin/Umakefile
$CP $CUTER/build/scripts/umakefile.mycuterpreclib    $MYCUTERPREC/lib/Umakefile
$CP $CUTER/build/scripts/umakefile.mycuterprecspecs  $MYCUTERPREC/cuter-specs/Umakefile
$CP $CUTER/build/scripts/umakefile.mycuterprecconfig $MYCUTERPREC/config/Umakefile
$CP $CUTER/build/scripts/makefile.cmds               $MYCUTERPREC/config/makefile.cmds

printf '[Ok]\n'

UMAKE_SEQUENCE=( "$UMAKE_COMMAND ${UMAKE_OPTIONS[@]} ${UMAKE_CPP_OPTIONS[@]}" )

printf '     %-30s\t' "Creating installation shortcut"

echo '#!/bin/bash'                                   >  $MYCUTER/install_mycuter
echo ''                                              >> $MYCUTER/install_mycuter

# Must re-generate the platform-specific commands

EXT=$$
$SED -f ${CUTER}/build/scripts/protect.sed ${CUTER}/build/scripts/makefile.cmds | ${CPP} -I${CUTER}/config -I. -DFromInstallCUTEr ${UMAKE_OPTIONS[@]} | sed -f ${CUTER}/build/scripts/gettabs.sed > /tmp/Makefile.cmds.${EXT}
cd /tmp
#[[  -f cmds  ]] && rm -f cmds
make -s -f Makefile.cmds.${EXT} EXT=.${EXT} 2> /dev/null
$RM Makefile.cmds.${EXT}
cd ${CUTER}

# Flush these into install_mycuter

$CAT /tmp/cmds.${EXT}                                >> $MYCUTER/install_mycuter
if [[  ! -f $MYCUTER/cmds.basic  ]]; then
    $CP /tmp/cmds.${EXT} $MYCUTER/cmds.basic
    $CHMOD -w $MYCUTER/cmds.basic
fi
$RM /tmp/cmds.${EXT}

echo ''                                              >> $MYCUTER/install_mycuter
$CAT $CUTER/build/prototypes/usage_mycuter.sh.pro    >> $MYCUTER/install_mycuter

echo "${UMAKE_SEQUENCE[@]} ${UMAKE_SPECIAL}"         >> $MYCUTER/install_mycuter
echo '$MAKE Makefile'                                >> $MYCUTER/install_mycuter
echo '$MAKE Makefiles'                               >> $MYCUTER/install_mycuter
$CAT $CUTER/build/prototypes/run_mycuter.sh.pro      >> $MYCUTER/install_mycuter
$CHMOD +x $MYCUTER/install_mycuter

printf '[Ok]\n'

printf '     %-30s\t' "Installing README files"

if [[  ! (-d $MYCUTER/single && -d $MYCUTER/double)  ]]; then
    $LN -s $CUTER/build/scripts/readme.mycuterbin    $MYCUTER/bin/README
    $LN -s $CUTER/build/scripts/readme.mycuter       $MYCUTER/README
    $LN -s $CUTER/README                             $MYCUTER/README_FIRST
fi
$LN -s $CUTER/build/scripts/readme.mycuterprecconfig $MYCUTERPREC/config/README
$LN -s $CUTER/build/scripts/readme.mycuterpreclib    $MYCUTERPREC/lib/README
$LN -s $CUTER/build/scripts/readme.mycuterprecbin    $MYCUTERPREC/bin/README

if [[  "${CMP}" == 'ifc'  ]]; then
    $LN -fs $CUTER/build/scripts/readme.ifc           $MYCUTER/README.ifc
fi

printf '[Ok]\n'

echo ' -------------------------------------------------------'
echo ' '
echo " Please do not forget to define the environment variable"
echo "    [1] CUTER to be"
echo "          $CUTER"
echo "    [2] MYCUTER to be"
echo "          $MYCUTER"
echo ''
echo " as well as to set an appropriate value for MASTSIF,"
echo ' which should point to your main SIF repository.'
echo " "
echo "      In addition, please update your MANPATH to include "
echo "         $CUTER/common/man"
echo "      your PATH to include"
echo "         $MYCUTER/bin"
echo "      your MATLABPATH to include"
echo "         $MYCUTERPREC/bin"
echo "       and"
echo "         $CUTER/common/src/matlab"
echo "      and your LIBPATH (or LD_LIBRARY_PATH) to include"
echo "         $MYCUTERPREC/lib"
echo ''
echo " (see $CUTER/common/doc/readme.cshrc"
echo "  and $CUTER/common/doc/readme.bashrc"
echo "  for examples on how to do this)"
echo ''
echo " cd $MYCUTER,"
echo " read the README_FIRST and README files"
if (( important )); then
    echo ''
    echo '-----------------------------------'
    echo " ---and read the file IMPORTANT---"
    echo '-----------------------------------'
    echo ''
fi
echo " then run install_mycuter."

#
# See if the user wants to continue
#

echo ''
echo "install_cuter : Do you want install_mycuter to be run in"
echo " $MYCUTER now (Y/n)?"
YESNO=""
while [[  $YESNO != 'Y' && $YESNO != 'N'  ]]; do
    read YESNO
    [[  $YESNO == ""  ]] && YESNO="Y"
    YESNO=`echo $YESNO | tr '[a-z]' '[A-Z]'`
done

case  $YESNO  in
    Y) cd $MYCUTER
       if [[  "$PRECISION" == 'double'  ]]; then
           ./install_mycuter -double
       else
           ./install_mycuter -single
       fi
       ;;
    N) echo ''
       echo ' To complete the installation, change directory to'
       echo "    $MYCUTER"
       echo ' and run install_mycuter there.'
       echo ''
       exit 0
       ;;
esac

# This version serves as a template for the Umakefiles
# We exit here. The rest will be done by umake.

exit 0
