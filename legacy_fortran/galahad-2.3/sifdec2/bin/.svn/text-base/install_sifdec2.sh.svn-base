#!/bin/bash
#
# Revised build script for SIFDEC
#
# syntax: install_sifdec
#
# N.Gould, D. Orban & Ph. Toint
#
###############################################################################

# Display synopsis

Display_Usage() {
    echo "Use: ${this}"
}

###############################################################################

# Initialize some variables and paths used during installation

Initialize() {
    SIFDEC=${PWD}
    ARCH=${SIFDEC}/build/arch

    UMAKE_COMMAND="${SIFDEC}/bin/umake.sh -I. -Iconfig"
    UMAKE_OPTIONS=( "" )
    UMAKE_CPP_OPTIONS=( "" )

    # Available platforms
    platform_list=( "Compaq Alpha" Cray "HP workstation" "IBM RS/6000" PC \
                    "SGI workstation" "SUN workstation" "MAC OS/X" Other )
    machines_list=( Compaq-Alpha CRAY-T3E HP-workstation IBM-RS6000   \
                    Intel-like-PC SGI-workstation Sun-workstation mac \
                    Generic-Unix )
    shortmch_list=( alp cry hp rs6 pc sgi sun mac gnx )
    machine_flags=( -D__osf__ -D_CRAY -Dhp -Dibm -Dpc -Dsgi -Dsun -Dmac -Dgen )

    # Choices of operating system if applicable
    alpha_os=( Linux Digital-Unix Tru64-Unix )
    alpha_os_nick=( lnx dux t64 )
    alpha_flags=( -DCompaqLinux -Ddux -DTru64 )
    
    cray_os=( Unicos )
    cray_os_nick=( unc )
    cray_flags=( -Dunc )
    
    hp_os=( Linux HPUX )
    hp_os_nick=( lnx hpu )
    hp_flags=( -DHpLinux -DHpUx )
    
    ibm_os=( Linux AIX )
    ibm_os_nick=( lnx aix )
    ibm_flags=( -DIbmLinux -DIbmAix )
    
    pc_os=( Linux MinGW )
    pc_os_nick=( lnx mgw )
    pc_flags=( -Dlinux -DMinGW )
    
    sgi_os=( Linux IRIX )
    sgi_os_nick=( lnx irx )
    sgi_flags=( -DSgiLinux -DSgiIrix )
    
    sun_os=( Linux Solaris )
    sun_os_nick=( lnx sol )
    sun_flags=( -DSunLinux -DSunSolaris )
    
    mac_os=( OS-X )
    mac_os_nick=( osx )
    mac_flags=( -Dosx )
    
    gen_os=( Generic-Unix )
    gen_os_nick=( gnx )
    gen_flags=( -Dgnx )
    
    os_list=( "${alpha_os[*]}" "${cray_os[*]}" "${hp_os[*]}" "${ibm_os[*]}" \
              "${pc_os[*]}" "${sgi_os[*]}" "${sun_os[*]}" "${mac_os[*]}"    \
              "${gen_os[*]}" )
    os_nicks=( "${alpha_os_nick[@]}" "${cray_os_nick[*]}" "${hp_os_nick[*]}" \
               "${ibm_os_nick[*]}" "${pc_os_nick[*]}" "${sgi_os_nick[*]}"    \
               "${sun_os_nick[*]}" "${mac_os_nick[*]}" "${gen_os_nick[*]}" )
    os_flags=( "${alpha_flags[@]}" "${cray_flags[*]}" "${hp_flags[*]}" \
               "${ibm_flags[*]}" "${pc_flags[*]}" "${sgi_flags[*]}"    \
               "${sun_flags[*]}" "${mac_flags[*]}" "${gen_flags[*]}" )
}

###############################################################################

# Prompt for a yes or no answer

Prompt_yesno() {
    yesno=''
    while [[  $yesno != 'Y' && $yesno != 'N'  ]]; do
        read yesno
        [[  $yesno == ""  ]] && YESNO="Y"
        yesno=`echo $yesno | tr '[a-z]' '[A-Z]'`
    done
}

###############################################################################

# Offer a selection
# The option selected is available as $REPLY and is value as $item

Make_a_Choice() {
    select item   # in list passed as argument
    do
        (( 1 <= REPLY && REPLY <= $# )) && break
    done
}

###############################################################################

this=`basename $0`

#  check input arguments (if any)
if [[  $# != 0  ]]; then
    Display_Usage
    exit 1
fi

# Initialize data
Initialize
new_flags=()

PS3='Select platform: '
Make_a_Choice "${platform_list[@]}"
printf "Temporary remark: item = $item, element number $REPLY\n"
(( REPLY-- ))
MCH=${shortmch_list[$REPLY]}
MACHINE=${machines_list[$REPLY]}
new_flags=( ${new_flags[@]} ${machine_flags[$REPLY]} )

PS3='Select operating system: '
applicable_os=( "${os_list[$REPLY]}" )
applicable_nicks=( "${os_nicks[$REPLY]}" )
applicable_flags=( "${os_flags[$REPLY]}" )
echo "applicable_os = ${applicable_os[@]}"
echo "applicable_nicks = ${applicable_nicks[@]}"
echo "applicable_flags = ${applicable_flags[@]}"

Make_a_Choice ${applicable_os[@]}
printf "Temporary remark: item = $item, element number $REPLY\n"
(( REPLY-- ))
OS=${applicable_nicks[$REPLY]}
new_flags=( ${new_flags[@]} ${applicable_flags[$REPLY]} )

echo "So far we have:"
echo "new_flags = ${new_flags[@]}"
echo "OS = $OS"
echo "MCH = $MCH"
echo "MACHINE = $MACHINE"
exit 0



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
            export OS="unc"
            UMAKE_OPTIONS=( "${UMAKE_OPTIONS[@]} -D_CRAY" )
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
            UMAKE_OPTIONS=( "${UMAKE_OPTIONS[@]} -Dsun" )

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
sed -f ${SIFDEC}/build/scripts/protect.sed ${SIFDEC}/build/scripts/makefile.cmds | ${CPP} -I${SIFDEC}/config -I. -DFromInstallSifDec ${UMAKE_OPTIONS[@]} | sed -f ${SIFDEC}/build/scripts/gettabs.sed > /tmp/Makefile.cmds.$$
cd /tmp
[[  -f "cmds"  ]] && rm -f cmds
make -s -f Makefile.cmds.$$ cmds 2> /dev/null
#unalias source
source /tmp/cmds
$RM Makefile.cmds.$$ cmds slash.sed
cd ${SIFDEC}

echo ' Select Fortran compiler'

COMP=(`$SED -e '/#.*$/d' ${ARCH}/f.arch | $SED -e '/^ *$/d' | $GREP ${MCH} | $GREP ${OS} | $AWK 'BEGIN{FS=";"}{print $5}'`)
SYMB=(`$SED -e '/#.*$/d' ${ARCH}/f.arch | $SED -e '/^ *$/d' | $GREP ${MCH} | $GREP ${OS} | $AWK 'BEGIN{FS=";"}{print $4}'`)
CMPEXT=(`$SED -e '/#.*$/d' ${ARCH}/f.arch | $SED -e '/^ *$/d' | $GREP ${MCH} | $GREP ${OS} | $AWK 'BEGIN{FS=";"}{print $3}'`)
OPSYS=(`$SED -e '/#.*$/d' ${ARCH}/f.arch | $SED -e '/^ *$/d' | $GREP ${MCH} | $GREP ${OS} | $AWK 'BEGIN{FS=";"}{print $2}'`)
SYS=(`$SED -e '/#.*$/d' ${ARCH}/f.arch | $SED -e '/^ *$/d' | $GREP ${MCH} | $GREP ${OS} | $AWK 'BEGIN{FS=";"}{print $1}'`)
ALLCOMP=(`$SED -e '/#.*$/d' ${ARCH}/f.arch | $SED -e '/^ *$/d' | $GREP 'all' | $GREP 'all' | $AWK 'BEGIN{FS=";"}{print $5}'`)
ALLSYMB=(`$SED -e '/#.*$/d' ${ARCH}/f.arch | $SED -e '/^ *$/d' | $GREP 'all' | $GREP 'all' | $AWK 'BEGIN{FS=";"}{print $4}'`)
ALLCMPEXT=(`$SED -e '/#.*$/d' ${ARCH}/f.arch | $SED -e '/^ *$/d' | $GREP 'all' | $GREP 'all' | $AWK 'BEGIN{FS=";"}{print $3}'`)
ALLOPSYS=(`$SED -e '/#.*$/d' ${ARCH}/f.arch | $SED -e '/^ *$/d' | $GREP 'all' | $GREP 'all' | $AWK 'BEGIN{FS=";"}{print $2}'`)
ALLSYS=(`$SED -e '/#.*$/d' ${ARCH}/f.arch | $SED -e '/^ *$/d' | $GREP 'all' | $GREP 'all' | $AWK 'BEGIN{FS=";"}{print $1}'`)

#echo "COMP = ${COMP[@]}"
#echo "ALLCOMP = ${ALLCOMP[@]}"

COMP=(${COMP[@]} ${ALLCOMP[@]})
SYMB=(${SYMB[@]} ${ALLSYMB[@]})
CMPEXT=(${CMPEXT[@]}  ${ALLCMPEXT[@]})
OPSYS=(${OPSYS[@]} ${ALLOPSYS[@]})
SYS=(${SYS[@]} ${ALLSYS[@]})

#echo "SYMB = ${SYMB[@]}"
#echo "CMPEXT = ${CMPEXT[@]}"
#echo "OPSYS = ${OPSYS[@]}"
#echo "SYS = ${SYS[@]}"

let NUMBER=${#COMP[@]}
#echo "Found $NUMBER compilers"
#exit 0

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
   while [[  ($i -lt $NUMBER) &&  ($CORRECT_COMPILER == "false")  ]]; do
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

export INSTALLPRECISION=""
export INSTALLSIZE=""

while [[  $INSTALLPRECISION != "S" && $INSTALLPRECISION != "D"  ]]; do
    echo " Select install precision (D=double, s=single): "
    read INSTALLPRECISION
    [[  $INSTALLPRECISION == ""  ]] && export INSTALLPRECISION="D"
    export INSTALLPRECISION=`echo $INSTALLPRECISION | tr '[a-z]' '[A-Z]'`
done

case $INSTALLPRECISION in
    S)
        export PRECISION="single"
        PRC="CS   "
        PRC90='!S '
    ;;
    D)
        export PRECISION="double"
        PRC="CD   "
        PRC90='!D '
    ;;
    *)
        echo "$0: error occured while setting precision. Aborting."
        exit 1
esac

while [[ $INSTALLSIZE != "L" && $INSTALLSIZE != "M" && \
        $INSTALLSIZE != "S" && $INSTALLSIZE != "C" ]]; do
    echo " Select install size (L=large, M=medium, S=small, C=customized): "
    read INSTALLSIZE
    [[  $INSTALLSIZE == ""  ]] && export INSTALLSIZE="L"
    export INSTALLSIZE=`echo $INSTALLSIZE | tr '[a-z]' '[A-Z]'`
done

case $INSTALLSIZE in
    C)
        export SIZE="custom"
        SIZ="CCUS"
        SIZ90='!CUS'
        UMAKE_OPTIONS=( "${UMAKE_OPTIONS[@]} -DCustomSize" )
    ;;
    L)
        export SIZE="large"
        SIZ="CBIG"
        SIZ90='!BIG'
        UMAKE_OPTIONS=( "${UMAKE_OPTIONS[@]} -DLargeSize" )
    ;;
    M)
        export SIZE="medium"
        SIZ="CMED"
        SIZ90='!MED'
        UMAKE_OPTIONS=( "${UMAKE_OPTIONS[@]} -DMediumSize" )
    ;;
    S)
        export SIZE="small"
        SIZ="CTOY"
        SIZ90='!TOY'
        UMAKE_OPTIONS=( "${UMAKE_OPTIONS[@]} -DSmallSize" )
    ;;
    *)
        echo "$0: error occured while setting size. Aborting."
        exit 1
esac

#
# Now select the directory to install to
#
# Default install directories have the form
# SifDec.large.sun.sol.f77 for large version on a
# Sun machine running Solaris.
#

YESNO=""
let chosen=0

while (( chosen == 0 )); do

    MYSIFDEC=$SIFDEC/SifDec.${SIZE}.${MCH}.${OS}.${CMP}

    echo
    echo ' By default, SifDec with your selections will be installed in'
    echo "  $MYSIFDEC"
    echo ' Is this OK (Y/n)?'

    while [[  $YESNO != 'Y' && $YESNO != 'N'  ]]; do
        read YESNO
        [[  $YESNO == ""  ]] && YESNO="Y"
        YESNO=`echo $YESNO | tr '[a-z]' '[A-Z]'`
    done

    case  $YESNO  in
        N) 
            echo " Enter alternative directory for the installation:"
            read MYSIFDEC
            [[ "$MYSIFDEC" == "$SIFDEC" ]] && echo '   Conflictual choice!'
            ;;
        Y) 
            ;;
    esac

    [[ ! -d $MYSIFDEC ]] && $MKDIR $MYSIFDEC 2> /dev/null || true
    if (( $? != 0 )); then
        echo " Er... Unable to create directory $MYSIFDEC"
        YESNO=""
    else
        export MYSIFDEC
        echo " Using directory $MYSIFDEC"
        let chosen=1
    fi

done

# The following variable is defined in usage_mysifdec.xxx.pro
#UMAKE_SPECIAL=''
UMAKE_SPECIAL='${x_umake_flags[@]}'
if [[  (-d $MYSIFDEC/single && $PRECISION == 'double') || ($PRECISION == 'single' &&  -d $MYSIFDEC/double)  ]]; then
#    UMAKE_SPECIAL='$1'
    $LN -sf $SIFDEC/build/scripts/important $MYSIFDEC/IMPORTANT
    let important=1
#elif [[  "$PRECISION" == 'double'  ]]; then
#    UMAKE_OPTIONS=( "${UMAKE_OPTIONS[@]} -DDoublePrecision" )
#else
#    UMAKE_OPTIONS=( "${UMAKE_OPTIONS[@]} -DSinglePrecision" )
fi

#
# Create a "working" directory for the version being installed
#

export MYSIFDECPREC=$MYSIFDEC/$PRECISION
if [[  -d $MYSIFDECPREC  ]]; then
    echo " Directory"
    echo "  $MYSIFDECPREC"
    echo " already exists. Proceeding will overwrite it." 
    echo " Do you wish to continue (Y/n)?"
    YESNO=""
    while [[  $YESNO != 'Y' && $YESNO != 'N'  ]]; do
        read YESNO
        [[  $YESNO == ""  ]] && YESNO="Y"
        YESNO=`echo $YESNO | tr '[a-z]' '[A-Z]'`
    done

    case  $YESNO  in
    Y) 
        ;;
    N) 
        echo "$0 : Aborting."
        exit 1
    esac
else
    $MKDIR $MYSIFDECPREC
fi

#
# Check for a similar distribution already installed
#

cd $SIFDEC
EXISTS="0"
SIFDECVERSION="$PRECISION $SIZE $MACHINE $OS $CMP"

echo ''
printf '     %-30s\t' "Logging operations"

if [[  -d log  ]]; then
   cd log
   for i  in  "*.log"; do
       EXISTS=`$GREP "$SIFDECVERSION" $i`
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

echo  "`date` ! $SIFDECVERSION ! $MYSIFDEC" >> $SIFDEC/log/install.log

[[ ! -d $MYSIFDEC/bin  ]] &&        $MKDIR $MYSIFDEC/bin
[[ ! -d $MYSIFDECPREC/config  ]] && $MKDIR $MYSIFDECPREC/config
[[ ! -d $MYSIFDECPREC/bin  ]] &&    $MKDIR $MYSIFDECPREC/bin
[[ ! -d $MYSIFDEC/config  ]] &&     $MKDIR $MYSIFDEC/config

printf '[Ok]\n'

# The following 'source' puts the timer routine in $CONFIG

#. $ARCH/compiler.${SYSTAG}.${OSTAG}.${CMP}
. $ARCH/size.$SIZE

#  Define some useful directories

CONFIG=$MYSIFDECPREC/config
BIN=$MYSIFDECPREC/bin
SRC=$SIFDEC/common/src

#
#  Further initializations
#

set +C
DATE=`date`
echo "#version for ${MACHINE} under ${SYSTEM} with ${CMP} ($DATE)" > $CONFIG/this_build

printf '     %-30s\t' "Copying Umakefiles"

# Put the Umake template and rules files in place
$CP $SIFDEC/config/* $MYSIFDEC/config

# Put the Umakefiles in place
$CP $SIFDEC/build/scripts/umakefile.mysifdecbin        $MYSIFDEC/bin/Umakefile
$CP $SIFDEC/build/scripts/umakefile.mysifdec           $MYSIFDEC/Umakefile
$CP $SIFDEC/build/scripts/umakefile.mysifdecprec       $MYSIFDECPREC/Umakefile
$CP $SIFDEC/build/scripts/umakefile.mysifdecprecbin    $MYSIFDECPREC/bin/Umakefile
$CP $SIFDEC/build/scripts/umakefile.mysifdecprecconfig $MYSIFDECPREC/config/Umakefile
$CP $SIFDEC/build/scripts/makefile.cmds                $MYSIFDECPREC/config/makefile.cmds

printf '[Ok]\n'

UMAKE_SEQUENCE=( "$UMAKE_COMMAND ${UMAKE_OPTIONS[@]} ${UMAKE_CPP_OPTIONS[@]}" )

printf '     %-30s\t' "Creating installation shortcut"

echo '#!/bin/bash'                                >  $MYSIFDEC/install_mysifdec
echo ''                                           >> $MYSIFDEC/install_mysifdec

# Re-generate the platform-specific commands if necessary

$SED -f ${SIFDEC}/build/scripts/protect.sed ${SIFDEC}/build/scripts/makefile.cmds | ${CPP} -I${SIFDEC}/config -I. -DFromInstallSifDec ${UMAKE_OPTIONS[@]} | sed -f ${SIFDEC}/build/scripts/gettabs.sed > /tmp/Makefile.cmds.$$
cd /tmp
[[  -f cmds  ]] && rm -f cmds
make -s -f Makefile.cmds.$$ cmds 2> /dev/null
$RM Makefile.cmds.$$ slash.sed
cd ${SIFDEC}

# Flush these into install_mysifdec

$CAT /tmp/cmds                                    >> $MYSIFDEC/install_mysifdec
$CP -f /tmp/cmds $MYSIFDEC/cmds.basic
$CHMOD -w $MYSIFDEC/cmds.basic
$RM /tmp/cmds


echo ''                                           >> $MYSIFDEC/install_mysifdec
$CAT $SIFDEC/build/prototypes/usage_mysifdec.sh.pro >> $MYSIFDEC/install_mysifdec
echo "${UMAKE_SEQUENCE[@]} $UMAKE_SPECIAL"        >> $MYSIFDEC/install_mysifdec
echo 'make Makefile'                              >> $MYSIFDEC/install_mysifdec
echo 'make Makefiles'                             >> $MYSIFDEC/install_mysifdec
$CAT $SIFDEC/build/prototypes/run_mysifdec.sh.pro >> $MYSIFDEC/install_mysifdec
chmod +x $MYSIFDEC/install_mysifdec

printf '[Ok]\n'

printf '     %-30s\t' "Installing README files"

if [[  ! (-d $MYSIFDEC/single && -d $MYSIFDEC/double)  ]]; then
    $LN -sf $SIFDEC/build/scripts/readme.mysifdecbin    $MYSIFDEC/bin/README
    $LN -sf $SIFDEC/build/scripts/readme.mysifdec       $MYSIFDEC/README
    $LN -sf $SIFDEC/README                              $MYSIFDEC/README_FIRST
fi
$LN -sf $SIFDEC/build/scripts/readme.mysifdecprecconfig $MYSIFDECPREC/config/README
$LN -sf $SIFDEC/build/scripts/readme.mysifdecprecbin    $MYSIFDECPREC/bin/README

printf '[Ok]\n'

echo ''
echo " Please do not forget to define the environment variable"
echo "    [1] SIFDEC to be"
echo "          $SIFDEC"
echo "    [2] MYSIFDEC to be"
echo "          $MYSIFDEC"
echo ''
echo " as well as setting an appropriate value for MASTSIF,"
echo ' which should point to your main SIF repository.'
echo ''
echo "      In addition, please update your MANPATH to include "
echo "         $SIFDEC/common/man"
echo "      and your PATH to include"
echo "         $MYSIFDEC/bin"
echo ''
echo " (see $SIFDEC/common/doc/readme.cshrc"
echo "  and $SIFDEC/common/doc/readme.bashrc"
echo "  for examples on how to do this)"
echo ''
echo " cd $MYSIFDEC,"
echo " read the README_FIRST and README files"
if (( important != 0  )); then
    echo ''
    echo '-----------------------------------'
    echo "--- and read the file IMPORTANT ---"
    echo '-----------------------------------'
    echo ''
fi
echo " then run install_mysifdec."

#
# See if the user wants to continue
#

echo ''
echo "install_sifdec : Do you want install_mysifdec to be run in"
echo " $MYSIFDEC now (Y/n)?"
YESNO=""
while [[  $YESNO != 'Y' && $YESNO != 'N'  ]]; do
    read YESNO
    [[  $YESNO == ""  ]] && YESNO="Y"
    YESNO=`echo $YESNO | tr '[a-z]' '[A-Z]'`
done

case  $YESNO  in
    Y) 
    cd $MYSIFDEC
    #if ((  important != 0 )); then
        if [[  "$PRECISION" == 'double'  ]]; then
            ./install_mysifdec -double
        else
            ./install_mysifdec -single
        fi
    #else
    #    ./install_mysifdec
    #fi
    ;;
    N)
    echo ''
    echo ' To complete the installation, change directory to'
    echo "    $MYSIFDEC"
    echo ' and run install_mysifdec there.'
    echo ''
    exit 0
    ;;
esac

# This version serves as a template for the Umakefiles
# We exit here. The rest will be done by umake.

exit 0

# N. Gould, D. Orban and Ph. Toint, dec 5, 2001.

