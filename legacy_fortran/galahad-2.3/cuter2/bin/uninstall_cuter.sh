#!/bin/bash

# Uninstall script for CUTEr
# N.Gould, D.Orban & Ph.Toint
# ( Last modified on 8 Feb 2001 at 20:11:40 )
#

# arg 1 should be the directory containing the CUTEr
# distribution to remove

# Check args first

name=`basename $0`

if [[  $# == 1 && ( $1 == '-h' || $1 == '--help' )  ]]; then
    echo "Usage: $name [directory]"
    echo "       If present, directory contains the distribution to remove."
    echo "       It should contain the *whole* path, e.g.:"
    echo "       /home/jack/CUTEr/CUTEr.large.alp.t64.f77"
    echo "       free of any tmp_mnt-like leading directory."
    echo ''
    echo "       Called with no argument, $name prompts you"
    echo "       with a choice."
    exit 1
fi

precIsDtm="false"
typeIsDtm="false"
lineIsDtm="false"

if [[  $# == 0  ]]; then

    let foundSystem=1
    if [[  `uname` == "SunOS" ]]; then
        UMAKE_OPTIONS=-Dsun
    elif [[  `uname` == "IRIX64"  ]]; then
        UMAKE_OPTIONS=-Dsgi
    elif [[  `uname` == "Linux" || `uname` == "linux"  ]]; then
        UMAKE_OPTIONS=-Dlinux
    elif [[ `uname` == 'Darwin' ]]; then
        UMAKE_OPTIONS=-Dmac
    elif [[  `uname` == "OSF1" ]]; then            # For both Compaq and Digital
        UMAKE_OPTIONS=-D__osf__
    elif [[  `uname` == "AIX"  ]]; then
        UMAKE_OPTIONS=-Dibm
    elif [[  `uname` == "HP"  ]]; then
        UMAKE_OPTIONS=-Dhpux
    elif [[  `uname` == "CRAY"  ]]; then
        UMAKE_OPTIONS=-D_CRAY
    else
        let foundSystem=0
        . ${CUTER}/build/arch/system.sh.all   # cross fingers
    fi

    if (( foundSystem == 1 )); then
        CPP='cpp -P'
        sed -f ${CUTER}/build/scripts/protect.sed ${CUTER}/build/scripts/makefile.cmds | ${CPP} -I${CUTER}/config -I. -DFromInstallCUTEr ${UMAKE_OPTIONS} | sed -f ${CUTER}/build/scripts/gettabs.sed > /tmp/Makefile.cmds.$$
        cd /tmp
        [[  -f "cmds"  ]] && rm -f cmds
        make -s -f Makefile.cmds.$$ cmds 2> /dev/null
        . /tmp/cmds
        cd ${CUTER}
        $RM /tmp/Makefile.cmds.$$ /tmp/cmds
    fi

    let nbselections=0
    if [[ ! -f ${CUTER}/log/install.log ]]; then
        let nblines=0
    else
        let nblines=`$WC -l ${CUTER}/log/install.log | $AWK '{print $1}'`
    fi

    # Lists initialisation
    let i=1
    dummyList=( "" )
    while (( i <= nblines )); do
        dummyList=( "${dummyList[@]} a" )
        (( i++ ))
    done
    dummyList=( ${dummyList[@]} )
    type=( ${dummyList[@]} )
    loc=( ${dummyList[@]} )
    line=( ${dummyList[@]} )

    let k=1
    let curline=1

    # Collect all possibilities
    # Pay attention to possible blank lines in install.log

    while (( curline <= nblines )); do
        if [[ ! -f ${CUTER}/log/install.log ]]; then
            contents=''
        else
            contents=`$HEAD -$curline ${CUTER}/log/install.log | $TAIL -1`
        fi
        if [[ "$contents" != "" ]]; then

            PREC=`echo "$contents" | $AWK -F! '{print $2}' | $AWK '{print $1}'`
            SIZE=`echo "$contents" | $AWK -F! '{print $2}' | $AWK '{print $2}'`
            MCH=`echo "$contents" | $AWK -F! '{print $2}' | $AWK '{print $3}'`
            OS=`echo "$contents" | $AWK -F! '{print $2}' | $AWK '{print $4}'`
            CMP=`echo "$contents" | $AWK -F! '{print $2}' | $AWK '{print $5}'`

            type[$k]="$PREC $SIZE $MCH $OS $CMP"
            loc[$k]=`echo "$contents" | $AWK -F! '{print $3}'`
            line[$k]=$curline
            (( k++ ))
        fi
        (( curline++ ))
    done

    (( nbselections = k - 1 ))
    correct_choice="false"

    if (( nbselections > 0 )); then
        while [[ $correct_choice == "false" ]]; do
            let k=1
            printf ' [ 0  ] : Abort\n'

            # Display possibilities

            while (( k <= nbselections )); do
                printf ' [ %-2d ] : %s,\n%s\n' $k "${type[$k]}" "          ${loc[$k]}"
                (( k++ ))
            done

            # Get user's choice

            read choice
            let i=0
            while [[ $i < $nbselections && $correct_choice == "false" ]]; do
                (( i++ ))
                [[ $choice == $i || $choice == "0" ]] && correct_choice="true"
            done
    
            if [[  $correct_choice == "true"  ]]; then
                if (( choice != 0 )); then
                    dir_to_erase=${loc[$choice]}
                    lines=${line[$choice]}
                    distrType="${type[$choice]}"

                    PREC=`echo ${type[$choice]} | $AWK '{print $1}'`
                    SIZE=`echo ${type[$choice]} | $AWK '{print $2}'`
                    MCH=`echo ${type[$choice]} | $AWK '{print $3}'`
                    OS=`echo ${type[$choice]} | $AWK '{print $4}'`
                    CMP=`echo ${type[$choice]} | $AWK '{print $5}'`

                    removeDir=${dir_to_erase}/$PREC
                    let isSingle=0
                    let isDouble=0
                    if [[ $PREC == "single" ]]; then
                        let isSingle=1
                        [[ -d ${dir_to_erase}/double ]] && let isDouble=1
                        newArg="-DDoublePrecision"
                    else
                        let isDouble=1
                        [[ -d ${dir_to_erase}/single ]] && let isSingle=1
                        newArg="-DSinglePrecision"
                    fi
                    precIsDtm="true"
                    typeIsDtm="true"
                    lineIsDtm="true"
                else
                    echo " Nothing was erased."
                    exit 2
                fi
            else
                echo " Please give an integer between 0 and $nbselections"
            fi
        done
    else
        echo " Could not find any instance installed."
        echo " Aborting."
        exit 0
    fi

else     # This is a call with at least one argument

    dir_to_erase=$1   # Ignore all other arguments

fi

dir_to_erase=`echo $dir_to_erase | sed -e 's/^[ ]+//;s/[ ]+$//'`

# Check if the specified directory exists

if [[ ! -d $dir_to_erase ]]; then
    echo "   $dir_to_erase "
    echo " does not seem to be a valid directory name."
    echo " Try re-running without argument."
    echo " Aborting."
    exit 1
fi

# Get rid of a possible trailing slash in the directory name

workOnDir=`echo $dir_to_erase | sed -e 's/\/$//'`

# Ask for which precision to remove

if [[ $precIsDtm == "false" ]]; then

    let isDouble=0
    let isSingle=0
    [[ -d "${workOnDir}/double" ]] && let isDouble=1
    [[ -d "${workOnDir}/single" ]] && let isSingle=1

    if (( isDouble && isSingle )); then
        echo " remove the single (s) or the double (d) distribution ?"
        read CHOICE
        CHOICE=`echo $CHOICE | tr a-z A-Z`
        while [[ $CHOICE != 'S' && $CHOICE != 'D' ]]; do
            echo 'I am expecting s or d'
            read CHOICE
            CHOICE=`echo $CHOICE | tr a-z A-Z`
        done

        # $MYCUTER/install_mycuter contains a '$1' argument to imake.
        # This needs to be updated according to the above choice.
        # -DSingleAndDouble will be replaced by $newArg

        if [[ $CHOICE == 'S' ]]; then
            PREC="single"
            newArg="-DDoublePrecision"
        else
            PREC="double"
            newArg="-DSinglePrecision"
        fi
    else
        if [[ $isSingle ]]; then
            PREC="single"
        elif [[ $isDouble ]]; then
            PREC="double"
        else
            echo "  Directory ${workOnDir} does not seem to contain a CUTEr distribution"
            #exit 1
        fi
    fi

    removeDir=${workOnDir}/$PREC

fi

# Get the UNIX commands right.
# If CUTEr was installed correctly, we shouldn't be in
# the 'else' case. But if we are, something went wrong
# and we don't know which system we're acting upon.

if [[ -e ${removeDir}/config/cmds ]]; then
    . ${removeDir}/config/cmds
else
    . ${CUTER}/build/arch/system.sh.all    # cross fingers
fi

if [[ $typeIsDtm == "false" ]]; then

# retrieve the size, machine, os and compiler name.
# For instance, MCH=Compaq-Alpha , OS=t64 and CMP=f77

    SIZE=`$GREP ${workOnDir} $CUTER/log/install.log | $GREP ${PREC} | "$AWK" -F! '{print $2}' | "$AWK" '{print $2}'`
    MCH=`$GREP ${workOnDir} $CUTER/log/install.log | $GREP ${PREC} | "$AWK" -F! '{print $2}' | "$AWK" '{print $3}'`
    OS=`$GREP ${workOnDir} $CUTER/log/install.log | $GREP ${PREC} | "$AWK" -F! '{print $2}' | "$AWK" '{print $4}'`
    CMP=`$GREP ${workOnDir} $CUTER/log/install.log | $GREP ${PREC} | "$AWK" -F! '{print $2}' | "$AWK" '{print $5}'`

    CUTERVERSION=\""${PREC} ${SIZE} ${MCH} ${OS} ${CMP}"\"

# If everything looks consistent, ask for confirmation
# and proceed

    distrType=`$GREP ${workOnDir} $CUTER/log/install.log | $GREP ${PREC} | awk -F! '{print $2}'`

fi

if [[ "${distrType}" == "" ]]; then
    echo 'There seems to be a conflict with your CUTEr environment variables'
    echo " CUTER = $CUTER"
    echo " MYCUTER  = $MYCUTER"
    exit 4
fi

echo " The $distrType CUTEr distribution stored in "
echo "  $removeDir "
echo " will be destroyed and lost."
echo " Do you wish to proceed (y/n)?"

read CHOICE
CHOICE=`echo $CHOICE | tr a-z A-Z`
while [[ $CHOICE != 'Y' && $CHOICE != 'N' ]]; do
      echo 'I am expecting y or n'
      read CHOICE
      CHOICE=`echo $CHOICE | tr a-z A-Z`
done

if [[ $CHOICE == 'N' ]]; then
 echo 'Nothing was erased'
 exit 0
fi

# Summon Mr. Clean

# We simply rm -r the specified directory
# and delete the corresponding entries in
# $CUTER/log/install.log

if [[ ! -d ${workOnDir}/bin ]]; then
    echo " Directory ${workOnDir} does not seem to contain"
    echo ' a valid/complete CUTEr instance...'
    echo ' Do you wish to remove this directory anyways?'
    read CHOICE
    CHOICE=`echo $CHOICE | tr a-z A-Z`
    while [[ $CHOICE != 'Y' && $CHOICE != 'N' ]]; do
      echo 'I am expecting y or n'
      read CHOICE
      CHOICE=`echo $CHOICE | tr a-z A-Z`
    done
    if [[ $CHOICE == 'N' ]]; then
    echo ' Aborting. Nothing was erased.'
    exit 1
    fi
fi

let cleanUp=1       # variable for debugging purposes

if (( cleanUp )); then

    $RM -r $removeDir

    # Unless single AND double precision were installed, we may
    # also remove the bin/ and ${workOnDir} directories,
    # otherwise, leave as is.

    if (( isSingle == 0 || isDouble == 0 )); then
        $RM -r ${workOnDir}
    else
        # Update $MYCUTER/install_mycuter if using the imake CUTEr

        if [[  -e ${workOnDir}/install_mycuter  ]]; then
            echo ' ...updating install_mycuter'
            script="s/"'$1'"/${newArg}/"
            $SED -e ${script} ${workOnDir}/install_mycuter > ${workOnDir}/install_mycuter.tmp
            $MV ${workOnDir}/install_mycuter.tmp ${workOnDir}/install_mycuter

            # Remove the 'important' message

            $RM ${workOnDir}/IMPORTANT
        fi
    fi
fi

if [[ $lineIsDtm == "false" ]]; then

    # Now clean the logfile, deleting lines according to their line number

    lines=`$GREP -n ${PREC} ${CUTER}/log/install.log | $GREP $SIZE | $GREP $MCH | $GREP $OS | $GREP $CMP | "$GREP" $workOnDir | $AWK -F: '{print $1}'`
    lines=`echo $lines | tr ' ' ','`

fi

if [[ $lines == "" ]]; then
    echo "Panic"
    exit 5
fi

$SED "$lines"d ${CUTER}/log/install.log > ${CUTER}/log/install.new.log

# If we are happy, summon Mr. Clean again

(( cleanUp )) && $MV $CUTER/log/install.new.log $CUTER/log/install.log

# End
exit 0
