#!/bin/bash
# sifdecode: script to decode a sif file
#  ( Last modified on 23 Dec 2000 at 17:29:56 )
#
#{version}
#
# Use: sifdecode [-s] [-h] [-k] [-o j] [-l secs] [-f] [-b] [-a j] [-show] [-param name=value[,name=value...]] [-force] [-debug] probname[.SIF]
#
# where: options -s     : run the single precision version
#                         (Default: run the double precision version)
#                -h     : print this help and stop execution
#                -k     : keep the load module after use
#                         (Default: delete the load module)
#                -o     : 0 for silent mode, 1 for brief description of
#                         the stages executed.
#                         (Default: -o 0)
#                -l     : sets a limit of secs second on the runtime
#                         (Default: 99999999 seconds)
#                -f     : use automatic differentiation in Forward mode
#                -b     : use automatic differentiation in Backward mode
#                -a     : 1 use the older HSL automatic differentiation package AD01
#                         2 use the newer HSL automatic differentiation package AD02
#                         (Default: -a 2)
#                -show  : displays possible parameter settings for
#                         probname[.SIF]. Other options are ignored
#                -param : cast probname[.SIF] against explicit parameter
#                         settings. Several parameter settings may be
#                         given as a comma-separated list following
#                         -param or using several -param flags.
#                         Use -show to view possible settings.
#                -force : forces setting of a parameter to the given value,
#                         even if this value is not specified in the file.
#                         This option should be used with care.
#                         (Default: do not enforce).
#                -debug : links all the libraries, creates the executable
#                         and stop to allow debugging. This option
#                         automatically disables -n and enables -k.
#
#       probname      probname.SIF is the name of the file containing
#                     the SIF file for the problem of interest.
#

#
#  N. Gould, D. Orban & Ph. Toint, November 7th, 2000
#

#{cmds}

set +C

#
# Environment check
#

let abort=0

if [[ ${SIFDEC+set} != 'set'  ]]; then
    echo ' SIFDEC is not set.'
    echo ' It should point to the directory where you installed SifDec.'
    echo ' Set it to the appropriate value and re-run.'
    let abort=1
fi

if [[ ${MYSIFDEC+set} != 'set'  ]]; then
    echo ' MYSIFDEC is not set.'
    echo ' It should point to your working instance of SifDec.'
    echo ' Set it to the appropriate value and re-run.'
    echo ''
    let abort=1
fi

if (( abort != 0 )); then
    echo ' Aborting.'
    exit 7
fi

#
#  variables for each option
#

#
#  LIMIT (maximum cputime)
#

let LIMIT=99999999

#
# PRECISION = single (single precision), = double (double precision)
#

PRECISION="double"
let prec=1

#
# KEEP = (null) (discard load module after use), = "-k" (keep it)
#

KEEP=""

#
# DEBUGARG = "" (normal execution)
#

DEBUGARG=""

#
# OUTPUT = 0 (summary output), = 1 (detailed output from decoder)
#

let OUTPUT=0

#
#   automatic = 0 (provided), = 1 (automatic forward), = 2 (automatic backward)
#

let automatic=0

#
#   ad0 = 1 (AD01 used), = 2 (AD02 used)
#

let ad0=2

#
# FORCE = 0 (do not enforce -param settings), = 1 (enforce settings)
#

let FORCE=0

#
#  interpret arguments
#

let METHOD=3

#
#  Check whether variable PAC has been set;
#  this allows sifdecode to be called from the command line.
#

if [[ $PAC == "" ]]; then
    export PAC=sifdecode
    export caller=sifdecode
    thiscaller=`basename $0`
else
    thiscaller=sd${PAC}
fi

let last=$#

if (( last == 0 )); then
  echo "Use: $thiscaller [-s] [-h] [-k] [-o j] [-l secs] [-f] [-b] [-a j] [-show] [-param name=value[,name=value...]] [-debug] probname[.SIF]"
  exit 1
fi

let i=1
let nbparam=0   # counts the number of parameters passed using -param
PARAMLIST=( "" )

while (( i <= last )); do
  opt=${!i}
  if [[ "$opt" == '-h' || "$opt" == '--help' ]]; then
    $MYSIFDEC/bin/helpmsg
    exit 0
  elif [[ "$opt" == '-s' ]]; then
    PRECISION="single"
    let prec=0
  elif [[ "$opt" == '-k' ]]; then
    KEEP="-k"
  elif [[ "$opt" == '-o' ]]; then
    (( i++ ))
    let OUTPUT=${!i}
  elif [[ "$opt" == '-l' ]]; then
    (( i++ ))
    let LIMIT=${!i}
  elif [[ "$opt" == '-m' ]]; then
    (( i++ ))
    let METHOD=${!i}
  elif [[ "$opt" == '-f' ]]; then
    let automatic=1
  elif [[ "$opt" == '-b' ]]; then
    let automatic=2
  elif [[ "$opt" == '-a' ]]; then
    (( i++ ))
    if [[  ${!i} == '1' || ${!i} == '2'  ]]; then
        let ad0=${!i}
        else
        echo " ${thiscaller} : error processing -a flag"
    exit 6
    fi
  elif [[ "$opt" == '-param' ]]; then
    # parameters can either be separated by -param flags
    # or by commas within a -param flag, e.g.:
    # -options... -param par1=val1,par2=val2 -otheroptions... -param par3=val3
    (( i++ ))
    # see if some parameters were given as a comma-separated list
    let nbparInThisGroup=`echo ${!i} | $AWK -F, '{print NF}'`
    if [[ "$nbparInThisGroup" == "0" ]]; then
        echo " ${thiscaller} : error processing -param flag"
        exit 5
    fi
    # check if what follows -param looks like a parameter setting
    chkparam=`echo ${!i} | $GREP =`
    if [[  $chkparam == ""  ]]; then
        echo " ${thiscaller} : error processing -param flag"
        exit 5
    fi
    # store the parameters in the PARAM array
    PARAMLIST=( "$PARAMLIST `echo ${!i} | $SED -e 's/,/ /g'`" )
    (( nbparam = nbparam + nbparInThisGroup ))
  elif [[ "$opt" == '-show' ]]; then
    let SHOWPARAMS=1
    (( i++ ))
  elif [[ "$opt" == '-force' ]]; then
    let FORCE=1
  elif [[ "$opt" == '-debug' ]]; then
    KEEP="-k"
    DEBUGARG="-debug"
  fi
  (( i++ ))
done

# Create problem names

let last=$#
PROBLEM=${!last}
PROBNAME=`echo $PROBLEM | $AWK -F/ '{print $NF}' | $AWK -F. '{print $1}'`
EXT=`echo $PROBLEM | $AWK -F/ '{print $NF}' | $AWK -F. '{print $2}'`

PROBDIR=`echo $PROBLEM | $AWK 'BEGIN{FS="\n"}{ l = split( $1, farray, "/" ); head = ""; if( l>1 ) { head = farray[1]; for ( i=2; i<l; i++ ) { head = head "/" farray[i]; } } print head;}'`

probDirNotGiven="false"

if [[  "$PROBDIR" == ""  ]]; then
    probDirNotGiven="true"
    PROBDIR='.'
fi


problemLoc=${PROBDIR}/${PROBNAME}.SIF

# Specify correct path for problemLoc

let lookInMastSif=0
if [[ (! -e "$problemLoc") && ("$probDirNotGiven" == "true") && (${MASTSIF+set} == 'set') ]]; then
    let lookInMastSif=1
    problemLoc=${MASTSIF}/${PROBNAME}.SIF
fi
#echo "Using source file ${problemLoc}"

# See whether the specified SIF file exists

if [[ ! -e "$problemLoc" ]]; then
  if [[ ${MASTSIF+set} != 'set' ]]; then
    echo ' MASTSIF is not set'
  fi
  echo "file $PROBNAME.SIF is not known in directories"
  [[  $PROBDIR == "$PROBLEM" ]] && PROBDIR=$PWD
  echo " $PROBDIR or "'$MASTSIF'
  exit 2
fi

# See if the -show flag is present

if (( SHOWPARAMS != 0 )); then
    SHOWDOTAWK=${MYSIFDEC}/bin/show.awk
    if (( OUTPUT == 1 )); then
    echo "possible parameter settings for $PROBNAME are:"
    fi
    $GREP '$-PARAMETER' ${problemLoc} | $AWK -f $SHOWDOTAWK
    exit 8
# Note: the 'exit 8' command is used if -show is given on the command-line
# of a sd* interface. If the exit code were 0, the interface would go on,
# invoking the solver it interfaces.
fi

# Check if -param arguments have been passed and process them

if ((  nbparam > 0  )); then
    # part parameter names from their value
    let p=0
    PARLIST=( ${PARAMLIST[@]} ) #`echo $PARAMLIST:q`
    PARNAMELIST=( ${PARAMLIST[@]} ) #`echo $PARAMLIST:q`
    PARVALUELIST=( ${PARAMLIST[@]} ) #`echo $PARAMLIST:q`
    while ((  p < nbparam  )); do
        PARNAMELIST[$p]=`echo ${PARLIST[$p]} | $AWK -F= '{print $1}'`
        PARVALUELIST[$p]=`echo ${PARLIST[$p]} | $AWK -F= '{print $2}'`
        (( p++ ))
    done

    PARAMDOTAWK=${MYSIFDEC}/bin/param.awk
    # substitute the chosen parameter values in the SIF file
    let p=0

    # the easiest looks like creating a sed script which operates
    # all the necessary changes at once. We fill this sed script
    # as each parameter setting is examined in turn.

    # I give it a name depending on the current pid
    # hopefully, this is a unique name.

    sedScript=$TMP/$$_${PROBNAME}.sed
    sifFile=$TMP/$$_${PROBNAME}.SIF
    echo '' > $sedScript

    while ((  p < nbparam  )); do
        # see if parameter number p can be set to value number p
        # if so, retrieve the number of the matching line in the file
        # and the numbers of the lines which should be commented out.
        matchingLines="0"
        nonMatchingLines="0"

        matchingLines=`$GREP -n '$-PARAMETER' ${problemLoc} | $AWK -F'$' '{print $1}' | $AWK -v pname=${PARNAMELIST[$p]} -v pval=${PARVALUELIST[$p]} -v doesmatch=1 -f ${PARAMDOTAWK}`

    nonMatchingLines=`$GREP -n '$-PARAMETER' ${problemLoc} | $AWK -F'$' '{print $1}' | $AWK -v pname=${PARNAMELIST[$p]} -v pval=${PARVALUELIST[$p]} -v doesmatch=0 -f ${PARAMDOTAWK}`

    if [[  "$nonMatchingLines" != ""  ]]; then
        if (( OUTPUT == 1 )); then
            echo "lines number $nonMatchingLines will be commented out"
        fi
        for l  in  $nonMatchingLines; do
            echo "$l s/^ /\*/" >> $sedScript
        done
    fi

    let failed=0

    if [[  "$matchingLines" == ""  ]]; then
        if (( FORCE == 0 )); then
            echo " ${thiscaller} : Failed setting ${PARNAMELIST[$p]} to ${PARVALUELIST[$p]} -- skipping"
            let failed=1
        else
            # get the number of the first line defining the parameter
            # and the parameter type (IE, RE, ...)
            fline=`$GREP -n '$-PARAMETER' ${problemLoc} | $GREP ${PARNAMELIST[$p]} | $HEAD -1 | $AWK -F: '{print $1}'`
            type=`$GREP '$-PARAMETER' ${problemLoc} | $GREP ${PARNAMELIST[$p]} | $HEAD -1 | $AWK '{print $1}' | $SED -e 's/^[ ]*//'`
            (( OUTPUT == 1 )) && echo "Forcing parameter ${PARNAMELIST[$p]} on line $fline of type $type to value ${PARVALUELIST[$p]}"
            echo "$fline s/^.*"'$'"/ ${type} ${PARNAMELIST[$p]}                   ${PARVALUELIST[$p]}/" >> $sedScript
        fi
    else
        (( OUTPUT == 1 )) && echo "${PARNAMELIST[$p]} will be set to ${PARVALUELIST[$p]} on line $matchingLines"
        # change the leading star to a whitespace on the matching lines
        # if there was no leading star, this has no effect.
        for ml  in  $matchingLines; do
            echo "$ml s/^\*/ /" >> $sedScript
        done
    fi

    (( p++ ))
    done

    # the sed script is ready; we now use it to cast the SIF file
    # note that the cast problem and the sed script have similar names
    if (( failed == 0 )); then
        $SED -f $sedScript $problemLoc > $sifFile
    else
        $LN -s $problemLoc $sifFile
    fi

fi

# If necessary, create a symbolic link between the current directory
# and the problem file
# Since the SIF decoder does not want file names longer than 10 chars,
# take the last 10 numbers of the process id, in an attempt to have a
# unique name.

if [[  "$probDirNotGiven" == "true" && $nbparam == 0  ]]; then
    TEMPNAME=$PROBNAME
    if ((  lookInMastSif == 0  )); then
        link="false"
    else
        link="true"
        $LN -s ${problemLoc} ./$TEMPNAME.SIF
    fi
elif ((  nbparam > 0  )); then
    link="true"
    TEMPNAME=`echo $$ | $AWK 'BEGIN{nb=10}{l=length($1); if(l<nb){start=1} else{start=l-nb+1} print substr($1,start,nb)}'`
#    TEMPNAME="${PROBNAME}__${TEMPNAME}"
    $RM $TEMPNAME.SIF
    $LN -s $sifFile ./$TEMPNAME.SIF
else
    link="true"
    TEMPNAME=`echo $$ | $AWK 'BEGIN{nb=10}{l=length($1); if(l<nb){start=1} else{start=l-nb+1} print substr($1,start,nb)}'`
    #    TEMPNAME="${PROBNAME}__${TEMPNAME}"
    $RM $TEMPNAME.SIF
    $LN -s $PROBDIR/${PROBNAME}.SIF ./$TEMPNAME.SIF
fi

# Define the path to the decoder

DECODER=$MYSIFDEC/double/bin/sifdec
if [[ ! -x $DECODER  ]]; then
  DECODER=$MYSIFDEC/single/bin/sifdec
  if [[ ! -x $DECODER  ]]; then
    echo ' '
    echo "No SIF decoder sifdec found in $MYSIFDEC/double/bin"
    echo "or $MYSIFDEC/single/bin"
    echo 'Terminating execution.'
    exit 4
  fi
fi

if (( OUTPUT == 1 )); then
  echo 'convert the sif file into data and routines suitable for optimizer...'
  echo ' '
  echo 'problem details will be given'
  echo ' '
fi

[[ -e EXTERN.f ]] && $RM EXTERN.f

# construct input file for the decoder

#sdinput=${thiscaller}.$$.input
sdinput='SIFDECODE.CNF'

echo $TEMPNAME   > $sdinput
echo $METHOD    >> $sdinput
echo $OUTPUT    >> $sdinput
echo $PROBNAME  >> $sdinput
echo $automatic >> $sdinput
echo $ad0       >> $sdinput
echo $prec      >> $sdinput

# Finally, decode the problem

#$DECODER < $sdinput
$DECODER

# Clean up

$RM $sdinput
if ((  nbparam > 0  )); then
    $RM $sedScript
    $RM $sifFile
fi


if [[  $link == "true"  ]]; then
    $RM $TEMPNAME.SIF
fi

if [[ ! -e OUTSDIF.d ]]; then
  echo ' '
  echo "error exit from decoding stage. terminating execution."
  exit 3
fi

#  Rename files

if (( automatic != 0  )); then
  [[ -e ELFUND.f  ]] && $MV ELFUND.f ELFUN.f
  [[ -e GROUPD.f  ]] && $MV GROUPD.f GROUP.f
  [[ -e EXTERA.f  ]] && $MV EXTERA.f EXTER.f
else
  [[ -e ELFUNS.f  ]] && $MV ELFUNS.f ELFUN.f
  [[ -e GROUPS.f  ]] && $MV GROUPS.f GROUP.f
  [[ -e EXTERN.f  ]] && $MV EXTERN.f EXTER.f
fi

[[  -e RANGES.f  ]] && $MV RANGES.f RANGE.f

if [[ -e EXTER.f  ]]; then
 [[ -z EXTER.f  ]] && $RM EXTER.f
fi

#  Record the type of derivatives used in the decoding

echo $automatic $ad0 > AUTOMAT.d

(( OUTPUT == 1 )) && echo ' '
exit 0
