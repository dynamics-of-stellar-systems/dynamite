#!/bin/bash

# sdgal: generic script to decode an SIF file and then run package
#        on architecture on the output
#  * version for Bourne/bash shell

# Use: sdgal architecture package [-s] [-h] [-k] [-f] [-b] [-a j]
#                            [-p] [-o j] [-l secs] probname[.SIF]
#
# where: options -s : run the single precision version
#                     (Default: run the double precision version)
#                -h : print this help and stop execution
#                -k : keep the load module after use
#                     (Default: delete the load module)
#                -f : use automatic differentiation in forward mode
#                -b : use automatic differentiation in backward mode
#                -a : 1 use the older HSL automatic differentiation package AD01
#                     2 use the newer HSL automatic differentiation package AD02
#                     (Default: -a 2)
#                -p : profile the code (when possible)
#                -o : 0 for silent mode, 1 for brief description of
#                     the stages executed.
#                     (Default: -o 0)
#                -l : limit the cputime used to secs seconds
#                     (Default: 99999999 seconds)
#
#       probname      probname[.SIF] is the name of the file containing
#                     the SIF file for the problem of interest.

#  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
#  Principal authors: Nick Gould, Dominique Orban and Philippe Toint

#  History -
#   originally released pre GALAHAD Version 1.0. August 20th, 1999 (for csh)
#   update released with GALAHAD Version 2.0. August 11th 2005 (for sh)

if [[ -z $GALAHAD ]]; then
  echo ' You have not the=true GALAHAD environment variable. '
  echo ' This needs to point to the main directory in which '
  echo ' you have installed the GALAHAD package. '
  exit 5
fi

# Obtain number of arguments
let last=$#
#(( last=last-1 ))

stop_early="false"
if (( last < 2 )); then
    stop_early="true"
elif [[ "`echo $1 | grep -e '^-'`" != "" || "`echo $2 | grep -e '^-'`" != ""  ]]; then
    stop_early="true"
fi

if [[ "$stop_early" == "true" ]]; then
    echo ' Use: sdgal architecture package [-s] [-h] [-k] [-f] [-b] [-a j]'
    echo '                            [-p] [-o j] [-l secs] probname[.SIF]'
    exit 1
fi

if [[ $2 == 'lanb'  ]]; then
  METHARG="-m 0"
else
  METHARG=""
fi

#  directory for temporary files

TMP=/tmp

#  variables for each option

#  LIMIT (maximum cputime for running object file)

LIMIT=99999999
LIMIT=1800

# PRECISION = 0 (single precision), = 1 (double precision)

PRECISION=1

#   AUTOMATIC = 0 (provided), = 1 (automatic forward), = 2 (automatic backward)

AUTOMATIC=0

#   AD0 = 1 (AD01 used), = 2 (AD02 used)

AD0=2

# KEEP = 0 (discard load module after use), = 1 (keep it)

KEEP=0

# PROFILE = 0 (do not profile the code), = 1 (profile it)

PROFILE=0

# OUTPUT = 0 (summary output), = 1 (detailed output from decoder)

OUTPUT=0

let i=3

while (( i <= last ))
do
  opt=${!i}
  echo " opt = $opt"
  if [[ "$opt" == '-h' || "$opt" == '--help' ]]; then
    echo ' Use: sdgal architecture package [-s] [-h|--help] [-k] [-f] [-b] [-a j]'
    echo '                            [-p] [-o j] [-l secs] probname[.SIF]'
    echo ' '
    echo ' where: options -s : run the single precision version'
    echo '                     (Default: run the double precision version)'
    echo '                -h : print this help and stop execution'
    echo '                -k : keep the load module after use '
    echo '                     (Default: delete the load module)'
    echo '                -f : use automatic differentiation in forward mode'
    echo '                -b : use automatic differentiation in backward mode'
    echo '                -a : 1 use the older HSL automatic differentiation'
    echo '                     package AD01, 2 use the newer HSL automatic '
    echo '                     differentiation package AD02 (Default: -a 2)'
    echo '                -p : profile the code (when possible)'
    echo '                -o : 0 for silent mode, 1 for brief description of'
    echo '                     the stages executed'
    echo '                     (Default: -o 0)'
    echo '                -l : limits the cputime to secs seconds'
    echo '                     (Default: unlimited cputime)'
    echo ' '
    echo '       probname      probname.SIF is the name of the file containing'
    echo '                     the SIF file for the problem of interest.'
    exit 0
  elif [[ "$opt" == '-s' ]]; then
    PRECISION=0
  elif [[ "$opt" == '-f' ]]; then
    AUTOMATIC=1
  elif [[ "$opt" == '-b' ]]; then
    AUTOMATIC=2
  elif [[ "$opt" == '-k' ]]; then
    KEEP=1
  elif [[ "$opt" == '-p' ]]; then
    PROFILE=1
  elif [[ "$opt" == '-o' ]]; then
    (( i++ ))
    OUTPUT=${!i}
  elif [[ "$opt" == '-m' ]]; then
    (( i++ ))
    METHARG=""
  elif [[ "$opt" == '-l' ]]; then
    (( i++ ))
    LIMIT=${!i}
  elif [[ "$opt" == '-a' ]]; then
    (( i++ ))
    AD0=${!i}
  fi
  (( i++ ))
done

# call with two argument allows user to choose minimization method

if [[ $PRECISION == "1" ]]; then
 p=""
 up="d"
 PRECIS=double
else
 p="-s"
 up="s"
 PRECIS=single
fi

#  machine-dependent bits

eval "`cat $GALAHAD/bin/sys/$1`"
#. $GALAHAD/bin/sys/$1

if [[ $KEEP == "1" ]]; then
 k="-k"
else
 k=""
fi

if [[ $PROFILE == "1" ]]; then
 pro="-p"
else
 pro=""
fi

probname=${!#}
#probname=${@[$last]}

if [[ $OUTPUT == "1" ]]; then
  echo 'convert the sif file into data and routines suitable for optimizer...'
  echo ' '
  echo 'problem details will be given'
  echo ' '
fi

[[ -e EXTER.f ]] && $RM EXTER.f

# decode the problem

if [[ ! -z $MYSIFDEC  ]]; then
  $MYSIFDEC/bin/sifdecode $METHARG $*
  [[ $? != 0  ]] && exit $status
else
  echo " sdgal : environment variable MYSIFDEC not set"
  echo "         Either SifDec is not installed or you"
  echo "         need to properly set MYSIFDEC"
  exit 4
fi

#  Check for decoding errors

[[ $OUTPUT == "1" ]] && echo ' '
if [[ ! -e OUTSDIF.d ]]; then
  echo ' '
  echo "error exit from decoding stage. terminating execution."
  exit 3
fi
[[ $OUTPUT == "1" ]] && echo ' '

#  Record the type of derivatives used in the decoding

[[ -f AUTOMAT.d ]] && $RM AUTOMAT.d
echo $AUTOMATIC $AD0 > AUTOMAT.d

#  run the program on the output

$GALAHAD/bin/gal.sh $1 $2 $p $k $pro -o $OUTPUT -l $LIMIT

$RM $TMP/sdgal.input

