#!/bin/sh

# gampl: generic script to pass an Ampl model to a GALAHAD package
# version for Bourne/bash shell

# Use: gampl architecture [-h] [-l secs] [galahad_ampl options]

# where: architecture is the required GALAHAD architecture
#        options       -h : displays usage message
#                      -l : limit the cputime used to secs seconds
#                           (Default: -l 99999999)
#    galahad_ampl options : if specified, options passed directly
#                           to the main Galahad/Ampl driver.

#  Nick Gould, Dominique Orban and Philippe Toint,
#  March 2003, for GALAHAD Productions.
# ( Last modified on 26 July 2005 at 20:00:00 GMT )

if [[ -z $GALAHAD ]]; then
  echo ' You have not the=true GALAHAD environment variable. '
  echo ' This needs to point to the main directory in which '
  echo ' you have installed the GALAHAD package. '
  exit 5
fi

# LIMIT = 0 (no cputime limit)

LIMIT=99999999

# Galahad-Ampl-specific options, to be passed to the main driver

GAMPL_OPT=''

# Default Galahad package

default='qpb'

#  interpret arguments

let last=$#
let i=1

if (( last == 0 )); then
    echo " Use: gampl architecture [-h] [-l secs] [galahad_ampl-specific options]"
    echo ' '
    echo ' where: architecture is the *required* GALAHAD architecture'
    echo '        options -h : print this help and stop execution'
    echo '                -l : limits the cputime to secs seconds'
    echo '                     (Default: -l 99999999)'
    echo ' galahad_ampl-specific options, if specified, are passed'
    echo ' as-is to the main Galahad/Ampl driver. See  man gampl  for details.'
    exit 2
fi

while (( i <= last ))
do
  opt=${!i}
  if [[ "$i" == '1' ]]; then
    AMPLBIN=$GALAHAD/ampl_bin/$opt
  elif [[ "$opt" == '-h' ]]; then
    echo " Use: gampl architecture [-h] [-l secs] [galahad_ampl-specific options]"
    echo ' '
    echo ' where: architecture is the *required* GALAHAD architecture'
    echo '        options -h : print this help and stop execution'
    echo '                -l : limits the cputime to secs seconds'
    echo '                     (Default: -l 99999999)'
    echo ' galahad_ampl-specific options, if specified, are passed'
    echo ' as-is to the main Galahad/Ampl driver. See  man gampl  for details.'
    exit 0
  elif [[ "$opt" == '-v' ]]; then
    cat $GALAHAD/version
    exit 0
  elif [[ "$opt" == '-l' ]]; then
    (( i++ ))
    LIMIT=${!i}
  else
    [[  "$opt" == 'qpa=1'   ]] && default='qpa'
    [[  "$opt" == 'qpb=1'   ]] && default='qpb'
    [[  "$opt" == 'qpc=1'   ]] && default='qpc'
    [[  "$opt" == 'pre=1'   ]] && default='pre'
    [[  "$opt" == 'lanb=1'  ]] && default='lanb'
    [[  "$opt" == 'filt=1'  ]] && default='filt'
    GAMPL_OPT="$GAMPL_OPT $opt"
  fi
  (( i++ ))
done

#  Ensure that the AMPL binary is present

if [[ ! -e $AMPLBIN ]]; then
  echo " AMPL binary for architecture $opt does not exist."
  echo " Please ensure that the binary has been installed"
  echo " in $GALAHAD/ampl_bin"
  exit 1
fi

#  Ensure that package-dependent specification files are present

spcset=""
if [[ $default == 'qpa'  ]]; then
  [[ ! -e RUNQPA.SPC ]] && spcset=RUNQPA.SPC
  [[ ! -e RUNQPA.SPC ]] && ln -s $GALAHAD/src/qpa/RUNQPA.SPC RUNQPA.SPC
elif [[ $default == 'qpb'  ]]; then
  [[ ! -e RUNQPB.SPC ]] && spcset=RUNQPB.SPC
  [[ ! -e RUNQPB.SPC ]] && ln -s $GALAHAD/src/qpb/RUNQPB.SPC RUNQPB.SPC
elif [[ $default == 'qpc'  ]]; then
  [[ ! -e RUNQPC.SPC ]] && spcset=RUNQPC.SPC
  [[ ! -e RUNQPC.SPC ]] && ln -s $GALAHAD/src/qpc/RUNQPC.SPC RUNQPC.SPC
elif [[ $default == 'pre'  ]]; then
  [[ ! -e RUNPRE.SPC ]] && spcset=RUNPRE.SPC
  [[ ! -e RUNPRE.SPC ]] && ln -s $GALAHAD/src/pre/RUNPRE.SPC RUNPRE.SPC
elif [[ $default == 'lanb'  ]]; then
  [[ ! -e RUNLANB.SPC ]] && spcset=RUNLANB.SPC
  [[ ! -e RUNLANB.SPC ]] && ln -s $GALAHAD/src/lanb/RUNLANB.SPC RUNLANB.SPC
elif [[ $default == 'filt'  ]]; then
  [[ ! -e RUNFILT.SPC ]] && spcset=RUNFILT.SPC
  [[ ! -e RUNFILT.SPC ]] && ln -s $GALAHAD/src/filtrane/RUNFILT.SPC RUNFILT.SPC
else 
  [[ ! -e RUNQPB.SPC ]] && spcset=RUNQPB.SPC
  [[ ! -e RUNQPB.SPC ]] && ln -s $GALAHAD/src/qpb/RUNQPB.SPC RUNQPB.SPC
fi

# Limit runtime, if required

ulimit -t $LIMIT

# Execute package

#echo $GAMPL_OPT
$AMPLBIN $GAMPL_OPT

#  Terminate gracefully

[[ ! -z $spcset ]] && rm -f $spcset
if [[ -e fort.6 ]]; then
  cat fort.6
  rm -r fort.6
fi
