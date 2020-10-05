#!/bin/sh

# dgal: generic script to run GALAHAD package on architecture on the
#       data input in a problem-data file

#  * version for Bourne/bash shell

# Use: dgal architecture package [-s] [-h] [--help] [-k] [-p] [-o j] 
#                                [-l secs] probname[.data]
#
# where: options -s : run the single precision version
#                     (Default: run the double precision version)
#                -h : print this help and stop execution
#                -k : keep the executable after use
#                     (Default: remove the executable)
#                -p : profile the code (when possible)
#                -o : 0 for silent mode, 1 for brief description of
#                     the stages executed.
#                     (Default: -o 0)
#                -l : limit the cputime used to secs seconds
#                     (Default: 99999999 seconds)
#
#       probname      probname[.data] is the name of the file containing
#                     the data for the problem of interest.

#  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
#  Principal authors: Nick Gould, Dominique Orban and Philippe Toint

#  History -
#   originally released with GALAHAD Version 2.0. January 19th 2006

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
    echo ' Use: dgal architecture package [-s] [-h] [--help][-p] [-o j]'
    echo '                                [-l secs] probname[.data]'
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
#LIMIT=1800

# PRECISION = 0 (single precision), = 1 (double precision)

PRECISION=1

# KEEP = 0 (do not keep the executable), = 1 (keep it)

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
    echo ' Use: dgal architecture package [-s] [-h] [--help] [-p] [-o j]'
    echo '                                [-l secs] probname[.data]'
    echo ' '
    echo ' where: options -s : run the single precision version'
    echo '                     (Default: run the double precision version)'
    echo '                -h : print this help and stop execution'
    echo '                -p : profile the code (when possible)'
    echo '                -o : 0 for silent mode, 1 for brief description of'
    echo '                     the stages executed'
    echo '                     (Default: -o 0)'
    echo '                -l : limits the cputime to secs seconds'
    echo '                     (Default: unlimited cputime)'
    echo ' '
    echo '       probname      probname.data is the name of the file containing'
    echo '                     the data for the problem of interest.'
    exit 0
  elif [[ "$opt" == '-s' ]]; then
    PRECISION=0
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
  fi
  (( i++ ))
done

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

if [[ $PROFILE == "1" ]]; then
 pro="-p"
else
 pro=""
fi

probname=${!#}
if [[ ! -e $probname ]]; then
  if [[ -e $probname.data ]]; then
    probname=$probname.data
  else
    echo " No problem-data file $probname or $probname.data in current"
    echo " directory $cwd"
    exit 3
  fi
fi

#  If there are compiled, library versions of the blas
#  (basic linear algebra subprograms), set BLAS to a list of
#  names of the object library suffix -lx, where the object library
#  libx.a contains relevant blas. For instance if the blas are
#  shared between object libraries libblas1.a and libblas2.a,
#  BLAS should be set to "-lblas1 -lblas2", noting that those in
#  libblas1.a will take precedence over those in libblas2.a.
#  If compiled blas are unavailable, BLAS should be set to ""

BLAS=""

#  If there are compiled, library versions of the LAPACK library
#  set LAPACK to a list of names of the object library suffix -lx, 
#  where the object library libx.a contains relevant lapack. For instance 
#  if LAPACK is shared between object libraries liblapack1.a and liblapack2.a,
#  LAPACK should be set to "-llapack1 -llapack2", noting that those in
#  liblapack1.a will take precedence over those in liblapack2.a.
#  If compiled lapack are unavailable, LAPACK should be set to ""

LAPACK=""

#  If there is a compiled, library version of the Harwell
#  Subroutine Library, set HSL to -lx, where the object library
#  libx.a contains the relevant Harwell Subroutine Library.
#  For instance if the Harwell Subroutine Library is contained
#  in the object library libhsl.a, HSL should be set to "-lhsl".
#  If a compiled version of the Harwell Subroutine Library is
#  unavailable, HSL should be set to ""

HSL=""

#  If there is a compiled, library version of the Metis graph partitioning 
#  package (http://www-users.cs.umn.edu/~karypis/metis/) , set METIS to -lx, 
#  where the object library libx.a contains Metis.  For instance if Metis 
#  is contained in the object library libmetis.a, HSL should be set to 
#  "-lmetis".  If a compiled version of Metis is unavailable, METIS should 
#  be set to "". 
#  N.B. Metis is only required if MA57 (version 2 or later) is to be used.

METIS=""

set +C

#  directory for the main executable file

EXEC=$PWD

#  name of executable module

galmin=$1.$2

#  minimizer object codes to link

if [[ $PRECISION == "0" ]]; then
   PRECIS=single
   DOUBLE="s"
else
   PRECIS=double
   DOUBLE="d"
fi

#  machine-dependent bits

eval "`cat $GALAHAD/bin/sys/$1`"

#  directory for object files

GALOBJ=$GALAHAD/objects/$1/$PRECIS

#  libraries for BLAS, LAPACK, HSL and METIS

if [[ "$BLAS" == "" ]]; then
  BLASLIB="-lgalahad_blas"
else
  BLASLIB="$BLAS"
fi

if [[ "$LAPACK" == "" ]]; then
  LAPACKLIB="-lgalahad_lapack"
else
  LAPACKLIB="$LAPACK"
fi

if [[ "$HSL" == "" ]]; then
  HSLLIB="-lgalahad_hsl"
else
  HSLLIB="$HSL"
fi

if [[ "$METIS" == "" ]]; then
  METISLIB="-lgalahad_metis"
else
  METISLIB="$METIS"
fi

# link all the object files together.

if [[ $OUTPUT == "1" ]]; then
  echo ' '
  echo 'linking all the object files together ... '
  echo ' '
fi

#  ensure that package-dependent specification files are present

if [[ $2 == 'qpa' ]]; then
  RUNMAIN=$GALOBJ/inqpa.o
  [[ ! -e RUNQPA.SPC ]] && ln -s $GALAHAD/src/qpa/RUNQPA.SPC RUNQPA.SPC
elif [[ $2 == 'qpb' ]]; then
  RUNMAIN=$GALOBJ/inqpb.o
  [[ ! -e RUNQPB.SPC ]] && ln -s $GALAHAD/src/qpb/RUNQPB.SPC RUNQPB.SPC
elif [[ $2 == 'qpc' ]]; then
  RUNMAIN=$GALOBJ/inqpc.o
  [[ ! -e RUNQPC.SPC ]] && ln -s $GALAHAD/src/qpc/RUNQPC.SPC RUNQPC.SPC
elif [[ $2 == 'lpb' ]]; then
  RUNMAIN=$GALOBJ/inrunlpb.o
  [[ ! -e RUNLPB.SPC ]] && ln -s $GALAHAD/src/lpb/RUNLPB.SPC RUNLPB.SPC
elif [[ $2 == 'eqp' ]]; then
  RUNMAIN=$GALOBJ/ineqp.o
  [[ ! -e RUNEQP.SPC ]] && ln -s $GALAHAD/src/eqp/RUNEQP.SPC RUNEQP.SPC
#elif [[ $2 == 'wcp' ]]; then
#  RUNMAIN=$GALOBJ/inwcp.o
#  [[ ! -e RUNWCP.SPC ]] && ln -s $GALAHAD/src/wcp/RUNWCP.SPC RUNWCP.SPC
#elif [[ $2 == 'pre' ]]; then
#  RUNMAIN=$GALOBJ/inpre.o
#  [[ ! -e RUNPRE.SPC ]] && ln -s $GALAHAD/src/pre/RUNPRE.SPC RUNPRE.SPC
#elif [[ $2 == 'filt' ]]; then
#  RUNMAIN=$GALOBJ/infiltrane.o
#  [[ ! -e RUNFILT.SPC ]] && ln -s $GALAHAD/src/filtrane/RUNFILT.SPC RUNFILT.SPC
#elif [[ $2 == 'superb' ]]; then
#  RUNMAIN=$GALOBJ/insuperb.o
#  [[ ! -e RUNSUPERB.SPC ]] && ln -s $GALAHAD/src/superb/RUNSUPERB.SPC RUNSUPERB.SPC
#elif [[ $2 == 'fastr' ]]; then
#  RUNMAIN=$GALOBJ/infastr.o
#  [[ ! -e RUNFASTR.SPC ]] && ln -s $GALAHAD/src/fastr/RUNFASTR.SPC RUNFASTR.SPC
#elif [[ $2 == 'trtn' ]]; then
#  RUNMAIN=$GALOBJ/intrtn.o
#  [[ ! -e RUNTRTN.SPC ]] && ln -s $GALAHAD/src/trtn/RUNTRTN.SPC RUNTRTN.SPC
#elif [[ $2 == 'lpsqpa' ]]; then
#  RUNMAIN=$GALOBJ/inlpsqpa.o
#  [[ ! -e RUNLPSQP.SPC ]] && ln -s $GALAHAD/src/superb/RUNLPSQP.SPC RUNLPSQP.SPC
#elif [[ $2 == 'lpsqp' ]]; then
#  RUNMAIN=$GALOBJ/inlpsqp.o
#  [[ ! -e RUNLPSQP.SPC ]] && ln -s $GALAHAD/src/superb/RUNLPSQP.SPC RUNLPSQP.SPC
#else
## Default to LANCELOT B if necessary - 
## CUTEr and LAPACK are not needed in this case
#  RUNMAIN=$GALOBJ/inlanb.o
#  [[ ! -e RUNLANB.SPC ]] && ln -s $GALAHAD/src/lanb/RUNLANB.SPC RUNLANB.SPC
#  CUTERLIB=""
#  LAPACKLIB=""
else
  echo " Unfortunately no data file version of $2 is available at presnt."
  exit 2
fi

#  create the executable

$FORTRAN $FFLAGS -o $galmin $RUNMAIN -L$GALOBJ -lgalahad $CUTERLIB -lgalahad \
         $HSLLIB $METISLIB $LAPACKLIB $BLASLIB 

[[ $PWD != $EXEC ]] && $MV $galmin $EXEC/$galmin

#  run $galmin on the current test problem.

if [[ $OUTPUT == "1" ]]; then
  echo ' '
  echo "running $2 on current test problem ... "
  echo ' '
fi

ulimit -t $LIMIT
if [[ $PROFILE == "1" ]]; then
  ISPIXIE=`whereis -b pixie`
  if [[ ${#ISPIXIE[@]} == 2 ]]; then
#   atom $EXEC/$galmin -tool pixie -w0 -toolargs="-quiet" >  2>&1/dev/null
    pixie -quiet $EXEC/$galmin > /dev/null 2>&1
    $EXEC/$galmin.pixie < $probname
    prof -pixie -lines $EXEC/$galmin > $EXEC/$galmin.pixie.out
    $RM $EXEC/$galmin.pixie $EXEC/$galmin.Counts $EXEC/$galmin.Addrs
  else
    if [[ $OUTPUT == "1" ]]; then
      echo 'no profiling available, sorry ... '
      echo ' '
    fi
    $EXEC/$galmin < $probname
  fi
else
  $EXEC/$galmin < $probname
fi

#  tidy up the current directory, deleting all junk if required

if [[ $KEEP == "0" ]]; then
  $RM $EXEC/$galmin
fi
