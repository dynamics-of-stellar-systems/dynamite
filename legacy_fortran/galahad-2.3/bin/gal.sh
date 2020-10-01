#!/bin/bash

# gal: generic script to apply package on architecture and delete the 
#      executable after use.
#  * version for Bourne/bash shell

# Use: gal architecture package [-e] [-h] [-r] [-s] [-k] [-p] [-o i] [-l secs]

# where: options -e : use the load module architecture.package if it exists
#                     (Default: create a new load module)
#                -h : print this help and stop execution
#                -r : do not recompile the problem functions
#                -s : run the single precision version
#                     (Default: run the double precision version)
#                -k : keep the load module after use
#                     (Default: delete the load module)
#                -p : profile the code (when possible)
#                -o : 0 for silent mode, 1 for brief description of
#                     the stages executed
#                     (Default: -o 0)
#                -l : limit the cputime used to secs seconds
#                     (Default: -l 99999999)

#  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
#  Principal authors: Nick Gould, Dominique Orban and Philippe Toint

#  History -
#   originally released pre GALAHAD Version 1.0. August 20th, 1999 (for csh)
#   update released with GALAHAD Version 2.0. May 11th 2006 (for sh)

if [[ -z $GALAHAD ]]; then
  echo ' You have not set the GALAHAD environment variable. '
  echo ' This needs to point to the main directory in which '
  echo ' you have installed the GALAHAD package. '
  exit 5
fi

let last=$#

stop_early="false"
if (( last < 2 )); then
    stop_early="true"
elif [[ "`echo $1 | grep -e '^-'`" != "" || "`echo $2 | grep -e '^-'`" != ""  ]]; then
    stop_early="true"
fi

if [[ "$stop_early" == "true" ]]; then
    echo " Use: gal architecture package [-e] [-h] [-r] [-s] [-k] [-p] [-o i] [-l secs]"
    exit 1
fi

set +C

#  directory for the main executable file

EXEC=$PWD

#  directory for temporary files

TMP=/tmp

#  variables for each option

# PRECISION = 0 (single precision), = 1 (double precision)

PRECISION=1

# RECOMPILE = 0 (use existing problem objects), = 1 (recompile problem objects)

RECOMPILE=0

#  AUTOMATIC = 0 (provided), = 1 (automatic forward), = 2 (automatic backward)
#  AD = 0 (none), =1 (AD01 used), =2 (AD02 used)

if [[ -e AUTOMAT.d ]]; then
  DUM=( `cat AUTOMAT.d` )
  AUTOMATIC=${DUM[1]}
  AD0=${DUM[2]}
else
  AUTOMATIC=0
  AD0=0
fi

# NEW = 0 (run existing f module), = 1 (build a new module)

NEW=1

# KEEP = 0 (discard f load module after use), = 1 (keep it)

KEEP=0

# PROFILE = 0 (do not profile the code), = 1 (profile it)

PROFILE=0

# OUTPUT = 0 (summary output), = 1 (detailed output from decoder)

OUTPUT=0

# LIMIT = 0 (no cputime limit)

LIMIT=99999999
#LIMIT=1800

#  name of executable module

galmin=$1.$2

#  interpret arguments

let i=3

while (( i <= last ))
do
  opt=${!i}
  if [[ "$opt" == '-e' ]]; then
    NEW=0
  elif [[ "$opt" == '-r' ]]; then
    RECOMPILE=1
  elif [[ "$opt" == '-s' ]]; then
    PRECISION=0
  elif [[ "$opt" == '-h' || "$opt" == '--help' ]]; then
    echo " Use: gal architecture package [-e] [-h] [-r] [-s] [-k] [-p] [-o i] [-l secs]"
    echo ' '
    echo " where: options -e : use the load module $galmin if it exists"
    echo '                     (Default: create a new load module)'
    echo '                -h : print this help and stop execution'
    echo '                -r : recompile the problem functions'
    echo '                -s : run the single precision version'
    echo '                     (Default: run the double precision version)'
    echo '                -k : keep the load module after use '
    echo '                     (Default: delete the load module)'
    echo '                -p : profile the code (when possible)'
    echo '                -o : 0 for silent mode, 1 for brief description of'
    echo '                     the stages executed'
    echo '                     (Default: -o 0)'
    echo '                -l : limits the cputime to secs seconds'
    echo '                     (Default: -l 99999999)'
    exit 0
  elif [[ "$opt" == '-k' ]]; then
    KEEP=1
  elif [[ "$opt" == '-p' ]]; then
    PROFILE=1
  elif [[ "$opt" == '-o' ]]; then
    (( i++ ))
    OUTPUT=${!i}
  elif [[ "$opt" == '-l' ]]; then
    (( i++ ))
    LIMIT=${!i}
  else
    echo " Use: gal architecture package [-e] [-h] [-r] [-s] [-k] [-p] [-o i] [-l secs]"
    exit 1
  fi
  (( i++ ))
done

#  minimizer object codes to link

if [[ $PRECISION == "0" ]]; then
   PRECIS=single
   DOUBLE="s"
else
   PRECIS=double
   DOUBLE="d"
fi

#  ----------------------------------------------------------------------------
#  -*- Default values that will be overridden if set in $GALAHAD/bin/sys/$1 -*-
#  ----------------------------------------------------------------------------

#  standard unix commands (rm, make, cat, sed, mv, ls)

RM="rm -f"
MAKE="make"
CAT="cat"
SED="sed"
MV="mv"
LS="ls"

#  the command that invokes the fortran 95 compiler

FORTRAN="f95"

#  compiler flags for linking

FFLAGS=""

#  flags for compiling the fortran 77 problem-dependent roiutines 

PROBFLAGS="-c -fixed"

#  the location of the version of CUTEr to be used

CUTERUSED="$MYCUTER"

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

#  ----------------------------------------------------------------------------
#  -*- end of default values that may be overridden in $GALAHAD/bin/sys/$1 -*-
#  ----------------------------------------------------------------------------

#  machine-dependent bits

#eval "`cat $GALAHAD/bin/sys/$1`"
. ${GALAHAD}/bin/sys/$1

#  run galmin without rebuilding it

if [[ $NEW == "0" ]]; then
  if [[ ! -e $EXEC/$galmin || ! -x $EXEC/$galmin ]]; then
    echo ' '
    echo 'load module gal not found/executable. Rerun with -e option'
    echo ' '
    exit 3
  fi
  if [[ $OUTPUT ]]; then
    echo ' '
    echo "running $2 on current test problem ... "
    echo ' '
  fi
  #limit cputime $LIMIT
  ulimit -t $LIMIT
  if [[ $PROFILE == "1" ]]; then
    ISPIXIE=`whereis -b pixie`
    if [[ ${#ISPIXIE[@]} == 2 ]]; then
      pixie -quiet $EXEC/$galmin > /dev/null 2>&1
      $EXEC/$galmin.pixie
      prof -pixie -lines $EXEC/$galmin > $EXEC/$galmin.pixie.out
      $RM $EXEC/$galmin.pixie $EXEC/$galmin.Counts $EXEC/$galmin.Addrs
    else
      if [[ $OUTPUT == "1" ]]; then
        echo 'no profiling available, sorry ... '
        echo ' '
      fi
      $EXEC/$galmin
    fi
  else
    $EXEC/$galmin
  fi

#  tidy up the current directory, deleting all junk.

  [[ $KEEP == "0" ]] && $RM $EXEC/$galmin
  exit 0
fi

#  check that CUTEr has been installed in the location indicated

if [[ $2 != 'lanb' ]]; then
  if [[ ! -e $CUTERUSED ]]; then
    echo ' The CUTERUSED environment variable does not point.to'
    echo ' an installed version of CUTER. Set the correct location'
    echo ' in  $GALAHAD/bin/sys/$1'
    echo " (CUTERUSED = $CUTERUSED)"
    exit 6
  fi
fi

#  build $galmin and tools

#  directory for object files

GALOBJ=$GALAHAD/objects/$1/$PRECIS

if [[ $PRECISION == "0" ]]; then
   CUTERLIB="-L$CUTERUSED/single/lib -lcuter"
else
   CUTERLIB="-L$CUTERUSED/double/lib -lcuter"
fi
PROBLIB=""

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

# ensure that the current test problem has been compiled.

if [[ $OUTPUT == "1" ]]; then
  echo 'compiling the current test problem, if that is necessary ... '
  echo ' '
fi

[[ -e RANGE.o && $RECOMPILE == '0' ]] && $RM RANGE.o
[[ -e ELFUN.o && $RECOMPILE == '0' ]] && $RM ELFUN.o
[[ -e GROUP.o && $RECOMPILE == '0' ]] && $RM GROUP.o
[[ -e EXTER.o && $RECOMPILE == '0' ]] && $RM EXTER.o

NSUB=( "ELFUN.o GROUP.o RANGE.o" )
[[ -e EXTER.f && ! -z EXTER.f ]] && NSUB=( "$NSUB EXTER.o" )

for i  in  $NSUB; do
  if [[ ! -e $i ]]; then
    j=`basename $i .o`
    cp ${j}.f ${j}.f90
    $FORTRAN $PROBFLAGS ${j}.f90
    if [[ $? != 0 ]]; then
      exit 1
    fi
    $RM ${j}.f90
  fi
done

# link all the tools files together.

if [[ $OUTPUT == "1" ]]; then
  echo ' '
  echo 'linking all the object files together ... '
  echo ' '
fi

#  ensure that package-dependent specification files are present

if [[ $2 == 'qpa' ]]; then
  RUNMAIN=$GALOBJ/runqpa_sif.o
  [[ ! -e RUNQPA.SPC ]] && ln -s $GALAHAD/src/qpa/RUNQPA.SPC RUNQPA.SPC
elif [[ $2 == 'qpb' ]]; then
  RUNMAIN=$GALOBJ/runqpb_sif.o
  [[ ! -e RUNQPB.SPC ]] && ln -s $GALAHAD/src/qpb/RUNQPB.SPC RUNQPB.SPC
elif [[ $2 == 'qpc' ]]; then
  RUNMAIN=$GALOBJ/runqpc_sif.o
  [[ ! -e RUNQPC.SPC ]] && ln -s $GALAHAD/src/qpc/RUNQPC.SPC RUNQPC.SPC
elif [[ $2 == 'lpb' ]]; then
  RUNMAIN=$GALOBJ/runlpb_sif.o
  [[ ! -e RUNLPB.SPC ]] && ln -s $GALAHAD/src/lpb/RUNLPB.SPC RUNLPB.SPC
elif [[ $2 == 'wcp' ]]; then
  RUNMAIN=$GALOBJ/runwcp_sif.o
  [[ ! -e RUNWCP.SPC ]] && ln -s $GALAHAD/src/wcp/RUNWCP.SPC RUNWCP.SPC
elif [[ $2 == 'lcf' ]]; then
  RUNMAIN=$GALOBJ/runlcf_sif.o
  [[ ! -e RUNLCF.SPC ]] && ln -s $GALAHAD/src/lcf/RUNLCF.SPC RUNLCF.SPC
elif [[ $2 == 'pre' ]]; then
  RUNMAIN=$GALOBJ/runpre_sif.o
  [[ ! -e RUNPRE.SPC ]] && ln -s $GALAHAD/src/pre/RUNPRE.SPC RUNPRE.SPC
elif [[ $2 == 'eqp' ]]; then
  RUNMAIN=$GALOBJ/runeqp_sif.o
  [[ ! -e RUNEQP.SPC ]] && ln -s $GALAHAD/src/eqp/RUNEQP.SPC RUNEQP.SPC
elif [[ $2 == 'lls' ]]; then
  RUNMAIN=$GALOBJ/runlls_sif.o
  [[ ! -e RUNLLS.SPC ]] && ln -s $GALAHAD/src/lls/RUNLLS.SPC RUNLLS.SPC
elif [[ $2 == 'filt' ]]; then
  RUNMAIN=$GALOBJ/runfiltrane_sif.o
  [[ ! -e RUNFILT.SPC ]] && ln -s $GALAHAD/src/filtrane/RUNFILT.SPC RUNFILT.SPC
elif [[ $2 == 'superb' ]]; then
  RUNMAIN=$GALOBJ/runsuperb_sif.o
  [[ ! -e RUNSUPERB.SPC ]] && ln -s $GALAHAD/src/superb/RUNSUPERB.SPC RUNSUPERB.SPC
elif [[ $2 == 'tru' ]]; then
  RUNMAIN=$GALOBJ/runtru_sif.o
  PROBLIB="-lgalahad_problem"
  [[ ! -e RUNTRU.SPC ]] && ln -s $GALAHAD/src/tru/RUNTRU.SPC RUNTRU.SPC
elif [[ $2 == 'aco' ]]; then
  RUNMAIN=$GALOBJ/runaco_sif.o
  PROBLIB="-lgalahad_problem"
  [[ ! -e RUNACO.SPC ]] && ln -s $GALAHAD/src/aco/RUNACO.SPC RUNACO.SPC
elif [[ $2 == 'acoc' ]]; then
  RUNMAIN=$GALOBJ/runacoc_sif.o
  PROBLIB="-lgalahad_problem"
  [[ ! -e RUNACOC.SPC ]] && ln -s $GALAHAD/src/acoc/RUNACOC.SPC RUNACOC.SPC
elif [[ $2 == 'fastr' ]]; then
  RUNMAIN=$GALOBJ/runfastr_sif.o
  [[ ! -e RUNFASTR.SPC ]] && ln -s $GALAHAD/src/fastr/RUNFASTR.SPC RUNFASTR.SPC
elif [[ $2 == 'funnel' ]]; then
  RUNMAIN=$GALOBJ/runfunnel_sif.o
  PROBLIB="-lgalahad_problem"
  [[ ! -e RUNFUNNEL.SPC ]] && ln -s $GALAHAD/src/funnel/RUNFUNNEL.SPC RUNFUNNEL.SPC
elif [[ $2 == 'sqp' ]]; then
  RUNMAIN=$GALOBJ/runsqp_sif.o
  [[ ! -e RUNSQP.SPC ]] && ln -s $GALAHAD/src/sqp/RUNSQP.SPC RUNSQP.SPC
elif [[ $2 == 'trimsqp' ]]; then
  RUNMAIN=$GALOBJ/runtrimsqp_sif.o
  PROBLIB="-lgalahad_problem"
  [[ ! -e RUNTRIMSQP.SPC ]] && ln -s $GALAHAD/src/trimsqp/RUNTRIMSQP.SPC RUNTRIMSQP.SPC
elif [[ $2 == 'trtn' ]]; then
  RUNMAIN=$GALOBJ/runtrtn_sif.o
  [[ ! -e RUNTRTN.SPC ]] && ln -s $GALAHAD/src/trtn/RUNTRTN.SPC RUNTRTN.SPC
elif [[ $2 == 'lpsqpa' ]]; then
  RUNMAIN=$GALOBJ/runlpsqpa.o
  [[ ! -e RUNLPSQP.SPC ]] && ln -s $GALAHAD/src/superb/RUNLPSQP.SPC RUNLPSQP.SPC
elif [[ $2 == 'lpsqp' ]]; then
  RUNMAIN=$GALOBJ/runlpsqp.o
  [[ ! -e RUNLPSQP.SPC ]] && ln -s $GALAHAD/src/superb/RUNLPSQP.SPC RUNLPSQP.SPC
else
# Default to LANCELOT B if necessary - 
# CUTEr and LAPACK are not needed in this case
  RUNMAIN=$GALOBJ/runlanb_sif.o
  [[ ! -e RUNLANB.SPC ]] && ln -s $GALAHAD/src/lanb/RUNLANB.SPC RUNLANB.SPC
  CUTERLIB=""
  LAPACKLIB=""
fi

# See if there already is a decoded problem in the current directory
# and make sure it is suitable for the required package

if [ -e OUTSDIF.d ]; then
  m=`head -2 OUTSDIF.d | tail -1 | ${SED} -e 's/^[ ]*//' | cut -c 1`
  if [[ $2 == 'lanb' ]]; then
    if [ "$m" == "3" ]; then
      echo 'The decoded files in the current directory are not suitable'
      echo 'for input to LANCELOT-B. Please re-run with sdgal'
      exit 10
    fi
  else
    if [ "$m" != "3" ]; then
      echo 'The decoded files in the current directory are only suitable'
      echo 'for input to LANCELOT-B. Please re-run with sdgal'
      exit 10
    fi
  fi
else
  echo 'There does not appear to be a decoded problem in the current directory'
  echo 'Please re-run with sdgal'
fi


#  create the executable

$FORTRAN $FFLAGS -o $galmin $RUNMAIN $NSUB \
         -L$GALOBJ -lgalahad $PROBLIB $CUTERLIB -lgalahad \
         $HSLLIB $METISLIB $LAPACKLIB $BLASLIB 

[[ $PWD != $EXEC ]] && $MV $galmin $EXEC/$galmin

#  run $galmin on the current test problem.

if [[ $OUTPUT == "1" ]]; then
  echo ' '
  echo "running $2 on current test problem ... "
  echo ' '
fi

#limit cputime $LIMIT
ulimit -t $LIMIT
if [[ $PROFILE == "1" ]]; then
  ISPIXIE=`whereis -b pixie`
  if [[ ${#ISPIXIE[@]} == 2 ]]; then
#   atom $EXEC/$galmin -tool pixie -w0 -toolargs="-quiet" >  2>&1/dev/null
    pixie -quiet $EXEC/$galmin > /dev/null 2>&1
    $EXEC/$galmin.pixie
    prof -pixie -lines $EXEC/$galmin > $EXEC/$galmin.pixie.out
    $RM $EXEC/$galmin.pixie $EXEC/$galmin.Counts $EXEC/$galmin.Addrs
  else
    if [[ $OUTPUT == "1" ]]; then
      echo 'no profiling available, sorry ... '
      echo ' '
    fi
    $EXEC/$galmin
  fi
else
  $EXEC/$galmin
fi

#  tidy up the current directory, deleting all junk.

[[ $KEEP == "0" ]] && $RM $EXEC/$galmin

