#!/bin/csh -f
#  ( Last modified on 23 Dec 2000 at 17:29:56 )
#
# tao: apply TAO to a problem and delete the executable after use.
#
#{version}
#

#{args}

#
#  N. Gould, D. Orban & Ph. Toint, November 7th, 2000
#

#{cmds}

#
# Environment check
#

envcheck
if( $status != 0 ) exit $status

#
# Check to make sure PETSc, TAO environment variables are ok, libraries exist
#

set MPI_LIBDIR = '/usr/local/lib'
set MPI_LIB = "-L$MPI_LIBDIR -lmpich"
if ( ! -e $MPI_LIBDIR/libmpich.a && ! -e $MPI_LIBDIR/libmpich.so) then
    echo Error: $MPI_LIBDIR/libmpich.a or $MPI_LIBDIR/libmpich.so not found.
    echo Either install mpich or alter tao script to describe mpi location
    exit 1
endif

if ( ! ( $?PETSC_DIR && $?BOPT && $?TAO_DIR && $?PETSC_ARCH ) )then
    echo "Error: You must set the environment variables PETSC_DIR, TAO_DIR,"
    echo "BOPT (=g_c++ or O_c++) and PETSC_ARCH (linux, solaris, etc.)"
    echo "PETSC_DIR = $PETSC_DIR"
    echo "PETSC_ARCH = $PETSC_ARCH"
    echo "TAO_DIR = $TAO_DIR"
    echo "BOPT = $BOPT\n"
    exit 1
endif

echo "PETSC_DIR = $PETSC_DIR"
echo "PETSC_ARCH = $PETSC_ARCH"
echo "TAO_DIR = $TAO_DIR"
echo "BOPT = $BOPT\n"

foreach i ( petsc petscfortran petscsnes petscsles petscdm petscmat petscvec )
    if ( ! -e $PETSC_DIR/lib/lib$BOPT/$PETSC_ARCH/lib$i.so && ! -e $PETSC_DIR/lib/lib$BOPT/$PETSC_ARCH/lib$i.a ) then
	echo "Error: petsc library -l$i not found in $PETSC_DIR/lib/lib$BOPT/$PETSC_ARCH"
	exit 1
    endif
end

foreach i ( tao taofortran )
    if ( ! -e $TAO_DIR/lib/lib$BOPT/$PETSC_ARCH/lib$i.so && ! -e $TAO_DIR/lib/lib$BOPT/$PETSC_ARCH/lib$i.a ) then
	echo "Error: tao library -l$i not found in $TAO_DIR/lib/lib$BOPT/$PETSC_ARCH"
	exit 1
    endif
end

#
#  define a short acronym for the package to which you wish to make an interface
#

setenv caller tao
setenv PAC tao

#
#  define the name of the subdirectory of $CUTERDIR/common/src/pkg
#  in which the package lies
#

setenv PACKAGE tao

#
#  Check the arguments
#

set PRECISION = "double"
@ decode_problem = 0
@ last=$#argv
@ i=1

while ($i <= $last)
  set opt=$argv[$i]
  if("$opt" == '-s')then
    echo "TAO is not available in single precision"
    echo "rerun without -s"
    exit 1
  else if("$opt" == '-decode' ) then
    @ decode_problem = 1
  else if("$opt" == '-h' || "$opt" == '--help' )then
    $MYCUTER/bin/helpmsg
    exit 0
  endif
  @ i++
end

#
#  define the system libraries needed by the package
#  using the format -lrary to include library.a
#

setenv SYSLIBS  "-Wl,-rpath,$TAO_DIR/lib/lib$BOPT/$PETSC_ARCH -Wl,-rpath,$PETSC_DIR/lib/lib$BOPT/$PETSC_ARCH -L$TAO_DIR/lib/lib$BOPT/$PETSC_ARCH -ltaofortran -ltao -L$PETSC_DIR/lib/lib$BOPT/$PETSC_ARCH -lpetscfortran -lpetscsnes -lpetscsles -lpetscdm -lpetscmat -lpetscvec -lpetsc $MPI_LIB"

#  define the name(s) of the object file(s) (files of the form *.o)
#  and/or shared-object libraries (files of the form lib*.so)
#  and/or libraries (files of the form lib*.a) for the GEN package.
#  Object files must lie in 
#    $MYCUTER/(precision)/bin  
#  while libraries and/or shared-object libraries must be in 
#    $MYCUTER/(precision)/lib 
#  where (precision) is either single (when the -s flag is present) or 
#  double (when there is no -s flag)
#

setenv PACKOBJ ""

#
#  define the name of the package specification file (if any)
#  (this possibly precision-dependent file must either lie in
#  the current directory or in $MYCUTER/common/src/pkg/$PACKAGE/ )
#

setenv SPECS ""

#  this is where package-dependent commands (using the commands defined
#  in the current system.(your system) file) may be inserted prior to 
#  decoding the problem file


#  decode the problem file

#{sifdecode}
if( $decode_problem == 1 ) then
  if( $?MYSIFDEC ) then
     $MYSIFDEC/bin/sifdecode $argv
  else
     echo " ${caller} : environment variable MYSIFDEC not set"
     echo "      Either SifDec is not installed or you"
     echo "      should properly set MYSIFDEC"
     exit 7
  endif
endif

if ( $status != 0 ) exit $status

#  this is where package-dependent commands (using the commands defined
#  in the current system.(your system) file) may be inserted prior to 
#  running the package


#  run the package, removing -decode option if present

if( $decode_problem == 1 ) then
  @ n = ${#argv} - 1
  @ i = 1
  set arguments = ''
  while( $i <= $n )
    if( "$argv[$i]" != '-decode' ) then
      set arguments = ( "$arguments" "$argv[$i]" )
    endif
    @ i++
  end
else
  set arguments = "$argv"
endif

set CURDIR = $PWD
cd ${CUTER}/common/src/pkg/tao
$MAKE -s taoma.o
cd $CURDIR

$MYCUTER/bin/runpackage $arguments -n

if ( $status != 0 ) exit $status

#  this is where package-dependent commands (using the commands defined
#  in the current system.(your system) file) may be inserted to clean
#  up after running the package
