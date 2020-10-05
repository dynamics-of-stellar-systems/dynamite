#!/bin/csh -f
# mx: compile MATLAB MEX-files for CUTEr tools
#
# Use: mx [-h|--help] [-o i] [-r] [-u] [-debug]
#
# where: options -h : print this help and stop execution
#                -o : 0 for silent mode, 1 for brief description of
#                     the stages executed
#                     (Default: -o 0)
#                -r : discourage recompilation of the test problem
#                     (Default: recompile the objects)
#                -u : Compile the gateway interface for unconstrained
#                     optimization. (Default: constrained optimization)
#            -debug : link all the libraries, create the executable
#                     and stop to allow debugging. This option
#                     automatically turns off all compiler options except -g.

#  I. Bongartz and A. R. Conn, December 1994.

#  updated for CUTEr, D. Orban, January 2002.

#
#  basic system commands
#

#{cmds}

#
# Environment check
#

envcheck
if( $status != 0 ) exit $status

#
#  define a short acronym for the package to which you wish to make an interface
#

setenv caller mx
setenv PAC mx

#
#  define the name of the subdirectory of $CUTER/common/src/pkg
#  in which the package lies
#

setenv PACKAGE mexfile

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
    set PRECISION = "single"
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

setenv SYSLIBS ""

#
#  If there are compiled, library versions of the level-1 blas
#  (basic linear algebra subprograms), set BLAS to a list of
#  names of the object library suffix -lx, where the object library
#  libx.a contains relevant blas. For instance if the blas are
#  shared between object libraries libblas1.a and libblas2.a,
#  BLAS should be set to "-lblas1 -lblas2", noting that those in
#  libblas1.a will take precedence over those in libblas2.a.
#  If compiled blas are unavailable, BLAS should be set to ""
#

#set BLAS=""
set BLAS="-lblas"

#
#  directory for the main executable file
#

set EXEC=$cwd

#
#  variables for each option
#

#
# OUTPUT = 0 (summary output), = 1 (detailed output from decoder)
#

set OUTPUT=0

#
# RECOMPILE = 1 (by default, recompile the test problem)
# RECOMPILE = 0 -> inhibit recompilation
#

set RECOMPILE = 1

#
# DEBUG = 0 (normal execution), 
# DEBUG = 1 (keep the load module and do *not* execute it)
#

set DEBUG = 0

#
# UNCONS=0 (constrained problem), =1 (unconstrained problem)
#

set UNCONS=0

#
#  interpret arguments
#

@ last=$#argv
@ i=1

while ($i <= $last)
  set opt=$argv[$i]
  if("$opt" == '-h' || "$opt" == '--help')then
    echo ' Use: mx [-h|--help] [-o i] [-r] [-u] [-debug]'
    echo ' '
    echo ' where: options -h : print this help and stop execution'
    echo '                -o : 0 for silent mode, 1 for brief description of'
    echo '                     the stages executed'
    echo '                     (Default: -o 0)'
    echo '                -r : inhibit recompilation of the test problem'
    echo '                     (Default: recompile)'
    echo '                -u : Compile the gateway interface for unconstrained'
    echo '                     optimization. (Default: constrained optimization)'
    exit 0
  else if("$opt" == '-o')then
    @ i++
    set OUTPUT=$argv[$i]
  else if("$opt" == '-r')then
    set RECOMPILE = 0
  else if("$opt" == '-u')then
    set UNCONS=1
  else if ("$opt" == '-show') then
    @ i++
  else if ("$opt" == '-param') then
    @ i = `expr $i + 2`
  else if("$opt" == '-debug') then
    set DEBUG = 1
    set MEXFFLAGS="-g"
    set FFLAGS="-g"
  else
#   discard unrecognized options
    echo "mx: command-line flag $opt not recognized -- skipping"
  endif
  @ i++
end

#
#  build MATLAB MEX-files
#

set LIBDIR=$MYCUTER/double/lib

if (! -e $LIBDIR/libcuter.a ) then
  echo ' '
  echo "Object library for tools not in"
  echo " $LIBDIR"
  echo "Terminating execution."
  exit 4
endif

#  if there are no library blas, include the default ones.

if ("$BLAS" == '') then
   set BLAS=${MYCUTER}/double/bin/linpac.o
endif

# ensure that the current test problem has been compiled.

if ($OUTPUT) then
  echo 'compiling the current test problem, if that is necessary ... '
  echo ' '
endif

# See if we wish to recompile the test problem

if ( $RECOMPILE == 1 ) then
    $RM ELFUN.o GROUP.o RANGE.o EXTER.o
endif

# if ELFUN, etc, haven't been compiled, compile them.

foreach i ( ELFUN GROUP RANGE )
   if ( ! -e $i.o ) then
	$COMPILE $FFLAGS $i.f
	if ($status != 0) then
	    exit 1
	endif
   endif
end

set EXTER = ""
if (-e EXTER.f && -z EXTER.f) $RM EXTER.f
if (-e EXTER.f ) then
  $COMPILE $FFLAGS EXTER.f
  if ($status != 0) then
    exit 1
  endif
  if (-e EXTER.o && -z EXTER.o) then
    $RM EXTER.o
  else
    set EXTER = "EXTER.o"
  endif
endif

# compile the MEX-files

if ($OUTPUT) then
  echo ' '
  echo 'compiling the MEX-files... '
  echo ' '
endif

if (-e EXTER.o && -z EXTER.o) $RM EXTER.o
if (-e EXTER.f && -z EXTER.f) $RM EXTER.f

set PROBDEP = ELFUN,GROUP,RANGE

#
# Compile and link MEX file and CUTEr tools
#

## General tools gateway interface

set general = ${MYCUTER}/double/bin/gtools.f

## Include appropriate CUTEr tools gateway interface

if( $UNCONS ) then
    set gateway = ${MYCUTER}/double/bin/utools.f
    set outputName = utools
else
    set gateway = ${MYCUTER}/double/bin/ctools.f
    set outputName = ctools
endif

# Set MEX file name; include process id to avoid overwritings

#set tempName   = `echo $$ | $AWK 'BEGIN{nb=10}{l=length($1); if(l<nb){start=1}else{start=l-nb+1} print substr($1,start,nb)}'`
#set outputName = "cuter_${tempName}_mx"

$MEXFORTRAN $MEXFFLAGS -output $outputName $general $gateway {$PROBDEP}.o $EXTER $SYSLIBS $BLAS -L${LIBDIR} -lcuter

#  this is where package-dependent commands (using the commands defined
#  in the current system.(your system) file) may be inserted to clean
#  up after running the package

