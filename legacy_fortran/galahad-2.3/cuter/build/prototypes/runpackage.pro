#!/bin/csh -f
# runpackage: apply package PAC to a problem and delete the executable 
# after use.
#  ( Last modified on 23 Dec 2000 at 17:29:56 )
#
#
#{version}
#
# Use: runpackage [-n] [-h] [-s] [-k] [-r] [-o i] [-l secs] [-debug]
#
# where: options -n : use the load module f if it exists
#                     (Default: create a new load module)
#                -h : print this help and stop execution
#                -s : run the single precision version
#                     (Default: run the double precision version)
#                -k : keep the load module after use
#                     (Default: delete the load module)
#                -r : discourage recompilation of the test problem
#                     (Default: recompile the objects)
#                -o : 0 for silent mode, 1 for brief description of
#                     the stages executed
#                     (Default: -o 0)
#                -l : limit the cputime used to secs seconds
#                     (Default: -l 99999999)
#                -debug : links all the libraries, creates the executable
#                         and stop to allow debugging. This option
#                         automatically disables -n, enables -k,
#                         and turns off all compiler options except -g.

#
#  N. Gould, D. Orban & Ph. Toint, November 7th, 2000
#

#{cmds}

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

set BLAS=""
#set BLAS="-lblas"
@ blas_set = 0

set LAPACK = ""

#
#  directory for the main executable file
#

set EXEC=$cwd

#
#  variables for each option
#

#
# PRECISION = single (single precision), = double (double precision)
#

set PRECISION="double"

#
# NEW = 0 (run existing f module), = 1 (build a new module)
#

set NEW = 0

#
# KEEP = 0 (discard f load module after use), = 1 (keep it)
#

set KEEP = 0

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
# OUTPUT = 0 (summary output), = 1 (detailed output from decoder)
#

set OUTPUT=0

#
# LIMIT = 0 (no cputime limit)
#

set LIMIT = 99999999

#
# Specify alternate library paths
#

set ALT_LIB_PATH = ''

#
#  interpret arguments
#

#
#  check for residual -show and/or -param arguments
#  and simply skip them
#

@ last=$#argv
@ i=1

while ($i <= $last)
  set opt=$argv[$i]
  if ("$opt" == '-n') then
    set NEW=1
  else if ("$opt" == '-s') then
    set PRECISION="single"
  else if ("$opt" == '-h' || "$opt" == '--help' ) then
    $MYCUTER/bin/helpmsg
    exit 0
  else if ("$opt" == '-k') then
    set KEEP=1
  else if ("$opt" == '-r') then
    set RECOMPILE=0
  else if ("$opt" == '-o') then
    @ i++
    set OUTPUT=$argv[$i]
  else if ("$opt" == '-l') then
    @ i++
    set LIMIT=$argv[$i]
  else if ("$opt" == '-c') then
    set FFLAGS = "$FFLAGS $SPECIALFLAGS"
  else if ("$opt" == '--blas' ) then
    # Note: This allows to set BLAS to nothing
    @ i++
    @ blas_set = 1
    if ( "$argv[$i]" != 'none' ) set BLAS = $argv[$i]
  else if ("$opt" == '--lapack' ) then
    @ i++
    if ( "$argv[$i]" != 'none' ) set LAPACK = $argv[$i]
  else if ("$opt" == '-debug') then
    set NEW=1
    set KEEP=1
    set DEBUG=1
    # remove all compiler flags and compile using -g
    set FFLAGS = "-g"
  else if ("$opt" == '-show') then
    # Skip this option
  else if ("$opt" == '-param') then
    # Skip this option and its argument
    @ i++
  else if ( `echo $opt | egrep '^\-L'` == "$opt" ) then
    set ALT_LIB_PATH = ( "$ALT_LIB_PATH" "$opt" )
  else
    echo " $0:t : option $opt not understood -- discarding"
#    echo " Use: ${PAC} [-n] [-h] [-s] [-k] [-r] [-o i] [-l secs] [-debug]"
#    exit 1
  endif
  @ i++
end

#
#  run PAC without rebuilding it
#

if (! $NEW) then
  if (! -e $EXEC/${PAC}min || ! -x $EXEC/${PAC}min) then
    echo ' '
    echo "load module ${PAC}min not found/executable. Rerun with -n option"
    echo ' '
    exit 3
  endif
  if ($OUTPUT) then
    echo ' '
    echo "running ${PAC}min on current test problem ... "
    echo ' '
  endif
  limit cputime $LIMIT
  $EXEC/${PAC}min

#  tidy up the current directory, deleting all junk.

  if (! $KEEP) $RM $EXEC/${PAC}min
  exit 0
endif

#
#  SUBR = Object file for the package subroutine(s)
#

set SUBR = ""
foreach i ( $PACKOBJ )
  if (-e $MYCUTER/$PRECISION/bin/$i) then
    set SUBR = "$SUBR $MYCUTER/$PRECISION/bin/$i"
  else if (-e $MYCUTER/$PRECISION/lib/$i) then
    set SUBR = "$SUBR -L$MYCUTER/$PRECISION/lib -l`echo $i:r | $SED 's/^lib//'`"
  else
    echo ' '
    echo "$i object file not in "
    echo " $MYCUTER/$PRECISION/bin/"
    echo "Terminating execution."
    exit 5
  endif
end

#
#  LIBDIR = CUTEr Library directory
#

set LIBDIR=$MYCUTER/$PRECISION/lib

if (! -e $LIBDIR/libcuter.a ) then
  echo ' '
  echo "Object library for tools not in"
  echo " $LIBDIR"
  echo "Terminating execution."
  exit 4
endif

#
#  DRIVER = main program driving the package
#

set DRIVER=$MYCUTER/$PRECISION/bin/${PAC}ma.o

#
#  If needed, find the correct BLAS and package spec file
#

if( $blas_set == 0 ) set BLAS=$MYCUTER/$PRECISION/bin/linpac.o
if ( $SPECS != "") then
  if (! -e $SPECS ) then
     if (-e $MYCUTER/$PRECISION/specs/$SPECS ) then
       $CP $MYCUTER/$PRECISION/specs/$SPECS $SPECS
     else
       $CP $CUTER/common/src/pkg/$PACKAGE/$SPECS $SPECS
     endif
  endif
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

# if ELFUN, etc, have not been compiled, compile them.

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

# link all the PAC and tools files together.

if ($OUTPUT) then
  echo ' '
  echo 'linking all the object files together ... '
  echo ' '
endif

$LOAD $FFLAGS -o ${PAC}min ELFUN.o GROUP.o RANGE.o $EXTER $DRIVER \
      $SUBR $ALT_LIB_PATH -L$LIBDIR $SYSLIBS $SPECIALLIBS -lcuter $BLAS $LAPACK

if ( $cwd != $EXEC ) $MV ${PAC}min $EXEC/${PAC}min

#  run PAC on the current test problem unless -debug is set.

if (! $DEBUG) then

    if ($OUTPUT) then
    echo ' '
    echo "running ${PAC}min on current test problem ... "
    echo ' '
    endif

    limit cputime $LIMIT
    $EXEC/${PAC}min

#     tidy up the current directory, deleting all junk.

    if (! $KEEP) then
	$RM $EXEC/${PAC}min
	$RM ELFUN.o GROUP.o RANGE.o EXTER.o
    endif

else
    
    echo '  debug enabled, load module is in '$EXEC'/'

endif

