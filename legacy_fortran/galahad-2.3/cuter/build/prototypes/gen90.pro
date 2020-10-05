#!/bin/csh -f
#  ( Last modified on 23 Dec 2000 at 17:29:56 )
#
# gen90: apply GEN90 to a problem and delete the executable after use.
#
#{version}
#

#{args}

#
#  N. Gould, D. Orban & Ph. Toint, November 7th, 2000
#

#{cmds}

#
# Compiler check
#

if( "${ISF9095}" == "no" ) then
  echo 'You are not currently using a Fortran 90/95 compiler.'
  echo 'Aborting.'
  exit 10
endif

#
# Environment check
#

envcheck
if( $status != 0 ) exit $status

#
#  define a short acronym for the package to which you wish to make an interface
#

setenv caller gen90
setenv PAC gen90

#
#  define the name of the subdirectory of $CUTER/common/src/pkg
#  in which the package lies
#

setenv PACKAGE gen

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

setenv PACKOBJ "cuterinter.o gen90.o"

#
#  define the name of the package specification file (if any)
#  (this possibly precision-dependent file must either lie in
#  the current directory or in $MYCUTER/common/src/pkg/$PACKAGE/ )
#

setenv SPECS "GEN.SPC"

#  this is where package-dependent commands (using the commands defined
#  in the current system.(your system) file) may be inserted prior to 
#  running the package


#  decode the problem file

#{sifdecode}
if( $decode_problem == 1 ) then
  if( $?MYSIFDEC ) then
     $MYSIFDEC/bin/sifdecode $argv
     if ( $status != 0 ) exit $status
  else
     echo " ${caller} : environment variable MYSIFDEC not set"
     echo "      Either SifDec is not installed or you"
     echo "      should properly set MYSIFDEC"
     exit 7
  endif
endif

#  See if we need to recompile the main program

set MYCUTERPREC = $MYCUTER/$PRECISION
set MYCUTERPRECBIN = ${MYCUTERPREC}/bin

set PACKOBJLIST = ( $PACKOBJ )
set PACKOBJLIST = `echo $PACKOBJLIST:q`
@ no = ${#PACKOBJLIST}
set PACKLAST = $PACKOBJLIST[$no]

setenv PACKOBJ "$PACKOBJ"

set PACKTYPE = $PACKLAST:e
if ( $PACKTYPE == "o" ) then
  set PACKFULL = $MYCUTERPRECBIN/$PACKLAST
else
  set PACKFULL = $MYCUTERPREC/lib/$PACKLAST
endif

if ( ! -e $PACKFULL ) then
  echo " There is no package binary file"
  echo "  $PACKFULL"
  echo " Please, (re)compile the $PACKAGE package and try again"
  echo " Terminating execution."
  exit 5
endif

if ( -e $MYCUTERPRECBIN/${PAC}ma.o ) then
  if ( -e $CUTER/common/src/tools/${PAC}ma.f90 && -e $PACKFULL ) then
    set newer = (`ls -ct $PACKFULL $MYCUTERPRECBIN/${PAC}ma.o`)
    if ( $newer[1] == "$PACKFULL" ) then
      set recompile = 'true'
    else
      set recompile = 'false'
    endif
  else
    if (! -e $CUTER/common/src/tools/${PAC}ma.f90 ) \\
      echo " File $CUTER/common/src/tools/${PAC}ma.f90 is missing"
    if (! -e $PACKFULL ) echo " File $PACKFULL is missing"
      exit 1
  endif
else
  set recompile = 'true'
endif

#  we need to recompile the main program

if ( $recompile == 'true' ) then
   echo " ... Recompiling ${PAC}ma.f90"
   echo " "
   set THISDIR = $PWD
   cd $CUTER/common/src/tools
   set FILE = ${PAC}ma
   set FILEF = "${FILE}.f90"
   $SED -f $MYCUTERPREC/config/cast90.sed $FILEF > $MYCUTERPRECBIN/$FILEF
   cd $MYCUTERPRECBIN
   $COMPILE9095 -o ./$FILE.o $FILEF
   if ( $status != 0 ) exit $status
   $RM $MYCUTERPREC/bin/$FILEF
   cd $THISDIR
endif

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

$MYCUTER/bin/runpackage $arguments -n

if ( $status != 0 ) exit $status

#  this is where package-dependent commands (using the commands defined
#  in the current system.(your system) file) may be inserted to clean
#  up after running the package
