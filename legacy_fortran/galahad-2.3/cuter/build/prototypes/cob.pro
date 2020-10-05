#!/bin/csh -f
# cob: script to decode a sif file and then run COBYLA on the output
#  ( Last modified on 23 Dec 2000 at 17:29:56 )
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
#  define a short acronym for the package to which you wish to make an interface
#

setenv caller cob
setenv PAC cob

#
#  define the name of the subdirectory of $CUTER/common/src/pkg
#  in which the package lies
#

setenv PACKAGE cobyla

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
#  Check that the requested precision is available
#

if ( $PRECISION == "double" ) then
    echo " ERROR: $PACKAGE is not available in $PRECISION precision. "
    echo "        Rerun with -s."
    exit 1
endif

#
#  define the name of the object file for the package which must lie in
#  $MYCUTER/(precision)/bin, where (precision) is either single or double
#

setenv PACKOBJ cobyla.o

#
#  define the system libraries needed by the package
#  using the format -lrary to include library.a
#

setenv SYSLIBS ""

#
#  define the name of the package specification file (if any)
#  (this possibly precision-dependent file must either lie in
#  the current directory or in $CUTER/common/src/pkg/$PACKAGE/ )
#

setenv SPECS "COBYLA.SPC"

#  this is where package-dependent commands (using the commands defined
#  in the current system.(your system) file) may be inserted prior to 
#  decoding the problem file


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

$MYCUTER/bin/runpackage $arguments -n

if ( $status != 0 ) exit $status

#  this is where package-dependent commands (using the commands defined
#  in the current system.(your system) file) may be inserted to clean
#  up after running the package


