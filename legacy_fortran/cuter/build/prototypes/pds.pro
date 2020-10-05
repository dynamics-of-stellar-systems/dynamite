#!/bin/csh -f
#  ( Last modified on 23 Dec 2000 at 17:29:56 )
#
# pds: apply PDS to a problem and delete the executable after use.
#
#{version}
#
# Use: pds [-n] [-h] [-s] [-k] [-r] [-o i] [-l secs] [-debug]
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
#                         automatically disables -n and enables -k.

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

setenv caller pds
setenv PAC pds

#
#  define the name of the subdirectory of $CUTER/common/src/pkg
#  in which the package lies
#

setenv PACKAGE pds

#
#  define the system libraries needed by the package
#  using the format -lrary to include library.a
#

setenv SYSLIBS ""

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
    echo 'IPOPT is only available in double precision'
    echo 'rerun without -s'
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
#  define the name(s) of the object file(s) for the package which must lie in
#  $MYCUTER/(precision)/bin, where (precision) is either single or double
#

setenv PACKOBJ "pds.o"

#
#  define the name of the package specification file (if any)
#  (this possibly precision-dependent file must either lie in
#  the current directory or in $MYCUTER/common/src/pkg/$PACKAGE/ )
#

setenv SPECS "PDS.SPC"

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
