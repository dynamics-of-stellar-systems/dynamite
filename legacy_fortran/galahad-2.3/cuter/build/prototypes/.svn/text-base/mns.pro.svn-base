#!/bin/csh -f
#  ( Last modified on 23 Dec 2000 at 17:29:56 )
#
# mns: apply MNS to a problem and delete the executable after use.
#
#{version}
#
# Use: mns [-n] [-h] [-s] [-k] [-r] [-o i] [-l secs] [-debug]
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

setenv caller mns
setenv PAC mns

#
#  define the name of the subdirectory of $CUTER/common/src/pkg
#  in which the package lies
#

setenv PACKAGE mns

#
#  define the system libraries needed by the package
#  using the format -lrary to include library.a
#

setenv SYSLIBS ""

#
#  Check the arguments
#

set PRECISION = "double"
@ last=$#argv
@ i=1

while ($i <= $last)
  set opt=$argv[$i]
  if("$opt" == '-s')then
    set PRECISION = "single"
  else if("$opt" == '-h' || "$opt" == '--help' )then
    $MYCUTER/bin/helpmsg
    exit 0
  endif
  @ i++
end

#
#  Check that the requested precision is available
#

if ( $PRECISION == "single" ) then
    echo " ERROR: $PACKAGE is not available in $PRECISION precision. "
    echo "        Rerun without -s."
    exit 1
endif

#
#  define the name of the object file for the package which must lie in
#  $MYCUTER/(precision)/bin, where (precision) is either single or double
#

setenv PACKOBJ minos.o

#
#  define the name of the package specification file (if any)
#  (this possibly precision-dependent file must either lie in
#  the current directory or in $CUTER/common/src/pkg/$PACKAGE/ )
#

setenv SPECS "MINOS.SPC"

#
#  open the MINOS specs and basis files, using the default specs if necessary
#

$LN -fs MINOS.SPC fort.4
if ( -e fort.9 ) $RM fort.9
if ( -e MINOS.BASIS ) $LN -sf MINOS.BASIS fort.11
if ( -e fort.12 ) $RM fort.12

#  run the package

$MYCUTER/bin/runpackage $argv

if ( $status != 0 ) exit $status

#  tidy up the current directory, deleting all junk.

if ( -e fort.4  ) $RM fort.4
if ( -e fort.9  ) $MV fort.9  MINOS.LIS
if ( -e fort.11 ) $RM fort.11
if ( -e fort.12 ) $MV fort.12 MINOS.NEWBASIS
