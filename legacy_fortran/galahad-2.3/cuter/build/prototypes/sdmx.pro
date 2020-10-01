#!/bin/csh -f
# sdmx: script to decode a SIF file and then create MATLAB MEX-files
#        for the CUTEr tools
#
# Use: sdmx [-h|--help] [-o j] [-r] [-u] [-debug] probname
#
# where: options -h : print this help and stop execution
#                -o : 0 for silent mode, 1 for brief description of
#                     the stages executed.
#                     (Default: -o 0)
#                -u : create MEX-file for unconstrained tools
#                     (Default: create MEX-file for constrained tools)
#            -debug : link all the libraries, create the executable
#                     and stop to allow debugging. This option
#                     automatically turns off all compiler options except -g.
#
#        probname     probname.SIF is the name of the file containing
#                     the SIF file for the problem of interest.
#

#  I. Bongartz and A. R. Conn, December 1994.

#  updated for CUTEr, D. Orban, January 2002.

#
#  basic system commands
#

#{cmds}

unset noclobber

#
# Environment check
#

envcheck
if( $status != 0 ) exit $status

#
#  define a short acronym for the package to which you wish to make an interface
#

setenv caller sdmx
setenv PAC mx

#
#  define the name of the subdirectory of $CUTER/common/src/pkg
#  in which the package lies
#

setenv PACKAGE mxfile

#
#  define the system libraries needed by the package
#  using the format -lrary to include library.a
#

setenv SYSLIBS ""

#
#  variables for each option
#

#
# OUTPUT=0 (summary output), = 1 (detailed output from decoder)
#

set OUTPUT=0

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

set METHOD=3

@ last = $#argv

if ($last == 0) then
  echo 'Use: sdmx [-h|--help] [-o j] [-u] [-show] [-param name=value[,name=value...]] [-debug] probname'
  exit 1
endif

@ i=1

while ($i <= $last)
  set opt=$argv[$i]
  if("$opt" == '-h' || "$opt" == "--help")then
    echo ' Use: sdmx [-h|--help] [-o j] [-u] [-show] [-param name=value[,name=value...]] [-force] [-debug] probname'
    echo ' '
    echo ' where: options -h : print this help and stop execution'
    echo '                -o : 0 for silent mode, 1 for brief description of'
    echo '                     the stages executed'
    echo '                     (Default: -o 0)'
    echo '                -u : create MEX-file for unconstrained tools'
    echo '                     (Default: create MEX-file for constrained tools)'
    echo '            -show  : displays possible parameter settings for'
    echo '                     probname[.SIF]. Other options are ignored'
    echo '            -param : cast probname[.SIF] against explicit parameter'
    echo '                     settings. Several parameter settings may be'
    echo '                     given as a comma-separated list following'
    echo '                     -param or using several -param flags. Use'
    echo '                     sifdecode -show probname  to view possible settings'
    echo '            -force : Forces the setting of the parameters named'
    echo '                     using -param to the given values even if'
    echo '                     those values are not predefined in the SIF file.'
    echo '            -debug : link all the libraries, create the executable'
    echo '                     and stop to allow debugging. This option'
    echo '                     automatically turns off all compiler options except -g.'
    echo ' '
    echo '        probname     probname.SIF is the name of the file containing'
    echo '                     the SIF file for the problem of interest.'
    exit 0
  else if("$opt" == '-o')then
    @ i++
    set OUTPUT=$argv[$i]
#  else if("$opt" == '-u')then
#    @ i++
#    set UNCONS=1
#  else if("$opt" == '-m') then
#    @ i++
#    set METHOD = $argv[$i]
#  else if("$opt" == '-debug') then
#    set DEBUG = 1
  endif
  @ i++
end

if ($OUTPUT) then
  echo 'Convert SIF file into data and routines before compiling MEX-files...'
  echo ' '
  echo 'Problem details will be given.'
  echo ' '
endif

if( $?MYSIFDEC ) then
	$MYSIFDEC/bin/sifdecode $argv
	if ( $status != 0 ) exit $status
else
	echo " ${caller} : environment variable MYSIFDEC not set"
	echo "      Either SifDec is not installed or you"
	echo "      should properly set MYSIFDEC"
	exit 7
endif

if ($OUTPUT) echo ' '

@ n = ${#argv} - 1
$MYCUTER/bin/mx $argv[1-$n]

$RM $TMP/sdmx.input

if ( $status != 0 ) exit $status

#  this is where package-dependent commands (using the commands defined
#  in the current system.(your system) file) may be inserted to clean
#  up after running the package


