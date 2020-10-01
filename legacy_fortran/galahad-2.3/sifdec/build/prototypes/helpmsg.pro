#!/bin/csh -f
# Displays a generic help message for all scripts
#
# N. Gould, D. Orban & Ph. Toint for SifDec

#{cmds}


# Comment on the first line of output:
# The -show option looks only useful with sifdecode. Moreover,
# you don't set parameters if the problem has already been decoded.
# Finally, if the problem has not been decoded yet, -n is useless.

if ( `echo ${caller} | $GREP ^sd` == "${caller}" ) then
    echo " ${caller} [-s] [-h] [-k] [-o i] [-l secs] [-param name=value[,name=value...]] [-force] [-debug] probname[.SIF]"
else if ( "${caller}" == 'sifdecode' ) then
    echo " ${caller} [-s] [-h] [-k] [-o i] [-l secs] [-f] [-b] [-a j] [-show] [-param name=value[,name=value...]] [-force] [-debug] probname[.SIF]"
else
    echo " ${caller} [-s] [-n] [-h] [-k] [-r] [-o i] [-l secs] [-debug]"
endif
echo ' '
echo ' where: options'
if ( `echo ${caller} | $GREP ^sd` != "${caller}" && "${caller}" != 'sifdecode' ) then
    echo '            -n     : use the load module f if it exists'
    echo '                     (Default: create a new load module)'
endif
echo '            -h     : print this help and stop execution'
echo '            -s     : run the single precision version'
echo '                     (Default: run the double precision version)'
echo '            -k     : keep the load module after use '
echo '                     (Default: delete the load module)'
if ( `echo ${caller} | $GREP ^sd` != "${caller}" && "${caller}" != 'sifdecode' ) then
    echo '                -r : discourage recompilation of the test problem'
    echo '                     (Default: recompile the objects)'
endif
echo '            -o     : 0 for silent mode, 1 for brief description of'
echo '                     the stages executed'
echo '                     (Default: -o 0)'
echo '            -l     : limits the cputime to secs seconds'
echo '                     (Default: -l 99999999)'
if ( ${caller} == 'sifdecode' ) then
    echo '            -f     : use automatic differentiation in Forward mode'
    echo '            -b     : use automatic differentiation in Backward mode'
    echo '            -a     : 1 use the older HSL automatic differentiation package AD01'
    echo '                     2 use the newer HSL automatic differentiation package AD02'
    echo '                     (Default: -a 2)'
endif
if( `echo ${caller} | $GREP ^sd` == "${caller}" || "${caller}" == 'sifdecode' ) then
    echo '            -show  : displays possible parameter settings for'
    echo '                     probname[.SIF]. Other options are ignored'
    echo '            -param : cast probname[.SIF] against explicit parameter'
    echo '                     settings. Several parameter settings may be'
    echo '                     given as a comma-separated list following'
    echo '                     -param or using several -param flags. Use'
    echo '                     sifdecode -show probname  to view possible settings'
    echo '            -force : forces recompilation of the test problem. The'
    echo '                     default is to reuse objects, if they exist.'
endif
echo '            -debug : links all the libraries, creates the executable'
echo '                     and stop to allow debugging. This option'
echo '                     automatically disables -n, enables -k,'
echo '                     and turns off all compiler options except -g.'
if ( `echo ${caller} | $GREP ^sd` == "${caller}" || "${caller}" == 'sifdecode' ) then
    echo '          probname   probname[.SIF] is the name of the file containing'
    echo '                     the SIF file for the problem of interest.'

endif
