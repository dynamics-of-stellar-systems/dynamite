#!/bin/csh -f
# Displays a generic help message for all scripts
#
# N. Gould, D. Orban & Ph. Toint for CUTEr

#{cmds}


# Comment on the first line of output:
# The -show option looks only useful with sifdecode. Moreover,
# you don't set parameters if the problem has already been decoded.
# Finally, if the problem has not been decoded yet, -n is useless.

set cmd_name = `echo $caller | $SED -e 's/^sd//'`

echo " [sd]${cmd_name} [-s] [-h|--help] [-k] [-o i] [-l secs] [-f] [-b]"
echo '               [-a j] [-show] [-param name=value[,name=value...]] [-force]'
echo '               [-debug] [-Lpath/to/lib] [--blas keyword] [--lapack keyword]'
echo '               probname[.SIF]'

echo ' '

echo ' where: options'
echo '            -n     : use the load module f if it exists'
echo '                     (Default: create a new load module)'
echo '            -h     : print this help and stop execution'
echo '            -s     : run the single precision version'
echo '                     (Default: run the double precision version)'
echo '            -k     : keep the load module after use '
echo '                     (Default: delete the load module)'
echo '            -r     : discourage recompilation of the test problem'
echo '                     (Default: recompile the objects)'
echo '            -o     : 0 for silent mode, 1 for brief description of'
echo '                     the stages executed'
echo '                     (Default: -o 0)'
echo '            -l     : limits the cputime to secs seconds'
echo '                     (Default: -l 99999999)'
echo '            -f     : use automatic differentiation in Forward mode'
echo '            -b     : use automatic differentiation in Backward mode'
echo '            -a     : 1 use the older HSL automatic differentiation package AD01'
echo '                     2 use the newer HSL automatic differentiation package AD02'
echo '                     (Default: -a 2)'
echo '         -show     : displays possible parameter settings for'
echo '                     probname[.SIF]. Other options are ignored'
echo '        -param     : cast probname[.SIF] against explicit parameter'
echo '                     settings. Several parameter settings may be'
echo '                     given as a comma-separated list following'
echo '                     -param or using several -param flags. Use'
echo '                     sifdecode -show probname  to view possible settings'
echo '        -force     : Forces the setting of the parameters named'
echo '                     using -param to the given values even if'
echo '                     those values are not predefined in the SIF file.'
echo '        -debug     : links all the libraries, creates the executable'
echo '                     and stop to allow debugging. This option'
echo '                     automatically disables -n, enables -k,'
echo '                     and turns off all compiler options except -g.'
echo '            -L     : Add the path path/to/lib to the library search path'
echo '        --blas     : --blas -lmyblas uses the BLAS library libmyblas.a'
echo '                     instead of the default linpack provided by CUTEr'
echo '                     --blas none does not link in any BLAS library'
echo '      --lapack     : --lapack -lmylapack uses the LAPACK library'
echo '                     libmylapack.a. --lapack none does not link in any'
echo '                     lapack library.'
echo '      probname     : probname[.SIF] is the name of the file containing'
echo '                     the SIF file for the problem of interest.'
