
##
## CUTEr
##
## If you are using the bourne again shell, (bash)
## what you need to add to your .bashrc file
## looks similar to the following.
##

##	This sets the environment variable CUTER
##	(uncomment as appropriate)
## export CUTER="/usr/local/Cuter"
## export CUTER="/usr/share/Cuter"
## export CUTER="/opt/Cuter"
export CUTER="${HOME}/Cuter"

##     This sets the environment variable MYCUTER
##     (uncomment as appropriate)
## export MYCUTER="${HOME}/CUTEr.large.sun.sol.f90"
## export MYCUTER="${HOME}/mycuter"
## export MYCUTER="${HOME}/Cuter4Solaris"
export MYCUTER="${CUTER}/CUTEr.large.sun.sol.f90"

##     This sets the environment variable MASTSIF
##     pointing to your source of SIF problems.
##     (uncomment as appropriate)
## export MASTSIF="${HOME}/Optimisation/Problems/mastsif"
## export MASTSIF="/usr/share/optimization/sif"
export MASTSIF="${CUTER}/common/sif"

##     This updates the environment variable MATLABPATH
export MATLABPATH="${MYCUTER}/double/bin:${MYCUTER}/single/bin:${CUTER}/common/src/matlab"

##     This updates the library path
export LD_LIBRARY_PATH="$MYCUTER/double/lib:$MYCUTER/single/lib:${LD_LIBRARY_PATH}"

##     This updates the search path for manual pages
export MANPATH="${CUTER}/common/man:${MANPATH}"

##     This updates the search path for executables
export PATH="${CUTER}:${MYCUTER}/bin:$PATH"
