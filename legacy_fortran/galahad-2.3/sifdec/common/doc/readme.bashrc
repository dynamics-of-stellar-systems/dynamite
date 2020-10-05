
##
## SifDec
##
## If you are using the bourne again shell, (bash)
## what you need to add to your .bashrc file
## looks similar to the following.
##

##	This sets the environment variable SIFDEC
##	(uncomment as appropriate)
## export SIFDEC="/usr/local/Sifdec"
## export SIFDEC="/usr/share/Sifdec"
## export SIFDEC="/opt/Sifdec"
export SIFDEC="${HOME}/Sifdec"

##     This sets the environment variable MYSIFDEC
##     (uncomment as appropriate)
## export MYSIFDEC="${HOME}/SifDec.large.sun.sol.f90"
## export MYSIFDEC="${HOME}/mysifdec"
## export MYSIFDEC="${HOME}/Sifdec4Solaris"
export MYSIFDEC="${SIFDEC}/SifDec.large.sun.sol.f90"

##     This sets the environment variable MASTSIF
##     pointing to your source of SIF problems.
##     (uncomment as appropriate)
## export MASTSIF="${HOME}/Optimisation/Problems/mastsif"
## export MASTSIF="/usr/share/optimization/sif"
export MASTSIF="${SIFDEC}/common/sif"

##     This updates the search path for manual pages
export MANPATH="${SIFDEC}/common/man:${MANPATH}"

##     This updates the search path for executables
export PATH="${SIFDEC}:${MYSIFDEC}/bin:$PATH"
