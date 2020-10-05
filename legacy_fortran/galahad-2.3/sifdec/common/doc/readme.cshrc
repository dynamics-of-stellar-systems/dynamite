
##
## SifDec
##
## If you are using the C shell (csh) or the
## enhanced C shell (tcsh),
## what you need to add to your .cshrc file
## looks similar to the following.
##

##	This sets the environment variable SIFDEC
##	(uncomment as appropriate)
## setenv SIFDEC /usr/local/Sifdec
## setenv SIFDEC /usr/share/Sifdec
## setenv SIFDEC /opt/Sifdec
setenv SIFDEC ${HOME}/Sifdec

##     This sets the environment variable MYSIFDEC
##     (uncomment as appropriate)
## setenv MYSIFDEC ${HOME}/SifDec.large.sun.sol.f90
## setenv MYSIFDEC ${HOME}/mysifdec
## setenv MYSIFDEC ${HOME}/Sifdec4Solaris
setenv MYSIFDEC  ${SIFDEC}/SifDec.large.sun.sol.f90

##     This sets the environment variable MASTSIF
##     pointing to your source of SIF problems.
##     (uncomment as appropriate)
## setenv MASTSIF ${HOME}/Optimisation/Problems/mastsif
## setenv MASTSIF /usr/share/optimization/sif
setenv MASTSIF ${SIFDEC}/common/sif

##     This updates the search path for manual pages
setenv MANPATH ${SIFDEC}/common/man:${MANPATH}

##     This updates the search path for executables
set path=(${SIFDEC} ${MYSIFDEC}/bin $path)
