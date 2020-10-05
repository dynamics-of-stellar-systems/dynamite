
##
## CUTEr
##
## If you are using the C shell (csh) or the
## enhanced C shell (tcsh),
## what you need to add to your .cshrc file
## looks similar to the following.
##

##	This sets the environment variable CUTER
##	(uncomment as appropriate)
## setenv CUTER /usr/local/Cuter
## setenv CUTER /usr/share/Cuter
## setenv CUTER /opt/Cuter
setenv CUTER ${HOME}/Cuter

##     This sets the environment variable MYCUTER
##     (uncomment as appropriate)
## setenv MYCUTER ${HOME}/CUTEr.large.sun.sol.f90
## setenv MYCUTER ${HOME}/mycuter
## setenv MYCUTER ${HOME}/Cuter4Solaris
setenv MYCUTER  ${CUTER}/CUTEr.large.sun.sol.f90

##     This sets the environment variable MASTSIF
##     pointing to your source of SIF problems.
##     (uncomment as appropriate)
## setenv MASTSIF ${HOME}/Optimisation/Problems/mastsif
## setenv MASTSIF /usr/share/optimization/sif
setenv MASTSIF ${CUTER}/common/sif

##     This updates the environment variable MATLABPATH
setenv MATLABPATH ${MYCUTER}/double/bin:${MYCUTER}/single/bin:${CUTER}/common/src/matlab

##     This updates the library path
setenv LD_LIBRARY_PATH $MYCUTER/double/lib:$MYCUTER/single/lib:${LD_LIBRARY_PATH}

##     This updates the search path for manual pages
setenv MANPATH ${CUTER}/common/man:${MANPATH}

##     This updates the search path for executables
set path=(${CUTER} ${MYCUTER}/bin $path)
