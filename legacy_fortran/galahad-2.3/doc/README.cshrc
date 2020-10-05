
##
## GALAHAD
##
## If you are using the C shell (csh) or the
## enhanced C shell (tcsh),
## what you need to add to your .cshrc file
## looks similar to the following.
##

##	This sets the environment variable GALAHAD
##	(uncomment as appropriate)
## setenv GALAHAD /usr/local/galahad
## setenv GALAHAD /usr/share/galahad
## setenv GALAHAD /opt/galahad
setenv GALAHAD ${HOME}/galahad

##     This updates the search path for manual pages
setenv MANPATH ${GALAHAD}/man:${MANPATH}

##     This updates the search path for executables
set path=(${GALAHAD}/bin $path)

##	This sets the environment variable AMPLDIR
##	(uncomment as appropriate)
## setenv AMPLDIR /usr/local/ampldir
## setenv AMPLDIR /usr/share/ampldir
## setenv AMPLDIR /opt/ampldir
setenv AMPLDIR ${HOME}/ampldir

