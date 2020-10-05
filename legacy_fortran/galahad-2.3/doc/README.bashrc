
##
## GALAHAD
##
## If you are using the bourne again shell, (bash)
## what you need to add to your .bashrc file
## looks similar to the following.
##

##	This sets the environment variable GALAHAD
##	(uncomment as appropriate)
## export GALAHAD="/usr/local/galahad"
## export GALAHAD="/usr/share/galahad"
## export GALAHAD="/opt/galahad"
export GALAHAD="${HOME}/galahad"

##     This updates the search path for manual pages
export MANPATH="${GALAHAD}/man:${MANPATH}"

##     This updates the search path for executables
export PATH="${GALAHAD}/bin:$PATH"

##	This sets the environment variable AMPLDIR
##	(uncomment as appropriate)
## export AMPLDIR=/usr/local/ampldir
## export AMPLDIR=/usr/share/ampldir
## export AMPLDIR=/opt/ampldir
export AMPLDIR=${HOME}/ampldir
