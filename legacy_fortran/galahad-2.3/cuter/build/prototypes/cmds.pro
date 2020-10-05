#
#  System specifications
#

set SYSTEM       =  "${SYSTEM}"
set KEYSYS       =  "${KEYSYS}"

#
#  Directory for temporary files
#

set TMP          =   "${TMP}"

#
#  basic system commands
#

set MKDIR        =   "${MKDIR}"
set CP           =   "${CP}"
set RM           =   "${RM}"
set MV           =   "${MV}"
set CAT          =   "${CAT}"
set CHMOD        =   "${CHMOD}"
set SED          =   "${SED}"
set LN           =   "${LN}"
set LS           =   "${LS}"
set AR           =   "${AR}"
set RMDIR        =   "${RMDIR}"
set GREP         =   "${GREP}"
set AWK          =   "${AWK}"
set HEAD         =   "${HEAD}"
set TAIL         =   "${TAIL}"
set WC           =   "${WC}"
set MAKE         =   "${MAKE}"

#
#  Fortran compilation and loading
#

set COMPILE      = "${COMPILE}"
set LOAD         = "${LOAD}"
set ISF9095      = "${ISF9095}"
set COMPILE9095  = "${COMPILE9095}"
set LOAD9095     = "${LOAD9095}"
set FFLAGS       = "${FFLAGS}"
set SPECIALFLAGS = "${SPECIALFLAGS}"

#
#  C compilation and loading
#

set CCOMPILE     = "${CCOMPILE}"
set CLOAD        = "${CLOAD}"
set CFLAGS       = "${CFLAGS}"
set SPECIALLIBS  = "${SPECIALLIBS}"

##
## Note: Recent versions of Matlab use a single script, mex,
##       to compile both C and Fortran code. This makes the
##       commands cmex and fmex obsolete. Change the
##       following lines according to your system specifications.
##
 
set MEXFORTRAN   = "${MEXFORTRAN}"
set MEXFFLAGS    = "${MEXFFLAGS}"
