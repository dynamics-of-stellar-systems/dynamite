#!/bin/bash
# Definitions for the DFO package
# D. Orban, June 3, 2009.

# Define a short acronym for the package
export PACK=dfo

# Subdirectory of ${CUTER}/common/src/pkg where the package lives
export PACKAGE=dfo
	
# Precision for which the package was written
# Valid values are "single", "double", "single double" and "double single"
PACK_PRECISION="double"

# Define the names of the object files for the package which must lie in
# ${MYCUTER}/(precision)/bin
export PACKOBJS="cuterinter.o ranlux.o"

# Define package and system libraries using -lrary to include library.a or .so.
# No need to mention the BLAS.
export PACKLIBS="-ldfo_ipopt `cat ${MYCUTER}/double/lib/ipopt.liblist`"

# Define the name of the package specification file if any. This possibly
# precision-dependent file must either lie in the current directory or in
# ${CUTER}/common/src/pkg/${PACKAGE}/ )
export SPECS="DFO.SPC"
