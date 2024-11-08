#!/bin/bash
# Definitions for the CG-descent package
# N. Gould, D. Orban & Ph. Toint, September 4th, 2004

# Define a short acronym for the package
export PACK=cgd

# Subdirectory of ${CUTER}/common/src/pkg where the package lives
export PACKAGE=cg_descent

# Precision for which the package was written
# Valid values are "single", "double", "single double" and "double single"
export PACK_PRECISION="double"

# Define the names of the object files for the package which must lie in
# ${MYCUTER}/(precision)/bin
export PACKOBJS="cg_descent.o"

# Define package and system libraries using -lrary to include library.a or .so
export PACKLIBS=""

# Define the name of the package specification file if any. This possibly
# precision-dependent file must either lie in the current directory or in
# ${CUTER}/common/src/pkg/${PACKAGE}/ )
export SPECS="CG_DESCENT.SPC"
