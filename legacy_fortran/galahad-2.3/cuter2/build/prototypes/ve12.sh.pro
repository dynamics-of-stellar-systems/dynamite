#!/bin/bash
# Definitions for the VE12 package
# N. Gould, D. Orban & Ph. Toint

# Define a short acronym for the package
export PACK=ve12

# Subdirectory of ${CUTER}/common/src/pkg where the package lives
export PACKAGE=hsl_ve12

# Precision for which the package was written
# Valid values are "single", "double", "single double" and "double single"
export PACK_PRECISION="single double"

# Define the name of the object files for the package which must lie in
# ${MYCUTER}/(precision)/bin
export PACKOBJS="cuterinter.o readin.o"

# Define package and system libraries using -lrary to include library.a or .so
export PACKLIBS="-lve12"

# Define the name of the package specification file if any. This possibly
# precision-dependent file must either lie in the current directory or in
# ${CUTER}/common/src/pkg/${PACKAGE}/ )
export SPECS="VE12.SPC"
