#!/bin/bash
# Definitions for the TAO package
# N. Gould, D. Orban & Ph. Toint

# Define a short acronym for the package
export PACK=tao

# Subdirectory of ${CUTER}/common/src/pkg where the package lives
export PACKAGE=tao

# Precision for which the package was written
# Valid values are "single", "double", "single double" and "double single"
export PACK_PRECISION="double"

# Define the name of the object files for the package which must lie in
# ${MYCUTER}/(precision)/bin
export PACKOBJS=""

# Define package and system libraries using -lrary to include library.a or .so
export PACKLIBS="-Wl,-rpath,$TAO_DIR/lib/lib$BOPT/$PETSC_ARCH -Wl,-rpath,$PETSC_DIR/lib/lib$BOPT/$PETSC_ARCH -L$TAO_DIR/lib/lib$BOPT/$PETSC_ARCH -ltaofortran -ltao -L$PETSC_DIR/lib/lib$BOPT/$PETSC_ARCH -lpetscfortran -lpetscsnes -lpetscsles -lpetscdm -lpetscmat -lpetscvec -lpetsc $MPI_LIB"

# Define the name of the package specification file if any. This possibly
# precision-dependent file must either lie in the current directory or in
# ${CUTER}/common/src/pkg/${PACKAGE}/ )
export SPECS=""
