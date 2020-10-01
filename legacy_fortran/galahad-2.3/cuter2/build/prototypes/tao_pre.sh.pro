#!/bin/bash
# tao_pre: Initialization script for TAO package
# N. Gould, D. Orban & Ph. Toint

# Make sure PETSc, TAO environment variables are ok, libraries exist
MPI_LIBDIR='/usr/local/lib'
MPI_LIB="-L$MPI_LIBDIR -lmpich"
if [[ ! -e $MPI_LIBDIR/libmpich.a && ! -e $MPI_LIBDIR/libmpich.so ]]; then
    echo "Error: $MPI_LIBDIR/libmpich.a or $MPI_LIBDIR/libmpich.so not found."
    echo "Either install mpich or alter tao script to describe mpi location"
    exit 1
fi

let vars_unset=0
if [[ ${PETSC_DIR+set} != 'set' || ${BOPT+set} != 'set' || ${TAO_DIR+set} != 'set' || ${PETSC_ARCH+set} != 'set' ]]; then
    echo "Error: You must set the environment variables PETSC_DIR, TAO_DIR,"
    echo "BOPT (=g_c++ or O_c++) and PETSC_ARCH (linux, solaris, etc.)"
    let vars_unset=1
endif

echo "PETSC_DIR = $PETSC_DIR"
echo "PETSC_ARCH = $PETSC_ARCH"
echo "TAO_DIR = $TAO_DIR"
echo "BOPT = $BOPT\n"
(( vars_unset )) && exit 1

plibdir=$PETSC_DIR/lib/lib$BOPT/$PETSC_ARCH
for i in ( petsc petscfortran petscsnes petscsles petscdm petscmat petscvec )
do
    if [[ ! -e ${plibdir}/lib$i.so && ! -e ${plibdir}/lib$i.a ]]; then
        echo "Error: petsc library -l$i not found in ${plibdir}"
        exit 1
    fi
done

tlibdir=$TAO_DIR/lib/lib$BOPT/$PETSC_ARCH
for i in ( tao taofortran )
do
    if [[ ! -e ${tlibdir}/lib$i.so && ! -e ${tlibdir}/lib$i.a ]]; then
        echo "Error: tao library -l$i not found in ${tlibdir}"
        exit 1
    fi
done
