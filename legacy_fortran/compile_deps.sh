#!/bin/bash

#FORTRAN=/opt/sw/spack-0.12.1/opt/spack/linux-centos7-x86_64/gcc-4.8.5/gcc-9.1.0-mj7s6dgfnhgi2n42fyxgmitnuslcyol3/bin/gfortran
#HSLARCHIVE=hslarchive-galahad-3.00000.0.tar.gz

[ -v DYNAMITE ] || { echo Need to set DYNAMITE environment variable, i.e. path where you unpacked dynamite ; exit 1 ; }
[ -v FORTRAN ] || { echo Need to set FORTRAN environment variable, e.g. /usr/bin/gfortran ; exit 1 ; }
[ -v HSLARCHIVE ] || { echo Need to set HSLARCHIVE environment variable, e.g. hslarchive-galahad-3.00000.0.tar.gz ; exit 1 ; }

cd ${DYNAMITE}/legacy_fortran

[ -a $HSLARCHIVE ] && tar -xf $HSLARCHIVE || { echo Did not find hslarchive-galahad in $HSLARCHIVE ; exit 1; }

[ -a galahad_makefile ] || { echo Did not find galahad_makefile ; exit 1 ; }

[ -d hslarchive-galahad ] || { echo Did not find extracted hslarchive-galahad ; exit 1; }

[ -a galahad_makefile ] || { echo Did not find galahad_makefile ; exit 1 ; }

[ -a cutest_makefile ] || { echo Did not find cutest_makefile ; exit 1 ; }

if [ ! -d archdefs ] ; then
    git clone https://github.com/ralna/ARCHDefs archdefs --branch v2.0.4
fi

if [ ! -d cutest ] ; then
    git clone https://github.com/ralna/CUTEst cutest --branch v2.0.3
fi

if [ ! -d sifdecode ] ; then
    git clone https://github.com/ralna/SIFDecode sifdecode --branch v2.0.3
fi

if [ ! -d galahad ] ; then
    git clone https://github.com/ralna/GALAHAD galahad --branch v4.0.0
fi

export CUTEST=${DYNAMITE}/legacy_fortran/cutest
export ARCHDEFS=${DYNAMITE}/legacy_fortran/archdefs
export GALAHAD=${DYNAMITE}/legacy_fortran/galahad

VERSION=pc.lnx.gfo

mkdir -p ${CUTEST}/makefiles
mkdir -p ${CUTEST}/modules/${VERSION}/{double,single}
mkdir -p ${CUTEST}/objects/${VERSION}/{double,single}
mkdir -p ${CUTEST}/packages/${VERSION}/{double,single}

cp -f ${DYNAMITE}/legacy_fortran/cutest_makefile ${CUTEST}/makefiles/${VERSION}

cd ${CUTEST}/src

make -f ../makefiles/${VERSION} FORTRAN=${FORTRAN} || exit 1

mkdir -p ${GALAHAD}/{versions,objects,modules,makefiles}
mkdir -p ${GALAHAD}/modules/${VERSION}/{double,single}
mkdir -p ${GALAHAD}/objects/${VERSION}/{double,single}

cp -f ${GALAHAD}/src/makedefs/packages.default ${GALAHAD}/src/makedefs/packages

cp -f ${DYNAMITE}/legacy_fortran/galahad_makefile ${GALAHAD}/makefiles/${VERSION}

cd ${GALAHAD}/src

make -f ../makefiles/${VERSION} qp_cutest FORTRAN=${FORTRAN} || exit 1

echo All dependencies compiled\! Make sure to set GALAHADDIR=${GALAHAD} before compiling dynamite



