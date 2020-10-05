#!/bin/bash

# Remco van den Bosch June 2007

# Script to make a Triaxial Schwarzschild model from Abel input
# This script is made for a dual processor machine. (hence the &'s)


#path to the code
dir="/home/lzhu/bin/"    
#dir="/Users/bosch/Desktop/Dropbox/TriaxSchwarzschild_forkMar2011_Lz/TriaxSchwarzschild_25102010/"
#dir=""

# copy additional input files to this direcotory
#cp ${dir}/Abel/infil_template/*.in infil/

export GFORTRAN_UNBUFFERED_ALL=y

mkdir datfil


#define this if you want to use the analytic abel potential
#undefine it to use the MGE potential
#abel="abel" 
  

# Generate the orbital startpoints
touch datfil/orbstart.log
tail -f datfil/orbstart.log &
${dir}orbitstart${abel} < infil/orbstart.in > datfil/orbstart.log

#export GFORTRAN_UNBUFFERED_ALL=n
# make orbitlibrary 
${dir}orblib${abel}   < infil/orblib.in > datfil/orblib.log  &


# make orbitlibrary with boxes
${dir}orblib${abel}   < infil/orblibbox.in > datfil/orblibbox.log

# calculate the intrinsic masses
${dir}triaxmass    < infil/triaxmass.in

# calculate the projected masses
${dir}triaxmassbin < infil/triaxmassbin.in

export GFORTRAN_UNBUFFERED_ALL=y
#solve with nnls 
${dir}triaxnnls    < infil/nn.in > datfil/nn.log  &
#triaxnnls    < infil/nnr4.in > datfil/nnr4.log




               
               
