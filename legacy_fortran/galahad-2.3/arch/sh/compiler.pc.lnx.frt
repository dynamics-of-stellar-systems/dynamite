# Fujitsu frt under linux (pc version)
#
#  Fortran compilation and loading
#

FORTRAN='frt'
BASIC='-c'
LIBCMD='-L /usr/i386-glibc21-linux/lib'
MODCMD='-Am -M $MOD -I $MOD'
MVMODS=':'
OPTIMIZATION='-O'
NOOPTIMIZATION='-O0'
DEBUG=
F77='-Fixed -f s -w -Wa,-W'
F90='-f s -Wa,-W'
F95='-f s -Wa,-W'
NOFMAIN=''
CCONDEF='-DFujitsu_frt'
USUAL=
SPECIAL=
F77SUFFIX=f90
F95SUFFIX=f95
TIMER=GEN
BLAS=
LAPACK=
HSL=
METIS=
NOT95=NOT95
BINSHELL=sh
