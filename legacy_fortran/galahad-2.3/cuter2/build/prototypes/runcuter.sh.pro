#!/bin/bash
###############################################################################
#
# runcuter: the new and improved runpackage
#
# N. Gould, D. Orban & Ph. Toint for GOT Productions, 2006.
#
###############################################################################

# Default settings

Initialize_Settings() {

    # Package to run
    PKG=''
    let package_set=0

    # Problem to decode
    PROBLEM=''
    let problem_set=0

    # If there are compiled, library versions of the level-1 blas (basic linear
    # algebra subprograms), set BLAS to a list of names of the object library
    # suffix -lx, where the object library libx.a contains relevant blas. For
    # instance if the blas are shared between object libraries libblas1.a and
    # libblas2.a, BLAS should be set to "-lblas1 -lblas2", noting that those in
    # libblas1.a will take precedence over those in libblas2.a. If compiled
    # blas are unavailable, BLAS should be set to ""
    BLAS=""
    let blas_set=0
    LAPACK=""

    # directory for the main executable file
    EXEC=${PWD}

    # PRECISION = single (single precision), = double (double precision)
    PRECISION="double"

    # KEEP = 0 (discard f load module after use), = 1 (keep it)
    let KEEP=0

    # RECOMPILE = 0 (by default, do not recompile the test problem)
    # RECOMPILE = 1 -> force recompilation
    let RECOMPILE=0

    # DEBUG = 0 (normal execution),
    # DEBUG = 1 (keep the load module and do *not* execute it)
    let DEBUG=0

    # PROFILE = 0 (normal execution)
    # PROFILE = 1 (instrument code for profiling)
    let PROFILE=0

    # OUTPUT = 0 (summary output), = 1 (detailed output from decoder)
    let OUTPUT=0

    # LIMIT = 0 (no cputime limit)
    let LIMIT=0

    # Specify alternate library paths
    ALT_LIB_PATH=' '

    # Specify extra libraries
    XTRALIBS=''

    # Specify constrained or unconstrained problem (for Matlab interface)
    # This option is now deprecated
    export constraint_flag=''

    # Specify whether or not to display compilation commands
    let show_commands=0

    # Specify whether or not to display current settings
    let show_config=0

    # Options to pass to sifdecode
    sd_opts=( )
}

###############################################################################

# Display help

Display_Usage() {

    printf "Syntax: ${thisprog} -p pkg [options...]\n"
    printf 'Here is a summary of the available options.\n'
    printf 'See the man page for further details.\n\n'
    printf '  -p  or --package \tspecify package name\n'
    printf '  -s  or --single  \trun single precision version of package\n'
    printf '  -h  or --help    \tdisplay this help message\n'
    printf '  -k  or --keep    \tkeep executable after run is over\n'
    printf '  -r  or --rebuild \tforce recompilation of problem subroutines\n'
    printf '  -o  or --output  \tset verbosity level\n'
    printf '  -t  or --limit   \tset a time limit on the run\n'
    printf '  -c  or --cfortran\tcompile mixed C/Fortran source code\n'
    printf '  -b  or --blas    \tspecify BLAS library\n'
    printf '  -K  or --lapack  \tspecify LAPACK library\n'
    printf '  -g  or --debug   \tcompile debug version of executable\n'
    printf '  -pg or --profile \tinstrument code for profiling\n'
    printf '  -D  or --decode  \tspecify problem to decode and compile\n'
    printf '  -u  or --uncons  \tspecify unconstrained problem (Matlab---deprecated)\n'
    printf '  -Lpath           \tadd path to libraries path\n'
    printf '  -lrary           \tlink in library.a or library.so\n'
    printf '        --command \tdisplay compilation commands\n'
    printf '        --config  \tdisplay current settings\n'
    printf 'All other options are passed unchanged to sifdecode\n'
}

###############################################################################

# Read command-line arguments

Parse_Arguments() {

    while [[ $# > 0 ]]; do
        case "$1" in
            -p|--package)  PKG="$2"
                           let package_set=1
                           shift
                           ;;
            -s|--single)   PRECISION="single"
                           ;;
            -h|--help)     Display_Usage
                           exit 0
                           ;;
            -k|--keep)     let KEEP=1
                           ;;
            -r|--rebuild)  let RECOMPILE=1
                           ;;
            -o|--output)   OUTPUT=$2
                           shift
                           ;;
            -t|--limit)    let LIMIT=$2
                           shift
                           ;;
            -c|--cfortran) FFLAGS="${FFLAGS[@]} ${SPECIALFLAGS[@]}"
                           # Used to mix C and Fortran code
                           ;;
            -b|--blas)     [[ "$2" != 'none' ]] && BLAS="$2"
                           # Note: This allows to set BLAS to nothing
                           let blas_set=1
                           shift
                           ;;
            -K|--lapack)   [[ "$2" != 'none' ]] && LAPACK="$2"
                           shift
                           ;;
            -g|--debug)    let NEW=1
                           let KEEP=1
                           let DEBUG=1
                           FFLAGS=${FDEBUGFLAGS[@]}
                           ;;
            -pg|--profile) let NEW=1
                           let PROFILE=1
                           FFLAGS="${FFLAGS[@]} ${FPROFILEFLAGS[@]}"
                           ;;
            -D|--decode)   PROBLEM="$2"
                           let problem_set=1
                           shift
                           ;;
            -u|--uncons)   export constraint_flag='unconstrained'
                           ;;
            -L*)           ALT_LIB_PATH=( "${ALT_LIB_PATH[@]}" "$1" )
                           ;;
            -l*)           XTRALIBS=( "${XTRALIBS[@]}" "$1" )
                           ;;
            --command)     let show_commands=1
                           ;;
            --config)      let show_config=1
                           ;;
            *)             # Pass this option unchanged to sifdecode
                           sd_opts=( ${sd_opts[@]} "$1" )
                           ;;
        esac
        shift
    done
}

###############################################################################

# Decode problem

Decode_Problem() {
    # Problem name is passed as argument #1
    ${RM} ${EXEC}/${PACK}min
    ${RM} ELFUN.f GROUP.f RANGE.f EXTER.f
    ${RM} ELFUN.o GROUP.o RANGE.o EXTER.o
    if [[ ${MYSIFDEC+set} == 'set' ]]; then
        echo "sifdecode ${sd_opts[@]} $1"
        ${MYSIFDEC}/bin/sifdecode ${sd_opts[@]} $1
        [[ $? != 0 ]] && exit $?
    else
        echo "${thisprog}: environment variable MYSIFDEC not set"
        echo "Either SifDec is not installed or MYSIFDEC is incorrectly set"
        exit 7
    fi
}

###############################################################################

# Run specified package

Run_Package() {

    (( OUTPUT )) && printf "\nRunning ${PACK}min on current test problem ...\n"
    (( LIMIT  )) && ulimit -t ${LIMIT}
    ${EXEC}/${PACK}min

    # Tidy up the current directory, deleting all junk.
    (( KEEP == 0 )) && ${RM} ${EXEC}/${PACK}min
}

###############################################################################

# Run pre script

Run_Pre() {

    # Run package-dependent commands before solving if any were defined
    [[ -e ${MYCUTER}/bin/${PACK}_pre ]] && . ${MYCUTER}/bin/${PACK}_pre
}

###############################################################################

# Run post script

Run_Post() {

    # Run post-solving package-dependent commands if any were defined
    [[ -e ${MYCUTER}/bin/${PACK}_post ]] && . ${MYCUTER}/bin/${PACK}_post
}

###############################################################################

# Clean up after run

Clean_Up() {

    # Tidy up the current directory, deleting all junk.
    if (( KEEP == 0 )); then
        ${RM} ${EXEC}/${PACK}min
        ${RM} ELFUN.o GROUP.o RANGE.o EXTER.o
    fi
}

###############################################################################

# Perform requested operations

thisprog=`basename $0`
WorkingDir=${PWD}

# Check that environment is correctly defined
envcheck
[[ $? != 0 ]] && exit $?

# Initialize default settings
Initialize_Settings

# Source commands
. ${MYCUTER}/${PRECISION}/config/cmds

# Parse command-line options
Parse_Arguments "$@"

if (( show_config )); then
    printf '========================\n'
    printf "  F77 compile:   ${COMPILE}\n"
    printf "  F77 link:      ${LOAD}\n"
    printf "  F95 compile:   ${COMPILE9095}\n"
    printf "  F95 link:      ${LOAD9095}\n"
    printf "  General flags: ${FFLAGS}\n"
    printf "  Special flags: ${SPECIALFLAGS}\n"
    printf "  Debug flags:   ${FDEBUGFLAGS}\n"
    printf "  Profile flags: ${FPROFILEFLAGS}\n"
    printf '\n'
    printf "  C compile:     ${CCOMPILE}\n"
    printf "  C link:        ${CLOAD}\n"
    printf "  General flags: ${CFLAGS}\n"
    printf "  Libraries:     ${SPECIALLIBS}\n"
    printf '\n'
    printf "  Mex compile:   ${MEXFORTRAN}\n"
    printf "  Mex flags:     ${MEXFFLAGS}\n"
    printf '========================\n'
    printf '\n'
fi

# Ensure that a package was specified
if (( package_set == 0 )); then
    echo "Please specify a package to run using -p"
    echo "Try $0 --help for more information"
    exit 1
fi

# Decode problem if required
(( problem_set )) && Decode_Problem ${PROBLEM}

# Source package definitions
. ${MYCUTER}/bin/${PKG}
[[ $? != 0 ]] && exit $?

# Main driver for the package
DRIVER=${PACK}ma.o

# Check that the required precision is available
PACK_PRECISION=( "${PACK_PRECISION}" )
let precision_ok=0
for prec in ${PACK_PRECISION}
do
    if [[ $prec == ${PRECISION} ]]; then
        let precision_ok=1
        break
    fi
done
if (( precision_ok == 0 )); then
    echo "Package ${PACKAGE} is not available in ${PRECISION} precision"
    exit 6
fi

Run_Pre

# Ensure CUTEr library is present
LIBDIR=${MYCUTER}/${PRECISION}/lib
if [[ ! -e ${LIBDIR}/libcuter.a  ]]; then
    cd ${LIBDIR}
    ${MAKE} -s libcuter.a
    cd ${WorkingDir}
fi

# Make sure package is up to date
cd ${MYCUTER}/${PRECISION}/bin
[[ ${PKG} != "mx" ]] && ${MAKE} -s ${DRIVER}
cd ${WorkingDir}

# If needed, find the correct BLAS and package spec file
(( blas_set == 0 )) && BLAS=${MYCUTER}/${PRECISION}/bin/linpac.o
if [[ ${SPECS} != "" ]]; then
  if [[ ! -e ${SPECS} ]]; then
     if [[ -e ${MYCUTER}/${PRECISION}/cuter-specs/${SPECS} ]]; then
       ${LN} ${MYCUTER}/${PRECISION}/cuter-specs/${SPECS} ${SPECS}
     elif [[ -e ${CUTER}/common/src/pkg/${PACKAGE}/${SPECS} ]]; then
       ${LN} ${CUTER}/common/src/pkg/${PACKAGE}/${SPECS} ${SPECS}
     else
        printf "\nCannot find spec file ${SPECS}---skipping\n"
     fi
  fi
fi

# Ensure that the current test problem has been compiled.
(( OUTPUT )) && printf '\nCompiling current test problem if necessary ...\n'
(( RECOMPILE )) && ${RM} ELFUN.o GROUP.o RANGE.o EXTER.o
for i in  ELFUN GROUP RANGE
do
    if [[ ! -e ${i}.o  ]]; then
        command="${COMPILE} ${FFLAGS} ${i}.f"
        (( show_commands )) && echo $command
        $command
        [[ $? != 0 ]] && exit $?
    fi
done

EXTER=""
[[ -e EXTER.f && -z EXTER.f ]] && ${RM} EXTER.f
if [[ -e EXTER.f ]]; then
    command="${COMPILE} ${FFLAGS} EXTER.f"
    (( show_commands )) && echo $command
    $command
    [[ $? != 0 ]] && exit $?
    [[ -e EXTER.o && -z EXTER.o ]] && ${RM} EXTER.o || EXTER="EXTER.o"
fi

# The package-dependent object files are in ${MYCUTER}/${PRECISION}/bin
PACKOBJS=( ${PACKOBJS} )
nobjs=${#PACKOBJS[@]}
cd ${MYCUTER}/${PRECISION}/bin
for (( i = 0; i < nobjs ; i++ ))
do
    ${MAKE} -s ${PACKOBJS[$i]}
    PACKOBJS[$i]=${MYCUTER}/${PRECISION}/bin/${PACKOBJS[$i]}
done
cd ${WorkingDir}

# If user requested Matlab, the procedure is special
if [[ ${PKG} == "mx" ]]; then
    (( OUTPUT )) && printf '\nBuilding MEX file ...\n'
    outputName=mcuter
    command="${MEXFORTRAN} -I${CUTER}/common/include -output ${outputName} ${CUTER}/common/src/tools/mcuter.c ELFUN.o GROUP.o RANGE.o ${EXTER} -L${LIBDIR} -lcuter ${BLAS} ${LAPACK} ${PACKLIBS}"
    (( show_commands )) && echo "$command"
    $command
    Run_Post
    Clean_Up
    exit 0
else
    # Link all the PACK and tools files together.
    (( OUTPUT )) && printf '\nLinking all the object files together ...\n'
    command="${LOAD} ${FFLAGS} -o ${PACK}min ELFUN.o GROUP.o RANGE.o ${EXTER} ${MYCUTER}/${PRECISION}/bin/${DRIVER} ${PACKOBJS[@]} ${ALT_LIB_PATH[@]} -L${LIBDIR} ${PACKLIBS} ${SPECIALLIBS} -lcuter ${BLAS} ${LAPACK} ${XTRALIBS[@]}"
    (( show_commands )) && echo "$command"
    $command
fi

[[ ${PWD} != ${EXEC} ]] && ${MV} ${PACK}min ${EXEC}/${PACK}min

#  run PACK on the current test problem unless -debug is set.

if (( DEBUG )); then
    echo '  debug enabled, load module is in '${EXEC}'/'
else
    Run_Package
    Run_Post
    Clean_Up
fi
exit 0
