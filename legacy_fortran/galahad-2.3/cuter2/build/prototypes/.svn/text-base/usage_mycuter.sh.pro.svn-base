# Environment check

if [[ ${CUTER+set} != 'set' ]]; then
    echo ' CUTER is not set.'
    echo ' Set it to the appropriate value and re-run.'
    echo ' Aborting.'
    exit 1
fi

if [[ -e ./install_mycuter ]]; then
    if [[ ${MYCUTER+set} == 'set' && "$MYCUTER" == "$PWD" ]]; then
        export MYCUTER=`echo ${MYCUTER} | ${SED} 's"/tmp_mnt""'`
    fi
else
    echo ' Launch install_mycuter from its home directory.'
    echo ' Aborting.'
    exit 2
fi

# End of environment check

# Calling sequence for install_mycuter

thisprog=`basename $0`
let doubleprec=1
x_umake_flags=( " " )
let debug=0
let profile=0

while [[ $# > 0 ]]; do
    case "$1" in
        -single)  let doubleprec=0
                  ;;
        -double)  let doubleprec=1
                  ;;
        -debug)   (( debug == 0 )) && x_umake_flags=( "${x_umake_flags[@]} -Ddebug" )
                  let debug=1
                  ;;
        -profile) (( profile == 0 )) && x_umake_flags=( "${x_umake_flags[@]} -Dprofile" )
                  let profile=1
                  ;;
        *)        echo "Use: ${thisprog} [arg]"
                  echo "where arg is either"
                  echo "   -single   install the single precision decoder"
                  echo "   -double   install the double precision decoder"
                  echo "   -debug    produce debugging information"
                  echo "   -profile  produce profiling information"
                  exit 0
                  ;;
    esac
    shift
done

if [[ ! -d ${MYCUTER}/double && "${doubleprec}" == "true" ]]; then
    echo "Please use  install_cuter  to create double precision version"
    exit 2
fi
if [[ ! -d ${MYCUTER}/single && "${singleprec}" == "true" ]]; then
    echo "Please use  install_cuter  to create single precision version"
    exit 3
fi

if [[ -d ${MYCUTER}/double && "${doubleprec}" == "false" ]]; then
    if [[ -d ${MYCUTER}/single && "${singleprec}" == "false" ]]; then
        echo "Please specify which precision to build."
        echo " Use $0 -h for more information."
        exit 1
    fi
fi

(( doubleprec )) && x_umake_flags=( "${x_umake_flags[@]} -DDoublePrecision" ) || \
                    x_umake_flags=( "${x_umake_flags[@]} -DSinglePrecision" )

# End of arg check

if [[ ${LOQODIR+set} == set ]]; then
  echo "You appear to have LOQO installed in $LOQODIR"
  if [[  -f ${LOQODIR}/loqo.h && -f ${LOQODIR}/myalloc.h && -f ${LOQODIR}/libloqo.a  ]]; then
    x_umake_flags=( "${x_umake_flags[@]} -Dhasloqo")
  else
    echo 
    echo "At least one of the files loqo.h, myalloc.h or libloqo.a"
    echo "seems to be missing in ${LOQODIR}... Please re-run install_mycuter"
    echo "when the problem is fixed."
    echo
  fi
else
  echo "You do not seem to have LOQO installed"
fi
if [[ ${KNITRODIR+set} == set ]]; then
  echo "You appear to have KNITRO installed in $KNITRODIR"
  if [[ -f ${KNITRODIR}/include/knitro.h && -d ${KNITRODIR}/lib  ]]; then
    x_umake_flags=("${x_umake_flags[@]} -Dhasknitro")
  else
    echo 
    echo "At least one of the files knitro.h or libknitro*"
    echo "seems to be missing in ${KNITRODIR}... Please re-run install_mycuter"
    echo "when the problem is fixed."
    echo
  fi
else
  echo "You do not seem to have KNITRO installed"
fi
if [[ ${TAO_DIR+set} == set ]]; then
  echo "You appear to have TAO installed in ${TAO_DIR}"
  x_umake_flags=("${x_umake_flags[@]} -Dhastao")
else
  echo "You do not seem to have TAO installed"
fi
