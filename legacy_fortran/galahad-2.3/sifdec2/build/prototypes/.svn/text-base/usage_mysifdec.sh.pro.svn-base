
# Environment check

if [[ ${SIFDEC+set} != 'set'  ]]; then
    echo ' SIFDEC is not set.'
    echo ' Set it to the appropriate value and re-run.'
    echo ' Aborting.'
    exit 1
fi

let mysifdecIsCWD=0

if [[  -e ./install_mysifdec  ]]; then
    if [[ ${MYSIFDEC+set} != 'set'  ]]; then
	    let mysifdecIsCWD=1
	    export MYSIFDEC=`dirs -l`
	    export MYSIFDEC=`echo $MYSIFDEC | $SED 's"/tmp_mnt""'`
    fi
else
    echo ' Launch install_mysifdec from its home directory.'
    echo ' Aborting.'
    exit 2
fi

# End of environment check

# Calling-sequence for install_mysifdec

singleprec="false"
doubleprec="false"
x_umake_flags=( " " )
let debug=0
let profile=0

let last=$#

let i=1
while (( i <= last ))
do
    opt=${!i}

    if [[ "$opt" == "-single" ]]; then
        singleprec="true"
        doubleprec="false"
    elif [[ "$opt" == "-double" ]]; then
        singleprec="false"
        doubleprec="true"
    elif [[ "$opt" == "-debug" ]]; then
        (( debug == 0 )) && x_umake_flags=( "${x_umake_flags[@]} -Ddebug" )
        let debug=1
    elif [[ "$opt" == "-profile" ]]; then
        (( profile == 0 )) && x_umake_flags=( "${x_umake_flags[@]} -Dprofile" )
        let profile=1
    else
        echo "Use: install_mysifdec [arg]"
        echo "where arg is either"
        echo "   -single   install the single precision decoder"
        echo "   -double   install the double precision decoder"
        echo "   -debug    produce debugging information"
        echo "   -profile  produce profiling information"
        exit 0
    fi

    (( i++ ))
done

if [[ ! -d ${MYSIFDEC}/double && "$doubleprec" == "true" ]]; then
    echo "Please use  install_sifdec  to create double precision version"
    exit 2
fi
if [[ ! -d ${MYSIFDEC}/single && "$singleprec" == "true" ]]; then
    echo "Please use  install_sifdec  to create single precision version"
    exit 3
fi

if [[ -d ${MYSIFDEC}/double && "$doubleprec" == "false" ]]; then
    if [[ -d ${MYSIFDEC}/single && "$singleprec" == "false" ]]; then
        echo "Please specifiy which precision to build."
        echo "Use  $0 -h  for more information."
        exit 1
    fi
fi
if [[ "$singleprec" == "true" ]]; then
    x_umake_flags=( "${x_umake_flags[@]} -DSinglePrecision" )
else
    x_umake_flags=( "${x_umake_flags[@]} -DDoublePrecision" )
fi

# End of arg check

