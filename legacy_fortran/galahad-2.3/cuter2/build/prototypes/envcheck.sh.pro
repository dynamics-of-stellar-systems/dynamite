#!/bin/bash -f
#  ( Last modified on Thu Jan 10 13:38:02 MET 2002 )
#
# envcheck: check CUTEr environment variables
#

let abort=0

if [[ ${CUTER+set} != 'set' ]]; then
    echo ' CUTER is not set.'
    echo ' It should point to the directory where you installed CUTEr.'
    echo ' Set it to the appropriate value and re-run.'
    echo ''
    let abort=1
fi

if [[ ${MYCUTER+set} != 'set' ]]; then
    echo ' MYCUTER is not set.'
    echo ' It should point to your working instance of CUTEr.'
    echo ' Set it to the appropriate value and re-run.'
    echo ''
    let abort=1
fi

if [[ ${MASTSIF+set} != 'set' ]]; then
    echo ' MASTSIF is not set.'
    echo ' It should point to your main SIF repository.'
    echo ' Set it to the appropriate value and re-run.'
    echo ''
    let abort=1
fi

if (( abort > 0 )); then
    echo ' Aborting.'
    exit 7
fi

exit 0
