#!/bin/csh -f
#  ( Last modified on Thu Jan 10 13:38:02 MET 2002 )
#
# envcheck: check CUTEr environment variables
#

@ abort = 0

if( ! $?CUTER ) then
    echo ' CUTER is not set.'
    echo ' It should point to the directory where you installed CUTEr.'
    echo ' Set it to the appropriate value and re-run.'
    echo ''
    @ abort = 1
endif

if( ! $?MYCUTER ) then
    echo ' MYCUTER is not set.'
    echo ' It should point to your working instance of CUTEr.'
    echo ' Set it to the appropriate value and re-run.'
    echo ''
    @ abort = 1
endif

if( ! $?MASTSIF ) then
    echo ' MASTSIF is not set.'
    echo ' It should point to your main SIF repository.'
    echo ' Set it to the appropriate value and re-run.'
    echo ''
    @ abort = 1
endif

if( $abort ) then
    echo ' Aborting.'
    exit 7
endif

exit 0
