#!/bin/csh
#  ( Last modified on 23 Dec 2000 at 17:29:56 )
#
#{version}
#
#  classify all SIF files in MASTSIF or in the directory given in an argument
#
#
#  I. Bongartz, A.R. Conn, Nick Gould and Ph. Toint, for CGT Productions.
#

#{cmds}

if ( $#argv == 1 ) then
  cd $1
else
  cd $MASTSIF
endif
if ( -e CLASSF.DB ) $RM CLASSF.DB
foreach i ( *.SIF )
  set input=`echo $i|sed "s/.SIF//"`
  echo "y"       > CLSF.DAT
  echo "$input" >> CLSF.DAT
  $MYCUTER/bin/clsf
  if( -e CLASSF.UDB) $MV CLASSF.UDB CLASSF.DB
end
if( -e CLASSF.DB) $CAT CLASSF.DB
$RM CLSF.DAT
