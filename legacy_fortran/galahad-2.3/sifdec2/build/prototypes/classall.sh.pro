#!/bin/bash
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

if [[  $# == 1  ]]; then
  cd $1
else
  cd $MASTSIF
fi
[[  -e CLASSF.DB  ]] && $RM CLASSF.DB
for i  in  *.SIF; do
  input=`echo $i|sed "s/.SIF//"`
  echo "y"       > CLSF.DAT
  echo "$input" >> CLSF.DAT
  $MYSIFDEC/bin/clsf
  [[  -e CLASSF.UDB ]] && $MV CLASSF.UDB CLASSF.DB
done
[[  -e CLASSF.DB ]] && $CAT CLASSF.DB
$RM CLSF.DAT
