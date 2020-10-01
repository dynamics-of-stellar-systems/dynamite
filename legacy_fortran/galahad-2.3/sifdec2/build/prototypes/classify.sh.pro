#!/bin/bash
#  ( Last modified on 23 Dec 2000 at 17:29:56 )
#
#{version}
#
#  add to the CLASSF.DB file the classification of the problem whose name is
#  given as an argument
#  If two arguments are given the first indicates the name of the directory
#  where this problem resides, otherwise the default is that pointed to by
#  the shell variable $MASTSIF.

#  A.R. Conn and Ph.L. Toint, October 1992, for CGT Productions.
#  modified by I. Bongartz, July 1994.
#

#{cmds}

if [[  $# == 1  ]]; then
  cd $MASTSIF
  echo "n"  > CLSF.DAT
  echo "$1">> CLSF.DAT
  echo ' '
  echo '   Your current classification file is : CLASSF.DB '
  echo ' '
  $MYSIFDEC/bin/clsf
  [[  -e CLASSF.UDB  ]] && $MV CLASSF.UDB CLASSF.DB
elif [[  $# == 2  ]]; then
  cd $1
  echo "n"   > CLSF.DAT
  echo "$2" >> CLSF.DAT
  echo ' '
  echo '   Your current classification file is : CLASSF.DB '
  echo ' '
  $MYSIFDEC/bin/clsf
  [[  -e CLASSF.UDB  ]] && $MV CLASSF.UDB CLASSF.DB
else
  echo ' Usage "classify name" or "classify directory name"'
  echo ''
  echo ' where name is the problem name without the .SIF appended,'
  echo '       directory is the name of the directory containing this problem'
fi
$RM CLSF.DAT
