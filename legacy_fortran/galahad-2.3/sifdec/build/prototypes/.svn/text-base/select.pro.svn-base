#!/bin/csh
#  ( Last modified on 23 Dec 2000 at 17:29:56 )
#
#{version}

#  Probes the classification file using slct.  By default,
#  the classification file is assumed to be $MYSIFDEC/CLASSF.DB,
#  but the slct program allows the user to give a full path name
#  for the classification file, thereby allowing both the filename
#  and the directory to be changed.
#
#
#  A.R. Conn, October 1992, for CGT Productions.
#  modified by I. Bongartz, May 1994.
#

#{cmds}

echo $MASTSIF > SLCT.DAT
$MYSIFDEC/bin/slct
$RM SLCT.DAT
