#!/bin/bash -f

# Remove the archive file FILE from the archive ARC if FILE is present

# syntax rmarfile ar-command grep-command arc-name file-name ,

# where ar-command   is the appropriate "ar" command,
#       grep-command is the appropriate "grep" command,
#       arc-name     is the name of the archive, and
#       file-name    is the name of the file to be removed

AR=$1
GREP=$2
ARC=$3
FILE=$4
#[[  `$AR -t $ARC | $GREP $FILE` == $FILE  ]] && echo "removing $FILE" || echo "no $FILE to remove"
if [[ -e $ARC ]] ; then
  [[  `$AR -t $ARC | $GREP -w $FILE` == $FILE  ]] && $AR -d $ARC $FILE
fi
#echo "done"

exit 0
