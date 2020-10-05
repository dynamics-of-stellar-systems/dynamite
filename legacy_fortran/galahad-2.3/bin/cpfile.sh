#!/bin/bash -f

# Copy file FILE to DESTFILE if FILE is present

# syntax cpfile cp-command file-name dest-file

# where cp-command   is the appropriate "cp" command,
#       file-name    is the name of the file to be copied
#       dest-name    is the name of the file after the copy

let f=$#
(( f-- ))

CMD="$@"
FILE="${!f}"

[[ -e $FILE ]] && $CMD

exit 0
