#!/bin/bash -f

# Remove the object file FILE from the directory FILEDIR if FILE is present

# syntax rmobfile rm-command file-dir file-name ,

# where rm-command   is the appropriate "rm" command,
#       file-dir     is the name of the directory from which the
#                    file to be removed
#       file-name    is the name of the file to be removed

f=$#
let d=f-1
let a=d-1

RM="${*:1:a}"
FILEDIR="${!d}"
FILE="${!f}"

[[ -e ${FILEDIR}/${FILE} ]] && $RM ${FILEDIR}/${FILE}

exit 0
