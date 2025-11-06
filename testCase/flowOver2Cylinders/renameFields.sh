#!/usr/bin/env bash
# Usage: ./rename_fields.sh OLD NEW [OLD NEW]...
# Example: ./rename_fields.sh U U1 p p1
#in case it is needed to rename fields. It is imporant to remember that paraview will not properly see the fields since internally they will have their original name. 
# The ROM codes instead, can be used correctly with the new field names. 
set -euo pipefail

# need pairs: OLD NEW ...
(( $# >= 2 && $# % 2 == 0 )) || { echo "Usage: $0 OLD NEW [OLD NEW]..."; exit 1; }

args=("$@")

for t in $(foamListTimes -noZero); do
  for ((i=0; i<${#args[@]}; i+=2)); do
    from="${args[i]}"
    to="${args[i+1]}"
    src="$t/$from"
    dst="$t/$to"

    if [ -e "$src" ]; then
      if [ -e "$dst" ]; then
        echo "skip: $src → $dst (target exists)"
      else
        mv -v -- "$src" "$dst"
      fi
    fi
  done
done
