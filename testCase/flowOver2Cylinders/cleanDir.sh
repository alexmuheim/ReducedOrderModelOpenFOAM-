#!/usr/bin/env bash
#make exectuable: chmod +x cleanDir.sh 
# Usage: ./clean_fields.sh FIELD [FIELD...]
#e.g ./cleanDir.sh uRom pRom U1 U2
# in case it is needed to remove files from multiple time folders

set -euo pipefail
shopt -s nullglob
(( $# )) || { echo "Usage: $0 FIELD [FIELD...]"; exit 1; }

for t in $(foamListTimes -noZero); do
  for pat in "$@"; do
    for path in "$t"/$pat; do
      rm -rf -- "$path"
      echo "removed: $path"
    done
  done
done
