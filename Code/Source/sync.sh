#!/bin/bash

if [[ $# -lt 1 ]]
then
  echo "Usage: $0 <tgt-dir>"
  exit 0
fi

TGT_DIR=$1
FILES=`find . -name '*.[hc]pp' -exec diff -q {} $TGT_DIR/{} \; 2>&1 | cut -d ' ' -f 2 | grep -o '\./.*pp'`
for f in $FILES
do
  echo "Copying $f to $TGT_DIR/$f"
  cp $f "$TGT_DIR/$f"
done
