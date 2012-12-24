#!/bin/bash

if [[ $# -lt 1 ]]
then
  echo "Usage: $0 <tgt-dir>"
  exit 0
fi

TGT_DIR=$1

SRC_FILES=`find Code/Source -name '*.[hc]pp' -exec diff -q {} $TGT_DIR/{} \; 2>&1 | cut -d ' ' -f 2 | grep -o 'Code/Source/.*pp'`
BLD_FILES=`find Code/Build -name 'CMakeLists.txt' -exec diff -q {} $TGT_DIR/{} \; 2>&1 | cut -d ' ' -f 2 | grep -o 'Code/Build/.*CMakeLists.txt'`

FILES="$SRC_FILES $BLD_FILES"

for f in $FILES
do
  echo "Copying $f to $TGT_DIR/$f"
  mkdir -p `dirname "$TGT_DIR/$f"`
  cp $f "$TGT_DIR/$f"
done
