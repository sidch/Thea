#!/bin/bash

cd G3D

REF_DIR="/Volumes/BOOTCAMP/Users/sid/Documents/GitHub/Fuse_OldWin/Code/Dependencies/Thea/Code/Dependencies/G3D/G3D/"
FILES=`find G3D.lib \( -name '*pp' -or -name '*.h' \) -exec diff -U3 {} $REF_DIR/{} \; | grep '^---' | cut -d ' ' -f 2 | cut -f 1`

echo "$FILES"

for f in $FILES
do
  diff -U3 $f $REF_DIR/$f > ../Patches/G3D-8.00/$f.patch
done

cd ..
