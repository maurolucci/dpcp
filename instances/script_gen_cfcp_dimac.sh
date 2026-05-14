#!/bin/bash

SRC1="dimacs2cfc.py"
SRC2="cfc2dpcp.py"
INPUT="dimac100/*"
OUT1="cfcp/open"
OUT2="cfcp/closed"

mkdir -p cfcp
mkdir -p cfcp/open
mkdir -p cfcp/closed

for file in $INPUT
do
    basename=$(basename -- "$file")
    basename=${basename%.*}
    name="${file%.*}"
    python3 $SRC1 $name $OUT1
    python3 $SRC1 $name $OUT2 --cerrada
    python3 $SRC2 "$OUT1/$basename"
    python3 $SRC2 "$OUT2/$basename"
    echo "$basename.dpcp" >> "$OUT1/instances.txt"
    echo "$basename.dpcp" >> "$OUT2/instances.txt"
done
