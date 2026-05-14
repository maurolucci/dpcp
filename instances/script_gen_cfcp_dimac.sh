#!/bin/bash

SRC="dimacs2dpcp.py"
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
    python3 $SRC $name $OUT1
    python3 $SRC $name $OUT2 --cerrada
    echo "$basename.dpcp" >> "$OUT1/instances.txt"
    echo "$basename.dpcp" >> "$OUT2/instances.txt"
done
