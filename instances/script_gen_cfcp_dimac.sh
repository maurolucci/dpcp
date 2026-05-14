#!/bin/bash

SRC="dimacs2dpcp.py"
OUT1="cfcp/open"
OUT2="cfcp/closed"

mkdir -p cfcp
mkdir -p cfcp/open
mkdir -p cfcp/closed

for file in dimac100/*.col:
    filename=$(basename -- "$file")
    name="${filename%.*}"
    python3 $SRC $name $OUT1
    python3 $SRC $name $OUT2 --cerrada
    echo "$name.dpcp" >> "$OUT1/instances.txt"
    echo "$name.dpcp" >> "$OUT2/instances.txt"
python3