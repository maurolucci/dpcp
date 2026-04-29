#!/bin/bash

declare -a N=("70" "35" "18")
declare -a R=("5" "10" "20")
declare num="5"

SRC="gen_circle_cfc.py"
CONV="cfc2dpcp.py"
OUT="cfcp/circle"

mkdir -p cfcp
mkdir -p cfcp/circle

for j in "${!N[@]}"
do
    n="${N[$j]}"
    r="${R[$j]}"
	python3 $SRC $n $r $num $OUT
	for i in $(seq 0 $((num - 1)))
	do
		python3 $CONV "$OUT/circle_n${n}_r${r}_i${i}" $OUT
		echo "circle_n${n}_r${r}_i${i}.dpcp" >> "$OUT/instances.txt"
	done
done
