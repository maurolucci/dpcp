#!/bin/bash

declare -a N=("18" "12" "9")
declare -a R=("20" "35" "50")
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
		python3 $CONV "$OUT/freq_c${n}_r${r}_i${i}"
		echo "freq_c${n}_r${r}_i${i}.dpcp" >> "$OUT/instances.txt"
	done
done
