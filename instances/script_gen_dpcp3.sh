#!/bin/bash

declare N="200"
declare NA="30" # 0.2*N
declare -a P=("0.25" "0.5" "0.75")
declare -a NB=("2" "4" "9")
declare num="5"

SRC="gen_random_dpcp1.py"
OUT="dpcp/er-infeas"

mkdir -p dpcp
mkdir -p dpcp/er-infeas

for i in "${!P[@]}"
do
    p="${P[$i]}"
	nb="${NB[$i]}"
	echo "n: $N, p: $p, nA: $NA, nB: $nb, i: $i"
	python3 $SRC $N $p $NA $nb $num $OUT
	for i in $(seq 0 $((num - 1)))
	do
		echo "r_n${N}_p${p}_nA${NA}_nB${nb}_i${i}.dpcp" >> "$OUT/instances.txt"
	done
done
