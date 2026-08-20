#!/bin/bash

declare N="150"
declare -a P=("0.25" "0.5" "0.75")
declare -a NA=("0.2" "0.3")
declare -a NB=("0.2" "0.3")
declare num="3"

SRC="gen_random_dpcp1.py"
OUT="dpcp/er-1"

mkdir -p dpcp
mkdir -p dpcp/er-1

for p in "${P[@]}"
do
    for na in "${NA[@]}"
	do
		for nb in "${NB[@]}"
	    do
		    naa=$(echo "scale=0; ($na * $N + 0.5)/1" | bc)
		    nbb=$(echo "scale=0; ($nb * $N + 0.5)/1" | bc)
		    echo "n: $N, p: $p, nA: $naa, nB: $nbb, i: $i"
		    python3 $SRC $N $p $naa $nbb $num $OUT
		    for i in $(seq 0 $((num - 1)))
		    do
		        echo "r_N${N}_p${p}_n${naa}_m${nbb}_i${i}.dpcp" >> "$OUT/instances.txt"
		    done
        done
    done
done
