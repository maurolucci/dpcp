#!/bin/bash

declare -a N=("90" "100" "110")
declare -a P=("0.25" "0.5" "0.75")
declare -a NA=("0.1" "0.2")
declare -a NB=("0.1" "0.2")
declare num="5"

SRC="gen_random_dpcp1.py"
OUT="dpcp/er2"

mkdir -p dpcp
mkdir -p dpcp/er2

for i in "${!N[@]}"
do
    n="${N[$i]}"
    p="${P[$i]}"
    for na in "${NA[@]}"
	do
	    if [[ $p == "0.25" && $na == "0.1" ]]; then
			continue
		fi
		for nb in "${NB[@]}"
	    do
		    naa=$(echo "scale=0; ($na * $n + 0.5)/1" | bc)
		    nbb=$(echo "scale=0; ($nb * $n + 0.5)/1" | bc)
		    echo "n: $n, p: $p, nA: $naa, nB: $nbb, i: $i"
		    python3 $SRC $n $p $naa $nbb $num $OUT
		    for i in $(seq 0 $((num - 1)))
		    do
		        echo "r_n${n}_p${p}_nA${naa}_nB${nbb}_i${i}.dpcp" >> "$OUT/instances.txt"
		    done
        done
    done
done
