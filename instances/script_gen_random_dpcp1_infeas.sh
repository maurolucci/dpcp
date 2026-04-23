#!/bin/bash

declare -a N=("130" "150" "170")
declare -a P=("0.25" "0.5" "0.75")
declare -a NA=("0.2")
declare -a NB=("1" "3" "8")
declare num="3"

SRC="gen_random_dpcp1.py"
OUT="dpcp/er-infeas"

mkdir -p dpcp
mkdir -p dpcp/er-infeas

for i in "${!N[@]}"
do
    n="${N[$i]}"
    p="${P[$i]}"
    for na in "${NA[@]}"
	do
		for nb in "${NB[@]}"
	    do
		    if [[ $p == "0.25" && $nb != "1" ]]; then
				continue
			fi
		    if [[ $p == "0.5" && $nb != "3" ]]; then
				continue
			fi
		    if [[ $p == "0.75" && $nb != "8" ]]; then
				continue
			fi
		    naa=$(echo "scale=0; ($na * $n + 0.5)/1" | bc)
		    echo "n: $n, p: $p, nA: $naa, nB: $nb, i: $i"
		    python3 $SRC $n $p $naa $nb $num $OUT
		    for i in $(seq 0 $((num - 1)))
		    do
		        echo "r_n${n}_p${p}_nA${naa}_nB${nb}_i${i}.dpcp" >> "$OUT/instances.txt"
		    done
        done
    done
done
