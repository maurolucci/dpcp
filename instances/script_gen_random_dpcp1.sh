#!/bin/bash

declare -a N=("100" "150")
declare -a P=("0.25" "0.5" "0.75")
declare -a NA=("0.1" "0.2")
declare -a NB=("0.1" "0.2")
declare num="3"

SRC="gen_random_dpcp1.py"
OUT="dpcp/random/"

for n in "${N[@]}"
do
    for p in "${P[@]}"
    do
        for na in "${NA[@]}"
        do
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
done