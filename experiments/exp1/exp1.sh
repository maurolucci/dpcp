#!/bin/bash

# exp1: Evaluate parameter ALPHA in the semigreedy heuristic

declare -a HEURS=("semigreedy2s")
declare -a VARIANTS=("2" "3")
declare -a ALPHA=("0.1" "0.15" "0.2" "0.25" "0.3" "0.35" "0.4")

declare INPUT="../../instances/dpcp/random"
declare INSTANCES="$INPUT/instances.txt"
declare BIN="../../dpcp"
declare OUT="out/"

echo "Running experiment #1"

# First, create output directories
for h in "${HEURS[@]}"
do
    mkdir -p "$OUT/$h"
    if [[ "$h" == "greedy2s" ]]; then
        for v in "${VARIANTS[@]}"
        do
            mkdir -p "$OUT/$h/v$v"
        done
    elif [[ "$h" == "semigreedy2s" ]]; then
        for v in "${VARIANTS[@]}"
        do
            mkdir -p "$OUT/$h/v$v"
            for a in "${ALPHA[@]}"
            do
                mkdir -p "$OUT/$h/v$v/a$a"
            done
        done
    fi
done

# Second, run experiments
while IFS= read -r LINE
do
    echo "Processing instance: $LINE"
    for h in "${HEURS[@]}"
    do
        if [[ "$h" == "greedy1s" ]]; then
            echo "Solving with heuristic: $h"
            time $BIN -s heur -f "$INPUT/$LINE" -o "$OUT/$h/" --heur-root 1 --preproc-off
        elif [[ "$h" == "greedy2s" ]]; then
            for v in "${VARIANTS[@]}"
            do
                echo "Solving with heuristic: $h, variant: $v"
                time $BIN -s heur -f "$INPUT/$LINE" -o "$OUT/$h/v$v/" --heur-root 2 --heur-2step-variant $v --preproc-off
            done
        elif [[ "$h" == "semigreedy2s" ]]; then
            for v in "${VARIANTS[@]}"
            do
                for a in "${ALPHA[@]}"
                do
                    echo "Solving with heuristic: $h, variant: $v, alpha: $a"
                    time $BIN -s heur -f "$INPUT/$LINE" -o "$OUT/$h/v$v/a$a/" --heur-root 3 --heur-2step-variant $v --heur-semigreedy-alpha $a --heur-semigreedy-iter 100 --preproc-off
                done
            done
        fi
    done
done < "$INSTANCES"
