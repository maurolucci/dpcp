#!/bin/bash

# exp1: Evaluate DPCP heuristics

declare -a HEURS=("greedy1s" "greedy2s" "semigreedy2s")
declare -a VARIANTS=("2" "3")
declare -a REPETITIONS=("10" "100" "1000")

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
            for r in "${REPETITIONS[@]}"
            do
                mkdir -p "$OUT/$h/v$v/r$r"
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
                for r in "${REPETITIONS[@]}"
                do
                    echo "Solving with heuristic: $h, variant: $v, repetition: $r"
                    time $BIN -s heur -f "$INPUT/$LINE" -o "$OUT/$h/v$v/r$r/" --heur-root 3 --heur-2step-variant $v --heur-root-iter $r --preproc-off
                done
            done
        fi
    done
done < "$INSTANCES"