#!/bin/bash

# exp3: Evaluate DPCP heuristics

declare -a HEURS=("greedy1s" "greedy2s" "semigreedy2s")
declare -a VARIANTS=("2" "3")
declare ALPHA="0.1"
declare REPETITIONS="500"

declare INPUT="../../instances/dpcp/random"
declare INSTANCES="$INPUT/instances.txt"
declare BIN="../../dpcp"
declare OUT="out/"

echo "Running experiment #3"

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
            	echo "Solving with heuristic: $h, variant: $v, repetition: $REPETITIONS"
            	time $BIN -s heur -f "$INPUT/$LINE" -o "$OUT/$h/v$v/" --heur-root 3 --heur-2step-variant $v --heur-semigreedy-iter $REPETITIONS --heur-semigreedy-alpha $ALPHA --preproc-off
            done
        fi
    done
done < "$INSTANCES"
