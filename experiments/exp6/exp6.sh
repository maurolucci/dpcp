#!/bin/bash

# exp6: Evaluation of pricing

# DPCP parameters:
declare TIME="900"

# Initial heuristic parameters:
declare HEURISTIC_INITIAL="0"

# DPCP heuristic parameters:
declare HEUR_NODES="2"
declare VARIANT="3"

# Pricing parameters:
declare ALPHA_PRI="0.2"
declare GREEDY_MAX_COLS="100"
declare MAX_COLS_PER_ITER="20"
declare TIME_PRI="900"

declare INPUT="../../instances/dpcp/er-1"
declare INSTANCES="$INPUT/instances.txt"
declare BIN="../../dpcp"
declare OUT="out/"

echo "Running experiment #6"

# First, create output directories
for p in "0" "1" "2" "3" "4" "5" "6"
do
    mkdir -p "$OUT/p$p"
done

# Second, run experiments
while IFS= read -r LINE
do
    echo "Processing instance: $LINE"
    # Baseline #0: exact, no initial heuristic
    time $BIN -s byp -f "$INPUT/$LINE" -o "$OUT/p0/" --relax --heur-initial $HEURISTIC_INITIAL --heur-root 0 --feas-root 0 --feas-nodes 0 --inherit-cols 0 --pricing-method 0 --pricing-exact-time $TIME_PRI --time $TIME
    # Pricing #0: exact
    time $BIN -s byp -f "$INPUT/$LINE" -o "$OUT/p1/" --relax --heur-initial $HEURISTIC_INITIAL --heur-nodes $HEUR_NODES --heur-2step-variant $VARIANT --feas-root 0 --feas-nodes 0 --inherit-cols 0 --pricing-method 0 --pricing-exact-time $TIME_PRI --time $TIME
    # Pricing #1: greedy + exact
    time $BIN -s byp -f "$INPUT/$LINE" -o "$OUT/p2/" --relax --heur-initial $HEURISTIC_INITIAL --heur-nodes $HEUR_NODES --heur-2step-variant $VARIANT --feas-root 0 --feas-nodes 0 --inherit-cols 0 --pricing-method 1 --pricing-greedy-alpha $ALPHA_PRI --pricing-greedy-max-cols $GREEDY_MAX_COLS --pricing-max-cols-per-iter $MAX_COLS_PER_ITER --pricing-exact-time $TIME_PRI --time $TIME
    # Pricing #2: greedy + P,Q-MWSSP + exact
    time $BIN -s byp -f "$INPUT/$LINE" -o "$OUT/p3/" --relax --heur-initial $HEURISTIC_INITIAL --heur-nodes $HEUR_NODES --heur-2step-variant $VARIANT --feas-root 0 --feas-nodes 0 --inherit-cols 0 --pricing-method 2 --pricing-greedy-alpha $ALPHA_PRI --pricing-greedy-max-cols $GREEDY_MAX_COLS --pricing-max-cols-per-iter $MAX_COLS_PER_ITER --pricing-exact-time $TIME_PRI --time $TIME
    # Pricing #3: greedy + P-MWSSP + exact
    time $BIN -s byp -f "$INPUT/$LINE" -o "$OUT/p4/" --relax --heur-initial $HEURISTIC_INITIAL --heur-nodes $HEUR_NODES --heur-2step-variant $VARIANT --feas-root 0 --feas-nodes 0 --inherit-cols 0 --pricing-method 3 --pricing-greedy-alpha $ALPHA_PRI --pricing-greedy-max-cols $GREEDY_MAX_COLS --pricing-max-cols-per-iter $MAX_COLS_PER_ITER --pricing-exact-time $TIME_PRI --time $TIME
    # Pricing #4: greedy + P,Q-MWSSP + P-MWSSP + exact
    time $BIN -s byp -f "$INPUT/$LINE" -o "$OUT/p5/" --relax --heur-initial $HEURISTIC_INITIAL --heur-nodes $HEUR_NODES --heur-2step-variant $VARIANT --feas-root 0 --feas-nodes 0 --inherit-cols 0 --pricing-method 4 --pricing-greedy-alpha $ALPHA_PRI --pricing-greedy-max-cols $GREEDY_MAX_COLS --pricing-max-cols-per-iter $MAX_COLS_PER_ITER --pricing-exact-time $TIME_PRI --time $TIME
    # Pricing #5: greedy + P-MWSSP + P,Q-MWSSP + exact
    time $BIN -s byp -f "$INPUT/$LINE" -o "$OUT/p6/" --relax --heur-initial $HEURISTIC_INITIAL --heur-nodes $HEUR_NODES --heur-2step-variant $VARIANT --feas-root 0 --feas-nodes 0 --inherit-cols 0 --pricing-method 5 --pricing-greedy-alpha $ALPHA_PRI --pricing-greedy-max-cols $GREEDY_MAX_COLS --pricing-max-cols-per-iter $MAX_COLS_PER_ITER--pricing-exact-time $TIME_PRI --time $TIME
done < "$INSTANCES"
