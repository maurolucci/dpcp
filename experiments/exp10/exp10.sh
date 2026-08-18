#!/bin/bash

# exp10: Evaluation of B&P in geometric CFCP instances

declare TIME_LIMIT="1200"
declare TREE_SEARCH="1"
declare VERBOSE="2"

# Initial heuristic parameters:
declare HEURISTIC_INITIAL="3"
declare HEURISTIC_INITIAL_ALPHA="0.15"
declare HEURISTIC_INITIAL_REPETITIONS="40000"
declare HEURISTIC_INITIAL_MAXTIME="300"

# Initialization parameters:
declare HEUR_NODES="2"
declare VARIANT="3"

# Feasibility parameters:
declare FEAS_ROOT="0"
declare FEAS_NODES="0"

# Columns inheritance parameters:
declare INHERIT_COLS="3"

# Pricing parameters:
declare PRICING_METHOD="4"
declare ALPHA_PRI="0.2"
declare GREEDY_MAX_COLS="10"
declare MAX_COLS_PER_ITER="20"
declare PRICING_EXACT_TIME="7200"

# Branching parameters:
declare BRANCHING_VARIABLE="1"

declare INPUT="../../instances/cfcp/circle"
declare INSTANCES="$INPUT/instances.txt"
declare BIN="../../dpcp"
declare OUT="out/"

echo "Running experiment #10"

# First, create output directories
mkdir -p "$OUT"

# Second, run experiments
while IFS= read -r LINE
do
    echo "Processing instance: $LINE"
    time $BIN -s byp -f "$INPUT/$LINE" -o "$OUT/" -t $TIME_LIMIT --tree-search $TREE_SEARCH --verbose $VERBOSE \
    --heur-initial $HEURISTIC_INITIAL --heur-nodes $HEUR_NODES --heur-2step-variant $VARIANT \
    --heur-semigreedy-alpha $HEURISTIC_INITIAL_ALPHA --heur-semigreedy-iter $HEURISTIC_INITIAL_REPETITIONS \
    --heur-semigreedy-time $HEURISTIC_INITIAL_MAXTIME \
    --feas-root $FEAS_ROOT --feas-nodes $FEAS_NODES --inherit-cols $INHERIT_COLS \
    --pricing-method $PRICING_METHOD --pricing-greedy-alpha $ALPHA_PRI \
    --pricing-greedy-max-cols $GREEDY_MAX_COLS --pricing-max-cols-per-iter $MAX_COLS_PER_ITER \
    --pricing-exact-time $PRICING_EXACT_TIME --branching-variable $BRANCHING_VARIABLE

    time $BIN -s compact -f "$INPUT/$LINE" -o "$OUT/" -t $TIME_LIMIT --verbose $VERBOSE \
    --heur-initial $HEURISTIC_INITIAL --heur-2step-variant $VARIANT \
    --heur-semigreedy-alpha $HEURISTIC_INITIAL_ALPHA --heur-semigreedy-iter $HEURISTIC_INITIAL_REPETITIONS \
    --heur-semigreedy-time $HEURISTIC_INITIAL_MAXTIME

    time $BIN -s gurobi -f "$INPUT/$LINE" -o "$OUT/" -t $TIME_LIMIT --verbose $VERBOSE \
    --heur-initial $HEURISTIC_INITIAL --heur-2step-variant $VARIANT \
    --heur-semigreedy-alpha $HEURISTIC_INITIAL_ALPHA --heur-semigreedy-iter $HEURISTIC_INITIAL_REPETITIONS \
    --heur-semigreedy-time $HEURISTIC_INITIAL_MAXTIME

done < "$INSTANCES"
