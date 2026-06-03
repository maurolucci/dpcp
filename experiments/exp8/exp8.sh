#!/bin/bash

# exp8: Infeasible instances

declare TIME_LIMIT="900"
declare TREE_SEARCH="1"
declare VERBOSE="2"

# DPCP heuristic parameters:
declare HEUR_ROOT="3"
declare HEUR_NODES="2"
declare VARIANT="2"
declare ALPHA_HEUR="0.1"
declare REPETITIONS_HEUR="40000"

# Feasibility parameters:
declare FEAS_ROOT="0"
declare FEAS_NODES="0"

# Columns inheritance parameters:
declare INHERIT_COLS="3"

# Pricing parameters:
declare PRICING_METHOD="6"
declare ALPHA_PRI="0.2"
declare GREEDY_MAX_COLS="100"
declare MAX_COLS_PER_ITER="10"
declare PRICING_EXACT_TIME="300"

# Branching parameters:
declare BRANCHING_VARIABLE="1"

declare INPUT="../../instances/infeasibles"
declare INSTANCES="$INPUT/instances.txt"
declare BIN="../../dpcp"
declare OUT="out/"

echo "Running experiment #8"

# First, create output directories
mkdir -p "$OUT"

# Second, run experiments
while IFS= read -r LINE
do
    echo "Processing instance: $LINE"

    # B&P  
    time $BIN -s byp -f "$INPUT/$LINE" -o "$OUT/" -t $TIME_LIMIT --tree-search $TREE_SEARCH --verbose $VERBOSE \
    --heur-root $HEUR_ROOT --heur-nodes $HEUR_NODES --heur-2step-variant $VARIANT \
    --heur-semigreedy-alpha $ALPHA_HEUR --heur-semigreedy-iter $REPETITIONS_HEUR \
    --feas-root $FEAS_ROOT --feas-nodes $FEAS_NODES --inherit-cols $INHERIT_COLS \
    --pricing-method $PRICING_METHOD --pricing-greedy-alpha $ALPHA_PRI \
    --pricing-greedy-max-cols $GREEDY_MAX_COLS --pricing-max-cols-per-iter $MAX_COLS_PER_ITER \
    --pricing-exact-time $PRICING_EXACT_TIME --branching-variable $BRANCHING_VARIABLE

    # ILP
    time $BIN -s feas-ilp -f "$INPUT/$LINE" -o "$OUT/" -t $TIME_LIMIT --verbose $VERBOSE

    # Maximum stable set
    timeout 15m $BIN -s feas-enum -f "$INPUT/$LINE" -o "$OUT/" --verbose $VERBOSE

done < "$INSTANCES"
