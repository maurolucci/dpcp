#!/bin/bash

# exp9: Infeasible instances

declare TIME_LIMIT="7200"
declare TREE_SEARCH="1"
declare VERBOSE="2"

# DPCP heuristic parameters:
declare HEUR_ROOT="0"
declare HEUR_NODES="0"

# Feasibility parameters:
declare FEAS_ROOT="0"
declare FEAS_NODES="0"

# Columns inheritance parameters:
declare INHERIT_COLS="3"

# Pricing parameters:
declare PRICING_METHOD="4"
declare ALPHA_PRI="0.2"
declare GREEDY_MAX_COLS="1000"
declare MAX_COLS_PER_ITER="10"
declare PRICING_EXACT_TIME="7200"

# Branching parameters:
declare BRANCHING_VARIABLE="1"

declare INPUT="../../instances/dpcp/er-infeas"
declare INSTANCES="$INPUT/instances.txt"
declare BIN="../../dpcp"
declare OUT="out/"

echo "Running experiment #9"

# First, create output directories
mkdir -p "$OUT"

# Second, run experiments
while IFS= read -r LINE
do
    echo "Processing instance: $LINE"

    # B&P  
    time $BIN -s byp -f "$INPUT/$LINE" -o "$OUT/" -t $TIME_LIMIT --tree-search $TREE_SEARCH --verbose $VERBOSE \
    --heur-root $HEUR_ROOT --heur-nodes $HEUR_NODES \
    --feas-root $FEAS_ROOT --feas-nodes $FEAS_NODES --inherit-cols $INHERIT_COLS \
    --pricing-method $PRICING_METHOD --pricing-greedy-alpha $ALPHA_PRI \
    --pricing-greedy-max-cols $GREEDY_MAX_COLS --pricing-max-cols-per-iter $MAX_COLS_PER_ITER \
    --pricing-exact-time $PRICING_EXACT_TIME --branching-variable $BRANCHING_VARIABLE --preproc-off

    # ILP
    time $BIN -s feas-ilp -f "$INPUT/$LINE" -o "$OUT/" -t $TIME_LIMIT --verbose $VERBOSE --heur-root $HEUR_ROOT

    # Maximum stable set
    timeout 120m $BIN -s feas-enum -f "$INPUT/$LINE" -o "$OUT/" --verbose $VERBOSE

done < "$INSTANCES"
