#!/bin/bash

# experiments/run_exp.sh: Run all experiments in the experiments folder

declare -a experiments=("exp11" "exp12" "exp10" "exp7" "exp6")

for exp in "${experiments[@]}"
do
    echo "Running experiment: $exp"
    cd $exp
    ./${exp}.sh
    cd ..
done