#!/bin/bash

## Usage:
#    test_sims.sh <options> <n trials>

## Test forward sim model.
if [ "$1" == "-new" ]
then
    rm output/means.out
    rm sims/alpha.dat
    rm sims/ch_params.csv
    rm output/quantiles.out
    rm sims/rsq.csv
    echo "REMOVED OLD SIM FILES"
fi

i=1
while [ $i -le $2 ]
do
    echo "Simulation $i of $2"
    cd sims
    R --no-save --slave < sim_dat.R

    cd ..
    ./gb -s 10 100
    
    R --no-save --slave < R/record_alpha_mean.R
    i=$(( $i + 1 ))
done

R --no-save --slave < R/plot_alpha_mean.R
# acroread pred_obs_alpha.pdf&

