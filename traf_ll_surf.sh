#!/bin/bash

## Likelihood suface of d and e GRAV params

for d in $(seq 0 0.1 3)
do
    for e in $(seq 0 0.033 1)
    do
        ll=`./gb -ll 10 100 $d $e`
        echo "$d    $e  $ll"
    done
done
exit 0;
