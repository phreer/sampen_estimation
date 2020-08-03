#!/bin/bash 
# Run experiment for the relation of sample size to convergent and length of 
# the data. 

# 2020-7-28 
set -o noclobber 
python script/run_ex_sample-size_n.py \
    -m 4 \
    -r 0.2 \
    --input-format multi-record \
    --sample-num 10 