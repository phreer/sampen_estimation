#!/bin/bash 
# Run experiment for the relation of sample size to convergent and length of 
# the data. 

# 2020-7-28 
set -o noclobber 
python script/run_ex_sample-size_n.py \
    -r 0.1 \
    -m 4 \
    --input-format multi-record \
    --sample-num 10 