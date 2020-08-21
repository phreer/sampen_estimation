#!/bin/bash 
# Run experiment for the relation of sample size to convergent and length of 
# the data. 

# 2020-7-28 
set -o noclobber 
label=20-08-03
result_file=result/linux/output_sample-size_n_m4_r03_sn10_ecg_$label.txt 
if [ -e $result_file ]; then 
    echo "Output file [$result_file] exists. " >&2 
    exit 1
fi 
python script/run_ex_sample-size_n.py \
    -m 4 \
    -r 0.3 \
    --input-format multi-record \
    --database-label $label \
    --sample-num 10 \
    > $result_file 