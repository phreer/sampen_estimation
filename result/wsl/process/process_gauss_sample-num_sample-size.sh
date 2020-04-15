#!/bin/bash
M=6
FILENAME="result/wsl/218_gauss_sample-num_sample-size_m"$M"_r30_2p17.txt"
OUTPUT=processed_gauss_sample-num_sample-size_m$M.txt
for method in "Quasi-random" "KD Tree Sampling"
do
    echo "Method: $method" >> $OUTPUT
    cat $FILENAME | grep -A 1 "$method" | grep "Error" >> $OUTPUT
done