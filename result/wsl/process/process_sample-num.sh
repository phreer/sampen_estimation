#!/bin/bash
M=6
FILENAME="result/wsl/218_sample-num_m"$M"_r30_sample-size2p11.txt"
OUTPUT=processed_sample-num_m$M"_sample-size2p10.txt"
for method in "Quasi-random" "KD Tree Sampling"
do
    echo "Method: $method" >> $OUTPUT
    cat $FILENAME | grep -A 1 "$method" | grep "Error" >> $OUTPUT
done