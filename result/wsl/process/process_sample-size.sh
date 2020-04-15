#!/bin/bash
M=6
FILENAME="result/wsl/219_sample-size_m"$M"_r30_sample-num2p10.txt"
OUTPUT=processed_sample-size_m$M"_sample-num2p10.txt"
for method in "Quasi-random" "KD Tree Sampling"
do
    echo "Method: $method" >> $OUTPUT
    cat $FILENAME | grep -A 1 "$method" | grep "Error" >> $OUTPUT
done