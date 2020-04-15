#!/bin/bash
M=6
FILENAME="result/wsl/217_time_m"$M"_r30_sample-size2p11.txt"
OUTPUT=processed_time_m$M.txt
for method in "Direct" "kd tree:" "kd tree (grid)" "Range Tree" "Quasi-random" "KD Tree Sampling"
do
    echo "Method: $method" >> $OUTPUT
    cat $FILENAME | grep -B 1 "$method" | grep "time" >> $OUTPUT
done