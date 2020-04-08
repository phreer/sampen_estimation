#!/bin/bash
FILENAME="data/power_two/ecg_114157_2p17.txt"
M=4
R=30
OUTPUT="result/wsl/219_sample-num_sample-size_m"$M"_r"$R"_2p17.txt"
for i in `seq 8 1 17`
do
    SAMPLE_SIZE=`python -c "print(2 ** $i)"`
    SAMPLE_NUM=`python -c "print(2 ** (17 - $i))"`
    build/bin/sampen \
        -m $M -r $R -sample_num $SAMPLE_NUM \
        -sample_size $SAMPLE_SIZE \
        -filename $FILENAME >> $OUTPUT
done