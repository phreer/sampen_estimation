#!/bin/bash
set x
M=5
R=30
FILENAME="data/power_two/ecg_114157_2p17.txt"
OUTPUT="result/wsl/219_sample-size_m"$M"_r"$R"_sample-num2p10".txt
SAMPLE_NUM=`python -c "print(2 ** 8)"`

for i in `seq 14 1 15`
do
    SAMPLE_SIZE=`python -c "print(2 ** $i)"`
    build/bin/sampen \
        -m $M -r $R -sample_num $SAMPLE_NUM \
        -sample_size $SAMPLE_SIZE \
        -filename $FILENAME >> $OUTPUT
done