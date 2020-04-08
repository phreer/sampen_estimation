#!/bin/bash
set x
M=3
R=30
FILENAME="data/power_two/ecg_114157_2p17.txt"
OUTPUT="result/wsl/218_sample-num_m"$M"_r"$R"_sample-size2p11".txt
SAMPLE_SIZE=`python -c "print(2 ** 11)"`

for i in `seq 0 1 15`
do
    SAMPLE_NUM=`python -c "print(2 ** $i)"`
    build/bin/sampen \
        -m $M -r $R -sample_num $SAMPLE_NUM \
        -sample_size $SAMPLE_SIZE \
        -filename $FILENAME >> $OUTPUT
done