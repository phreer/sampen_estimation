#!/bin/bash
set x
M=3
R=30
OUTPUT="result/wsl/217_time_m"$M"_r"$R"_sample-size2p11".txt
for i in `seq 19 1 19`
do
    FILENAME="data/power_two/ecg_114157_2p"$i".txt"
    SAMPLE_SIZE=`python -c "print(2 ** 11)"`
    SAMPLE_NUM=`python -c "print(2 ** ($i - 11))"`
    build/bin/sampen \
        -m $M -r $R -sample_num $SAMPLE_NUM \
        -sample_size $SAMPLE_SIZE \
        -filename $FILENAME >> $OUTPUT
done