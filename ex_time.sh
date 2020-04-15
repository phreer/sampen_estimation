#!/bin/bash
set x
M=6
R=30
OUTPUT="result/wsl/414/time_m"$M"_r"$R"_sample-size2p11_sample-num2p10".txt
for i in `seq 13 1 19`
do
    FILENAME="data/power_two/ecg_114157_2p"$i".txt"
    SAMPLE_SIZE=`python -c "print(2 ** 11)"`
    # SAMPLE_NUM=`python -c "print(2 ** ($i - 11))"`
    SAMPLE_NUM=1024
    build/bin/sampen \
        -m $M -r $R -sample_num $SAMPLE_NUM \
        -sample_size $SAMPLE_SIZE \
        -filename $FILENAME >> $OUTPUT
done