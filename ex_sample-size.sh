#!/bin/bash
set x
DATE=`date +"%Y-%m-%d"`
OUTPUT_PATH="result/linux/$DATE"
if [ ! -d $OUTPUT_PATH ]; then
  mkdir $OUTPUT_PATH
fi

FILENAME="data/power_two/ecg_114157_2p17.txt"
M=6
R=30
SAMPLE_NUM_P=8
SAMPLE_NUM=`python -c "print(2 ** 8)"`
OUTPUT="$OUTPUT_PATH/sample-size_m"$M"_r"$R"_sample-num2p$SAMPLE_NUM_P".txt

for i in `seq 8 1 14`
do
    SAMPLE_SIZE=`python -c "print(2 ** $i)"`
    build_linux/bin/sampen \
        -m $M -r $R -sample_num $SAMPLE_NUM \
        -sample_size $SAMPLE_SIZE \
        -filename $FILENAME >> $OUTPUT
done