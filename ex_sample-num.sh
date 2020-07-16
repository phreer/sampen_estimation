#!/bin/bash
set x
FILENAME="data/power_two/ecg_114157_2p17.txt"
DATE=`date +"%Y-%m-%d"`
OUTPUT_PATH="result/linux/$DATE"
if [ ! -d $OUTPUT_PATH ]; then
  mkdir $OUTPUT_PATH
fi

SAMPLE_SIZE=`python -c "print(2 ** 11)"`
M=4
R=30
OUTPUT="$OUTPUT_PATH/sample-num_m"$M"_r"$R"_sample-size2p11".txt
for i in `seq 12 1 14`
do
    SAMPLE_NUM=`python -c "print(2 ** $i)"`
    build_linux/bin/sampen \
        -m $M -r $R -sample_num $SAMPLE_NUM \
        -sample_size $SAMPLE_SIZE \
        -filename $FILENAME >> $OUTPUT
done