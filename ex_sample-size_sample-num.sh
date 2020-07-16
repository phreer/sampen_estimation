#!/bin/bash
DATE=`date +"%Y-%m-%d"`
OUTPUT_PATH="result/linux/$DATE"
if [ ! -d $OUTPUT_PATH ]; then
  mkdir $OUTPUT_PATH
fi
FILENAME="data/power_two/ecg_114157_2p17.txt"

M=6
R=30
OUTPUT="$OUTPUT_PATH/sample-num_sample-size_m"$M"_r"$R"_2p17.txt"
for i in `seq 8 1 17`
do
    SAMPLE_SIZE=`python -c "print(2 ** $i)"`
    SAMPLE_NUM=`python -c "print(2 ** (17 - $i))"`
    build_linux/bin/sampen \
        -m $M -r $R -sample_num $SAMPLE_NUM \
        -sample_size $SAMPLE_SIZE \
        -filename $FILENAME >> $OUTPUT
done