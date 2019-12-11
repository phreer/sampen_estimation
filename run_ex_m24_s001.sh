#!/bin/bash
set +x

RECORD=14157
EXEC_PATH=./build/bin
FSAMPEN=$EXEC_PATH/sampen
DATA_PATH=./data
r=30
s=0.1

for m in {2..4}
do
    for i in {0..19}
    do
        N=$[ (i % 5 * 2 + 1) * (10 ** (i / 5 + 2)) ]
        echo $FSAMPEN -m $m -r $r "$DATA_PATH/$RECORD$N"1.txt
        $FSAMPEN -m $m -r $r -s $s "$DATA_PATH/$RECORD$N"1.txt
    done
done