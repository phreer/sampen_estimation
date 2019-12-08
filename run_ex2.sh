#!/bin/bash
RECORD=14157
EXEC_PATH=./build/bin
FSAMPEN=$EXEC_PATH/sampen
DATA_PATH=./data
r=30
m=5
for i in {0..19}
    do
    N=$[ (i % 5 * 2 + 1) * (10 ** (i / 5 + 2)) ]
    $FSAMPEN -m $m -r $r "$DATA_PATH/$RECORD$N"1.txt
done