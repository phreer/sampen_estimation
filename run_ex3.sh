#!/bin/bash
RECORD=14157
EXEC_PATH=/Volumes/Data/workspace/sampen_range_tree_random
FSAMPEN=$EXEC_PATH/sampen_range_tree
DATA_PATH=/Volumes/Data/workspace/sampen_kdtree/data
r=30
m=5
for i in {17..19}
    do
    N=$[ (i % 5 * 2 + 1) * (10 ** (i / 5 + 2)) ]
    $FSAMPEN -m $m -r $r "$DATA_PATH/$RECORD$N"1.txt
done