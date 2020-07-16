#!/bin/sh

python script/ex_error.py -m 2 -r 0.1 -sample_size 1024 -sample_num 32 &
python script/ex_error.py -m 3 -r 0.1 -sample_size 1024 -sample_num 32 &
python script/ex_error.py -m 4 -r 0.1 -sample_size 1024 -sample_num 32 &
python script/ex_error.py -m 5 -r 0.1 -sample_size 1024 -sample_num 32 &
python script/ex_error.py -m 6 -r 0.1 -sample_size 1024 -sample_num 32 &

python script/ex_error.py -m 2 -r 0.1 -sample_size 2048 -sample_num 32 &
python script/ex_error.py -m 3 -r 0.1 -sample_size 2048 -sample_num 32 &
python script/ex_error.py -m 4 -r 0.1 -sample_size 2048 -sample_num 32 &
python script/ex_error.py -m 5 -r 0.1 -sample_size 2048 -sample_num 32 &
python script/ex_error.py -m 6 -r 0.1 -sample_size 2048 -sample_num 32 &

python script/ex_error.py -m 2 -r 0.1 -sample_size 4096 -sample_num 32 &
python script/ex_error.py -m 3 -r 0.1 -sample_size 4096 -sample_num 32 &
python script/ex_error.py -m 4 -r 0.1 -sample_size 4096 -sample_num 32 &
python script/ex_error.py -m 5 -r 0.1 -sample_size 4096 -sample_num 32 &
python script/ex_error.py -m 6 -r 0.1 -sample_size 4096 -sample_num 32 &