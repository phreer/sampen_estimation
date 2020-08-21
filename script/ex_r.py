import os 
import math 
import argparse

from ex_error import experiment 

parser = argparse.ArgumentParser(description='Do an experiment for the effect '
                                 'of the parameter `r`')
parser.add_argument('-n', type=int, default=0, help='The length of data to '
                    'compute, 0 for total file. ')
parser.add_argument('-m', type=int, required=True, help='The template length m. ')
parser.add_argument('-ss', '--sample-size', type=int, required=True, 
                    help='The number of points to sample. ')
parser.add_argument('-sn', '--sample-num', type=int, required=True, 
                    help='The number of results to be averaged. ')
parser.add_argument('--interval', type=float, required=True, 
                    help='The interval of r. ')
def main():
    filenames = ['data.PhysioNet/chfdb/chf01.txt', 
                 'data.PhysioNet/ltafdb/00.txt', 
                 'data.PhysioNet/ltstdb/s20011.txt', 
                 'data.PhysioNet/mghdb/mgh001.txt', 
                 'data.PhysioNet/mit-bih-long-term-ecg-database-1.0.0/14046.txt', 
                 'data.PhysioNet/pink/pink_noise-2000000.txt', 
                 'data.PhysioNet/gaussian/gaussian_noise-2000000.txt', 
                 ('data.PhysioNet/surrogate-data-with-correlations-trends-and-'
                 'nonstationarities-1.0.0/tns/d2h4pd050918s_2.txt')]
    args = parser.parse_args()
    interval = args.interval
    for filename in filenames: 
        for i in range(math.floor(1 / interval) + 1): 
            r = i * interval
            experiment(filename, args.n, args.m, r, args.sample_size, 
                       args.sample_num, input_format='multi-record')

if __name__ == '__main__': 
    main() 