import os
import sys
import math
import glob
import argparse
import time

import sampen
from utils import readfile, variance

parser = argparse.ArgumentParser('ex_error.py')
parser.add_argument('-f', '--filename', type=str, required=True, 
                    help='The filename of the input signal. ')
parser.add_argument('-n', type=int, default=0, help='The length of data to compute, 0 for total file. ')
parser.add_argument('-m', type=int, required=True, help='The template length m. ')
parser.add_argument('-r', type=float, default=0.1, help='The threshold r. ')
parser.add_argument('-ss', '--sample-size', type=int, required=True, help='The number of points to sample. ')
parser.add_argument('-sn', '--sample-num', type=int, required=True, help='The number of results to be averaged. ')
parser.add_argument('-if', '--input-format', type=str, default='simple', 
                    choices=['simple', 'multi-record'], 
                    help='The format of the input file. If set to be simple, then there is one record '
                    'in the file, and each line contains a single column. If set to be multi-record, '
                    'then multiple records may be contained in the file, and each there are NUM_RECORD + 1 '
                    'columns per line, of which the first is the line number and the rest NUM_RECORD '
                    'ones are signals. ')


def experiment(filename, n, m, r, sample_size, sample_num, 
               input_format='simple'):
    data = readfile(filename, n=n, input_format=input_format)[0]

    var = variance(data)
    r = int(r * math.sqrt(var))
    n = len(data)

    print('=' * 72)
    print('{:<26}: {}'.format('filename', filename)) 
    print('{:<26}: {}'.format('data length', len(data)))
    print('{:<26}: {}'.format('template length', m))
    print('{:<26}: {}'.format('threshold (r)', r))
    print('{:<26}: {}'.format('sample size', sample_size))
    print('{:<26}: {}'.format('sample num', sample_num))
    print('{:<26}: {:.4}'.format('variance', var))
    print('-' * 72)

    t = time.time()
    sampen_d, a_d, b_d = sampen.compute_sampen_direct(data, m, r)
    t = time.time() - t
    normalizer_d = (n - m) ** 2
    print('method: direct')
    print('\t{:<30}: {:.4f}'.format('time', t))
    print('\t{:<30}: {:.6}'.format('sample entropy', sampen_d))
    a_d_normalized = a_d / normalizer_d
    b_d_normalized = b_d / normalizer_d
    print('\t{:<30}: {:.2e} ({:.2e})'.format('a', a_d, a_d_normalized))
    print('\t{:<30}: {:.2e} ({:.2e})'.format('b', b_d, b_d_normalized))

    t = time.time()
    sampen_q1, a_q1, b_q1 = sampen.compute_sampen_qr(
        data, m, r, sample_size, sample_num, False)
    t = time.time() - t
    normalizer = (sample_size - 1) ** 2
    print('method: quasi-Monte Carlo')
    print('\t{:<30}: {:.4f}'.format('time', t))
    print('\t{:<30}: {:.6f}'.format('sample entropy', sampen_q1))
    err = (sampen_q1 - sampen_d) / (sampen_d + 1e-8)
    print('\t{:<30}: {:.2e}'.format('absolute error', sampen_q1 - sampen_d))
    print('\t{:<30}: {:.2e}'.format('relative error', err))
    a_q1_normalized = a_q1 / normalizer
    b_q1_normalized = b_q1 / normalizer 
    print('\t{:<30}: {:.2e} ({:.2e})'.format('a', a_q1, a_q1_normalized))
    print('\t{:<30}: {:.2e} ({:.2e})'.format('b', b_q1, b_q1_normalized))
    print('\t{:<30}: {:.2e}'.format('error_a', a_q1_normalized - a_d_normalized))
    print('\t{:<30}: {:.2e}'.format('error_b', b_q1_normalized - b_d_normalized))
    err_a_rel = (a_q1_normalized - a_d_normalized) / abs(a_d_normalized + 1e-10)
    err_b_rel = (b_q1_normalized - b_d_normalized) / abs(b_d_normalized + 1e-10)
    print('\t{:<30}: {:.2e}'.format('error_a (relative)', err_a_rel))
    print('\t{:<30}: {:.2e}'.format('error_b (relative)', err_b_rel))

    t = time.time()
    sampen_q2, a_q2, b_q2 = sampen.compute_sampen_qr(data, m, r, sample_size, sample_num, True)
    t = time.time() - t
    print('method: quasi-Monte Carlo (presort)')
    print('\t{:<30}: {:.4f}'.format('time', t))
    print('\t{:<30}: {:.6f}'.format('sample entropy', sampen_q2))
    err = (sampen_q2 - sampen_d) / (sampen_d + 1e-8)
    print('\t{:<30}: {:.2e}'.format('absolute error', sampen_q2 - sampen_d))
    print('\t{:<30}: {:.2e}'.format('relative error', err))
    a_q2_normalized = a_q2 / normalizer
    b_q2_normalized = b_q2 / normalizer 
    print('\t{:<30}: {:.2e} ({:.2e})'.format('a', a_q2, a_q2_normalized))
    print('\t{:<30}: {:.2e} ({:.2e})'.format('b', b_q2, b_q2_normalized))
    print('\t{:<30}: {:.2e}'.format('error_a', a_q2_normalized - a_d_normalized))
    print('\t{:<30}: {:.2e}'.format('error_b', b_q2_normalized - b_d_normalized))
    err_a_q2_rel = (a_q2_normalized - a_d_normalized) / abs(a_d_normalized + 1e-10)
    err_b_q2_rel = (b_q2_normalized - b_d_normalized) / abs(b_d_normalized + 1e-10)
    print('\t{:<30}: {:.2e}'.format('error_a (relative)', err_a_q2_rel))
    print('\t{:<30}: {:.2e}'.format('error_b (relative)', err_b_q2_rel))

    t = time.time()
    sampen_u, a_u, b_u = sampen.compute_sampen_uniform(data, m, r, sample_size, sample_num)
    t = time.time() - t
    print('method: Monte Carlo')
    print('\t{:<30}: {:.4f}'.format('time', t))
    print('\t{:<30}: {:.6f}'.format('sample entropy', sampen_u))
    err = (sampen_u - sampen_d) / (sampen_d + 1e-8)
    print('\t{:<30}: {:.2e}'.format('absolute error', sampen_u - sampen_d))
    print('\t{:<30}: {:.2e}'.format('relative error', err))
    a_u_normalized = a_u / normalizer
    b_u_normalized = b_u / normalizer 
    print('\t{:<30}: {:.2e} ({:.2e})'.format('a', a_u, a_u_normalized))
    print('\t{:<30}: {:.2e} ({:.2e})'.format('b', b_u, b_u_normalized))
    print('\t{:<30}: {:.2e}'.format('error_a', a_u_normalized - a_d_normalized))
    print('\t{:<30}: {:.2e}'.format('error_b', b_u_normalized - b_d_normalized))
    err_a_u_rel = (a_u_normalized - a_d_normalized) / abs(a_d_normalized + 1e-10)
    err_b_u_rel = (b_u_normalized - b_d_normalized) / abs(b_d_normalized + 1e-10)
    print('\t{:<30}: {:.2e}'.format('error_a (relative)', err_a_u_rel))
    print('\t{:<30}: {:.2e}'.format('error_b (relative)', err_b_u_rel))
    print('=' * 72) 

def main():
    args = parser.parse_args()
    experiment(args.filename, 
               args.n, 
               args.m, 
               args.r, 
               args.sample_size, 
               args.sample_num, 
               args.input_format)

if __name__ == '__main__':
    main()
