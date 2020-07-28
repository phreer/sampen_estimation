import os
import sys
import math
import glob
import argparse
import time

import sampen

parser = argparse.ArgumentParser('ex_error.py')
parser.add_argument('-n', type=int, default=0, help='The length of data to compute, 0 for total file. ')
parser.add_argument('-m', type=int, required=True, help='The template length m. ')
parser.add_argument('-r', type=float, default=0.1, help='The threshold r. ')
parser.add_argument('-sample_size', type=int, required=True, help='The number of points to sample. ')
parser.add_argument('-sample_num', type=int, required=True, help='The number of results to be averaged. ')
parser.add_argument('-input_format', type=str, default='simple', choices=['simple', 'multi-record'], 
                    help='The format of the input file. If set to be simple, then there is one record '
                    'in the file, and each line contains a single column. If set to be multi-record, '
                    'then multiple records may be contained in the file, and each there are NUM_RECORD + 1 '
                    'columns per line, of which the first is the line number and the rest NUM_RECORD '
                    'ones are signals. ')

def readfile(filename, n=None, input_format='simple'):
    """
    Read a file to a list or lists of integers with length `n`, according to 
    the `input_format`.. 

    @note: The file is supposed to be lines of integers. 
    """
    with open(filename) as f:
        lines = f.readlines()

    if input_format == 'simple': 
        result = list()
        count = 0
        for line in lines:
            if (line.strip()):
                result.append(int(line))
                count += 1
                if count == n: 
                    break 
    elif input_format == 'multi-record':
        result = list() 
        count = 0
        for line in lines:
            if line.strip():
                result.append(line.split()[1:])
                count += 1
                if count == n:
                    break 
        result = list(map(list, zip(*result)))
    
    return result

def variance(data):
    m = sum(data) / len(data)
    data_ = [(d - m) ** 2 for d in data]
    return sum(data_) / len(data)

def experiment(filename, n, m, r, sample_size, sample_num):
    data = readfile(filename)
    var = variance(data)
    r = int(r * math.sqrt(var))
    if n > 0: data = data[:n]
    n = len(data)

    print('=' * 70)
    print('{:<16}: {}'.format('filename', filename)) 
    print('{:<16}: {}'.format('data length', len(data)))
    print('{:<16}: {}'.format('template length', m))
    print('{:<16}: {}'.format('threshold (r)', r))
    print('{:<16}: {}'.format('sample size', sample_size))
    print('{:<16}: {}'.format('sample num', sample_num))
    print('{:<16}: {:.4}'.format('variance', var))
    print('=' * 70)

    t = time.time()
    sampen_d, a_d, b_d = sampen.compute_sampen_direct(data, m, r)
    t = time.time() - t
    print('method: direct')
    print('\t{:<30}: {:.4f}'.format('time', t))
    print('\t{:<30}: {:.6}'.format('sample entropy (direct)', sampen_d))
    print('\t{:<30}: {:.2e}'.format('a', a_d))
    print('\t{:<30}: {:.2e}'.format('b', b_d))

    t = time.time()
    sampen_q1, a_q1, b_q1 = sampen.compute_sampen_qr(data, m, r, sample_size, sample_num, False)
    t = time.time() - t
    print('method: quasi-Monte Carlo')
    print('\t{:<30}: {:.4f}'.format('time', t))
    print('\t{:<30}: {:.6f}'.format('sample entropy (qr1)', sampen_q1))
    err = (sampen_q1 - sampen_d) / (sampen_d + 1e-8)
    print('\t{:<30}: {:.2e}'.format('absolute error', sampen_q1 - sampen_d))
    print('\t{:<30}: {:.2e}'.format('relative error', err))
    print('\t{:<30}: {:.2e}'.format('a', a_q1))
    print('\t{:<30}: {:.2e}'.format('b', b_q1))

    print('\t{:<30}: {:.2e}'.format('error of a', a_q1 / sample_size / sample_size - a_d / n / n))
    print('\t{:<30}: {:.2e}'.format('error of b', b_q1 / sample_size / sample_size - b_d / n / n))

    t = time.time()
    sampen_q2, a_q2, b_q2 = sampen.compute_sampen_qr(data, m, r, sample_size, sample_num, True)
    t = time.time() - t
    print('method: quasi-Monte Carlo (presort)')
    print('\t{:<30}: {:.4f}'.format('time', t))
    print('\t{:<30}: {:.6f}'.format('sample entropy', sampen_q2))
    err = (sampen_q2 - sampen_d) / (sampen_d + 1e-8)
    print('\t{:<30}: {:.2e}'.format('absolute error', sampen_q2 - sampen_d))
    print('\t{:<30}: {:.2e}'.format('relative error', err))
    print('\t{:<30}: {:.2e}'.format('a', a_q2))
    print('\t{:<30}: {:.2e}'.format('b', b_q2))
    print('\t{:<30}: {:.2e}'.format('error of a', a_q2 / sample_size / sample_size - a_d / n / n))
    print('\t{:<30}: {:.2e}'.format('error of a', b_q2 / sample_size / sample_size - b_d / n / n))

    t = time.time()
    sampen_u, a_u, b_u = sampen.compute_sampen_uniform(data, m, r, sample_size, sample_num)
    t = time.time() - t
    print('method: Monte Carlo')
    print('\t{:<30}: {:.4f}'.format('time', t))
    print('\t{:<30}: {:.6f}'.format('sample entropy', sampen_u))
    err = (sampen_u - sampen_d) / (sampen_d + 1e-8)
    print('\t{:<30}: {:.2e}'.format('absolute error', sampen_u - sampen_d))
    print('\t{:<30}: {:.2e}'.format('relative error', err))
    print('\t{:<30}: {:.2e}'.format('a', a_u))
    print('\t{:<30}: {:.2e}'.format('b', b_u))
    print('\t{:<30}: {:.2e}'.format('error of a', a_u / sample_size / sample_size - a_d / n / n))
    print('\t{:<30}: {:.2e}'.format('error of a', b_u / sample_size / sample_size - b_d / n / n))

def main():
    signals = ['AF', 'CHF', 'Health']
    args = parser.parse_args()

    for signal in signals:
        input_dir = os.path.join('data_int', signal)
        output_dir = os.path.join('result/2020-6-11/', 
                                  'ss{}_sn{}'.format(args.sample_size, args.sample_num), 
                                  'r{:.1f}'.format(args.r), 
                                  'm{}'.format(args.m), signal)
        os.makedirs(output_dir, exist_ok=True)

        filenames = glob.glob(os.path.join(input_dir, '*.txt'))
        for filename in filenames:
            output_filename = os.path.join(os.path.basename(filename))
            f = open(os.path.join(output_dir, output_filename), 'w')
            # sys.stdout = f
            experiment(filename, args.n, args.m, args.r, args.sample_size, args.sample_num)
            f.close()

if __name__ == '__main__':
    main()
