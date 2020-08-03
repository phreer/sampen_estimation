import argparse 

import ex_error_file 

parser = argparse.ArgumentParser() 
parser.add_argument('-m', type=int, required=True, help='The template length m. ')
parser.add_argument('-r', type=float, default=0.1, help='The threshold r. ')
parser.add_argument('-sn', '--sample-num', type=int, required=True, help='The number of results to be averaged. ')
parser.add_argument('-if', '--input-format', type=str, default='simple', choices=['simple', 'multi-record'], 
                    help='The format of the input file. If set to be simple, then there is one record '
                    'in the file, and each line contains a single column. If set to be multi-record, '
                    'then multiple records may be contained in the file, and each there are NUM_RECORD + 1 '
                    'columns per line, of which the first is the line number and the rest NUM_RECORD '
                    'ones are signals. ')

def ex_sample_size_length(record, record_name, m, r, sample_num):
    sample_num = 10
    sample_sizes = list() 
    for i in range(1, 4 + 1):
        sample_size = 400 * i 
        sample_sizes.append(sample_size)
    for i in range(1, 20 + 1):
        sample_size = 2000 * i 
        sample_sizes.append(sample_size) 
    for sample_size in sample_sizes:
        ex_error_file.experiment(record, record_name, m, r, sample_size, sample_num)

info = (
'This file contains the result of an experiment that shows the relation \n'
'between the error and the sample size for different signals and different \n' 
'length of the signal. Note that the sample rates are variable. '
)
def main():
    recordnames = ['data.PhysioNet/chfdb/chf01.txt', 
                   'data.PhysioNet/ltafdb/00.txt', 
                   'data.PhysioNet/ltstdb/s20011.txt',
                   'data.PhysioNet/mit-bih-long-term-ecg-database-1.0.0/14046.txt']
    args = parser.parse_args() 
    print(info)
    for recordname in recordnames: 
        records = ex_error_file.readfile(recordname, input_format=args.input_format)
        length_record = len(records[0])
        record = records[0]
        for i in range(7):
            n = 10000 * (2 ** i)
            if n < length_record:
                ex_sample_size_length(record[: n], recordname, args.m, args.r, args.sample_num)

if __name__ == '__main__':
    main()
