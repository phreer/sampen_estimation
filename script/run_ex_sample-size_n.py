import os 
import argparse 
from sqlalchemy import create_engine 
from sqlalchemy.orm import sessionmaker 

import ex_error_file 
import database as db 

parser = argparse.ArgumentParser() 
parser.add_argument('-m', type=int, required=True, help='The template length m. ')
parser.add_argument('-r', type=float, default=0.1, help='The threshold r. ')
parser.add_argument('-sn', '--sample-num', type=int, required=True, 
                    help='The number of results to be averaged. ')
parser.add_argument('-if', '--input-format', type=str, default='simple', 
                    choices=['simple', 'multi-record'], 
                    help='The format of the input file. If set to be simple, then there is one record '
                    'in the file, and each line contains a single column. If set to be multi-record, '
                    'then multiple records may be contained in the file, and each there are NUM_RECORD + 1 '
                    'columns per line, of which the first is the line number and the rest NUM_RECORD '
                    'ones are signals. ')
parser.add_argument('-d', '--database-label', type=str, required=True, 
                    help='The tag to be added to the name of the database file to store the result. ')
parser.add_argument('-l', '--input-offset', type=int, default=0, 
                    help='The offset where the signal starts. ')

def ex_sample_size_length(record, record_name, m, r, sample_num, sess):
    sample_sizes = list() 
    for i in range(1, 4 + 1):
        sample_size = 400 * i 
        sample_sizes.append(sample_size)
    for i in range(1, 10 + 1):
        sample_size = 2000 * i 
        sample_sizes.append(sample_size) 
    for sample_size in sample_sizes:
        ex_error_file.experiment(
            record, record_name, m, r, sample_size, sample_num, sess)

def print_exp_setting(args): 
    print('=' * 76) 
    print('Experiment Settings: ')
    print('\t{:<24}: {}'.format('m', args.m))
    print('\t{:<24}: {}'.format('r', args.r))
    print('\t{:<24}: {}'.format('input-offset', args.input_offset))
    print('\t{:<24}: {}'.format('sample-num', args.sample_num))
    print('=' * 76)

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
    print_exp_setting(args)

    # Initialize database 
    db_file = 'result/linux/ex_sample-size_n_m%d_r%2f_sn%d_l%d_%s.db' % (
        args.m, args.r, args.sample_num, args.input_offset, args.database_label
    )
    if os.path.exists(db_file): 
        print('Database exists. ') 
        exit(1) 
    
    db_uri = 'sqlite:///' + db_file 
    engine = create_engine(db_uri) 
    db.Base.metadata.create_all(engine) 

    Session = sessionmaker(bind=engine) 
    sess = Session() 
    for recordname in recordnames: 
        records = ex_error_file.readfile(
            recordname, n=args.input_offset + 800000, 
            input_format=args.input_format)
        length_record = len(records[0])
        record = records[0]
        for i in range(4):
            n = 10000 * (4 ** i)
            if n < length_record:
                ex_sample_size_length(
                    record[args.input_offset: args.input_offset + n], 
                    recordname, args.m, args.r, args.sample_num, sess)

if __name__ == '__main__':
    main()
