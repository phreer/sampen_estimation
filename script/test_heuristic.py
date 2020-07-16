import os 
import glob 
from math import sqrt 

from matplotlib import pyplot as plt
from sampen import * 
from ex_error import readfile, variance 


signals = ['AF', 'CHF', 'Health']
input_dir = '../data_int/'
output_dir = 'err'

sample_sizes = [i * 1000 for i in range(1, 30)]
sample_num = 8 
m = 3 
r = 0.2 

for signal in signals:
  filenames = glob.glob(os.path.join(input_dir, signal, '*.txt'))
  os.makedirs(os.path.join(output_dir, signal), exist_ok=True)
  for filename in filenames:
    data = readfile(filename) 
    std = sqrt(variance(data))
    r_scaled = std * r 
    n = len(data)

    sampen_d, B_d, A_d = compute_sampen_direct(data, m, r_scaled)
    B_d_norm = B_d / n / n 
    A_d_norm = A_d / n / n

    A_errs = []
    B_errs = []
    for sample_size in sample_sizes:
      sampen_qr, B_qr, A_qr = compute_sampen_qr(
          data, m, r_scaled, sample_size, sample_num, True)

      B_qr_norm = B_qr / sample_size / sample_size 
      A_qr_norm = A_qr / sample_size / sample_size 
      B_err = B_qr_norm - B_d_norm 
      A_err = A_qr_norm - A_d_norm 
      B_errs.append(abs(B_err))
      A_errs.append(abs(A_err))
      sampen_err = sampen_qr - sampen_d 

      print('Filename: ', filename)
      print('Sampen (m=%d, r=%f, n=%d): %f, A: %f, B: %f'\
          % (m, r, n, sampen_d, A_d_norm, B_d_norm))
      print('Sample size: %f, error (A): %f, error (B): %f, error (sampen): %f'\
          % (sample_size, A_err, B_err, sampen_err))
    save_path = os.path.join(output_dir, signal, os.path.basename(filename) + '.pdf')
    fig = plt.figure(figsize=[8, 6], dpi=100)
    ax = fig.add_subplot()
    ax.plot(sample_sizes, B_errs, label='Error (B)')
    ax.plot(sample_sizes, A_errs, label='Error (A)')
    ax.set_title('%s: %s' % (signal, os.path.basename(filename)))
    ax.set_yscale('log')
    fig.savefig(save_path)
    plt.close(fig)
