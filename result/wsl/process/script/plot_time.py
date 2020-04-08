import csv
import numpy as np

from matplotlib import pyplot as plt

m = 6
if __name__ == '__main__':
    with open('result/wsl/process/219_time_m{}.csv'.format(m)) as f:
        lines = f.readlines()
    reader = csv.reader(lines)
    methods = ['Direct', 'KD-Tree', 'KD-Tree (grid)', 'Range-Tree', 
        'Quasi-Random Sampling', 'KD-Tree Sampling']
    reader = list(reader)
    N = np.array([int(item) for item in reader[0][1:]])
    times = dict()
    for i, method in enumerate(methods):
        times[method] = list()
        for item in reader[i+1][1:]:
            if item != 'None':
                times[method].append(float(item))
        times[method] = np.array(times[method])
        plt.loglog(N[:len(times[method])], times[method])
    plt.xlabel('Data Length N')
    plt.ylabel('Time (seconds)')
    plt.legend(methods)
    plt.savefig('result/wsl/process/fig/219_time_m{}.pdf'.format(m))
    plt.show()