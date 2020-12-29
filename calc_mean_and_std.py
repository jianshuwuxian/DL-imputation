#!/usr/bin/env python

import sys
import numpy as np
import math


def calc_mean_std(file_name):
    data = []
    with open(file_name, 'r') as f:
        lines = f.readlines()
        for line in lines:
            data.append(float(line.strip('\n')))
    SD_P = np.std(data, ddof=1)
    n = len(data)
    res = SD_P/math.sqrt(n)
    mean_res = np.mean(data)
    with open(file_name, 'a') as f:
        f.write('\n\nmean:'+str(mean_res)+'\n')
        f.write(str(res)+'\n')

if __name__ == '__main__':
    f_name = sys.argv[1]
    calc_mean_std(file_name=f_name)

