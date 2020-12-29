#!/usr/bin/env python
# coding=utf-8
# author: Tianfeng Shi<jianshuwuxian@163.com>
# 2020/2/26 11:43
#
from data import get_raw_data
import numpy as np
import time
import sys
from contextlib import contextmanager
@contextmanager
def timer(name="time", level="normal"):
    start_time = time.time()
    yield
    end_time = time.time()
    if level == "important":
        print('\033[1;31;40m{}:{}\033[0m'.format(name, end_time - start_time))
    else:
        print('{}:{}'.format(name, end_time - start_time))


def estimate_imputation(snp_num, data_dir):

    '''
    评估填充准确度
    '''
    # 无缺失
    reference_file = data_dir + 'reference.vcf'
    return_data = get_raw_data(file=reference_file)  # 函数返回的三个对象组成一个元组
    reference_data = np.array(return_data[0])

    # 有缺失
    miss_file = data_dir + 'miss_raw.vcf'
    return_data = get_raw_data(file=miss_file)
    with_miss_data = np.array(return_data[0])

    # 填充结果
    imputed_file = data_dir + 'result.vcf'
    return_data = get_raw_data(file=imputed_file)
    imputed_data = np.array(return_data[0])

    total_reference = []
    total_imputed = []

    g_way = []
    
    # 找出缺失位点索引
    for index in range(0, snp_num):
        single = with_miss_data[index]

        # 每个snp的填充准确度 使用对比基因型法 对比字符串
        result = np.where(reference_data[index][np.where(single == 5)[0]] == imputed_data[index][np.where(single == 5)[0]], 1, -1)

        equal = (result == 1)

        accuracy = result[equal].size / float(len(np.where(single == 5)[0]))
        g_way.append(str(accuracy)+'\n')

        total_reference = np.append(total_reference, reference_data[index][np.where(single == 5)[0]])
        total_imputed = np.append(total_imputed, imputed_data[index][np.where(single == 5)[0]])

    # 基因型
    r = np.where(total_reference == total_imputed, 1, -1)
    e = (r == 1)
    total = r[e].size / float(len(total_reference))

    #total = round(total + 0.0000001, 6)

    # 填充准确率写到文件
    accuracy_file = data_dir + 'accuracy.txt'
    with open(accuracy_file, 'w') as f:
        f.writelines(g_way)
        f.write(str(total) + '\n')


if __name__ == '__main__':
    snp_num = int(sys.argv[1])
    data_dir = sys.argv[2]
    estimate_imputation(snp_num=snp_num, data_dir=data_dir)
