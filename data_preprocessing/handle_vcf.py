#!/usr/bin/env python
# coding=utf-8
# author: Tianfeng Shi<jianshuwuxian@163.com>
# 2020/2/4 15:36
#
import collections
import numpy as np
import os
import sys
import time
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


def remove_vcf_header(data_dir):
    """
    去掉注释行以及每个标记的前9列信息
    :return:
    """

    # 文件名
    file_name = data_dir + "miss_raw.vcf"
    target_file_name = data_dir + "miss_raw_without_header.txt"
    header_file_name = data_dir + "header.txt"

    header_line_num = 0
    with open(file_name, 'r') as f:
        # VCF文件注释部分
        header_list = []

        line = f.readline()
        while line.startswith("##"):
            header_list.append(line)
            line = f.readline()
        header_list.append(line)

        header_line_num = len(header_list)

    with open(header_file_name, 'w') as f:
        f.writelines(header_list)

    # os.system("sed -i '1," + str(header_line_num) + "d' " + file_name)
    with timer():
        os.system("tail -n +" + str(header_line_num+1) + " " + file_name + " > " + target_file_name)


if __name__ == '__main__':
    data_dir = sys.argv[1]
    remove_vcf_header(data_dir=data_dir)
