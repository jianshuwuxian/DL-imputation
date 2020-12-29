#!/usr/bin/env python
# coding=utf-8
# author: Tianfeng Shi<jianshuwuxian@163.com>
# 2020/4/29 11:13
#
import os
import time
import sys


def timer(func):
    def wrapper(data_dir):
        start = time.time()
        func(data_dir)
        end = time.time()
        print("{}: {}".format(func.__name__, end-start))
    return wrapper


def remove_header_and_info(data_dir):
    """
    去掉注释行以及每个标记的前9列信息
    :return:
    """

    # 文件名
    file_name = data_dir + "miss_raw.vcf"
    target_file_name = data_dir + "miss_raw_data.txt"
    # header_file_name = data_dir + "header.txt"

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

    # with open(header_file_name, 'w') as f:
    #     f.writelines(header_list)

    os.system("tail -n +" + str(header_line_num + 1) + " " + file_name + " | cut -f10- > " + target_file_name)


def encode2emma(data_dir, sub_num):
    file_name = data_dir + 'sub_noinfo_geno_file' + sub_num + '.txt'
    new_file_name = data_dir + 'sub_miss_encoded_emma' + sub_num + '.txt'

    with open(file_name, 'r') as f:
        lines = f.readlines()
        new_lines = []
        for line in lines:
            data = line.strip().split()
            new_data = []
            for d in data:
                if d == '0|0':
                    new_data.append(0)
                elif d == '1|1':
                    new_data.append(1)
                elif d == '1|0' or d == '0|1':
                    new_data.append(0.5)
                else:
                    new_data.append('NA')

            new_lines.append("\t".join(list(map(str, new_data)))+'\n')

    with open(new_file_name, 'w') as f:
        f.writelines(new_lines)


def use_emma_script(script_dir, data_dir, sub_num):
    command = 'Rscript ' + script_dir + 'calc_gmatrix.R ' + data_dir + ' ' + sub_num
    os.system(command)


if __name__ == '__main__':
    script_dir = sys.argv[1]

    data_dir = sys.argv[2]

    sub_num = sys.argv[3]

    # remove_header_and_info(data_dir)

    encode2emma(data_dir, sub_num)

    use_emma_script(script_dir, data_dir, sub_num)

