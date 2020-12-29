#!/usr/bin/env python
# coding=utf-8
# author: Tianfeng Shi<jianshuwuxian@163.com>
# 2020/2/6 23:29
#
import os
import sys
import subprocess


def job(data_dir, sub_file_num, script_dir):
    for i in range(sub_file_num):
        subprocess.run(["python3", script_dir + "multi_process_data_generator.py", data_dir, str(i+1)])


if __name__ == '__main__':
    data_dir = sys.argv[1]
    sub_file_num = sys.argv[2]
    script_dir = sys.argv[3]
    job(data_dir=data_dir, sub_file_num=sub_file_num, script_dir=script_dir)
