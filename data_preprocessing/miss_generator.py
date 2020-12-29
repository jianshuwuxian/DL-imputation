#!/usr/bin/env python
# coding=utf-8
# author: Tianfeng Shi<jianshuwuxian@163.com>
# 2019/11/14 15:00
#
import collections
import random
import time
from contextlib import contextmanager
import numpy


# @contextmanager
# def timer():
#     start_time = time.time()
#     yield
#     end_time = time.time()
#     print('{}'.format(end_time - start_time))


def make_missing_data(file_name, isdiv=True):
    if isdiv:
        f = open(file_name, 'r')

        # VCF文件注释部分
        header_list = []

        line = f.readline()
        while line.startswith("##"):
            header_list.append(line)
            line = f.readline()
        header_list.append(line)

        line = line.strip('#\n')
        title_list = line.split('\t')

        data_list = []
        for l in f:
            data_line = l.strip('\n').split('\t')
            line_dict = collections.OrderedDict()

            for index, title in enumerate(title_list[:9]):
                line_dict[title] = data_line[index]
            line_dict['samples'] = data_line[9:]
            data_list.append(line_dict)

        f.close()

        # 总位点数
        net_num = len(title_list[9:]) * len(data_list) * 2

        # 缺失率
        random_rate = 0.1

        # random_pos_list = []
        random_pos_list = random.sample(range(0, net_num), int(random_rate * net_num))
        # for i in range(int(random_rate * net_num)):
        #     random_pos = random.randint(0, net_num - 1)
        #     random_pos_list.append(random_pos)

        allele_list = []
        for data in data_list:
            for allele_pair in data['samples']:
                allele_list += str(allele_pair).split('|')

        # 产生缺失
        allele_mat = numpy.mat(allele_list)
        allele_mat[0, random_pos_list] = 5
        allele_list = allele_mat.tolist()[0]

        # for index, allele_base in enumerate(allele_list):
        #     if index in random_pos_list:
        #         allele_list[index] = '.'

        temp_list = []
        for i in range(0, len(allele_list), 2):
            temp_list.append('|'.join(allele_list[i:i+2]))

        new_allele_list = []
        for i in range(0, len(temp_list), len(title_list[9:])):
            new_allele_list.append(temp_list[i:i+len(title_list[9:])])

        for index, data in enumerate(data_list):
            data['samples'] = new_allele_list[index]

        f = open('./data_set/miss_raw.vcf', 'w')

        # 写入注释部分
        f.writelines(header_list)

        new_data_list = []
        for data in data_list:
            one_line = ''
            for v in data.values():

                if isinstance(v, list):
                    for vv in v:
                        one_line += vv + '\t'
                else:
                    one_line += str(v) + '\t'
            one_line = one_line.strip("\t")
            one_line += '\n'

            new_data_list.append(one_line)

        f.writelines(new_data_list)
        f.close()
    else:
        # file_name = './target.vcf'

        f = open(file_name, 'r')

        # VCF文件注释部分
        header_list = []

        line = f.readline()
        while line.startswith("##"):
            header_list.append(line)
            line = f.readline()
        header_list.append(line)

        line = line.strip('#\n')
        title_list = line.split('\t')

        data_list = []
        for l in f:
            data_line = l.strip('\n').split('\t')
            line_dict = collections.OrderedDict()

            for index, title in enumerate(title_list[:9]):
                line_dict[title] = data_line[index]

            line_dict['samples'] = data_line[9:]
            data_list.append(line_dict)
        f.close()

        # 总位点数
        net_num = len(title_list[9:]) * len(data_list)

        # 缺失率
        random_rate = 0.3

        # random_pos_list = []
        random_pos_list = random.sample(range(0, net_num), int(random_rate * net_num))
        # for i in range(int(random_rate * net_num)):
        #     random_pos = random.randint(0, net_num - 1)
        #     random_pos_list.append(random_pos)
        allele_list = []
        for data in data_list:
            for allele_pair in data['samples']:
                allele_list.append(allele_pair)

        # 产生缺失
        # for index, allele_base in enumerate(allele_list):
        #     if index in random_pos_list:
        #         allele_list[index] = '.|.'

        allele_mat = numpy.mat(allele_list)
        allele_mat[0, random_pos_list] = '.|.'
        allele_list = allele_mat.tolist()[0]

        # temp_list = []
        # for i in range(0, len(allele_list), 2):
        #     temp_list.append('|'.join(allele_list[i:i + 2]))

        new_allele_list = []
        for i in range(0, len(allele_list), len(title_list[9:])):
            new_allele_list.append(allele_list[i:i + len(title_list[9:])])

        for index, data in enumerate(data_list):
            data['samples'] = new_allele_list[index]

        f = open('./new.target.vcf', 'w')

        # 写入注释部分
        f.writelines(header_list)

        new_data_list = []
        for data in data_list:
            one_line = ''
            for v in data.values():

                if isinstance(v, list):
                    for vv in v:
                        one_line += vv + '\t'
                else:
                    one_line += str(v) + '\t'
            one_line += '\n'

            new_data_list.append(one_line)

        f.writelines(new_data_list)
        f.close()


if __name__ == '__main__':

    make_missing_data(file_name='./data_set/123.vcf', isdiv=True)
    print('done')


