#!/usr/bin/env python
# coding=utf-8
# author: Tianfeng Shi<jianshuwuxian@163.com>
# 2020/2/5 22:48
#
from data_new import get_raw_data, imputed_file_writer, data_generator, data_generator_using_matrix, getLDMatrix, getGMatrix
import numpy as np
import time
from miss_generator import make_missing_data
import pickle
import os
import sys
import multiprocessing as mp
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


seq_size = 100000
num_epochs = 50
train_rate = 0.8
val_rate = 0.1
test_rate = 0.1

# window_size = 2 * N
# N = 10
#N_cross_snp = 10  # the number of different snp alleles besides the loci
#N_cross_ind = 3  # the number of different ind alleles besides the loci
# N = 5


def matrixToT(data):
    data_array = np.array(data)
    # 元素是每个ind
    return data_array.T


def main_func(data_dir, geno_file_num, matrix_info=False, N_cross_snp=0, N_cross_ind=0):
    """
    数据预处理
    """

    # 带有缺失的geno文件
    file_name = data_dir + "sub_geno_file" + str(geno_file_num) + ".txt"

    ld_file = data_dir + "plink.ld"

    gmatrix_file = data_dir + "sub_gmatrix" + str(geno_file_num) + ".txt"

    # 读取数据
    with timer(name='get geno data', level='important'):
        # raw_geno_data_with_miss, header_data, prev_data = get_raw_data(file=file_name)  # 得到0125序列 5代表缺失
        with open(file_name, 'r') as f:
            data_list = []

            for l in f:
                data_line = l.strip('\n').split('\t')
                data_list.append(data_line[9:])

        raw_geno_data_with_miss = []
        for data in data_list:
            line = []
            for allele in data:
                final_allele = sum(list(map(int, allele.split("|"))))
                line.append(5 if final_allele >= 5 else final_allele)
            raw_geno_data_with_miss.append(line)

    # # 矩阵转置
    geno_data_with_miss_array_inds = matrixToT(raw_geno_data_with_miss)  # 得到 num_ind行 num_snp列
    geno_data_with_miss_array_snps = np.array(raw_geno_data_with_miss)

    snp_num = len(geno_data_with_miss_array_snps)
    ind_num = len(geno_data_with_miss_array_inds)

    train_seq = []
    train_label_seq = []
    test_seq = []
    
    LDMatrix = []
    
    ## 计算LD矩阵
    #with timer(name='计算LD矩阵', level='important'):
    #    LDMatrix = getLDMatrix(ld_file, snp_num=snp_num)

    # 计算亲缘关系矩阵
    with timer(name='计算亲缘关系矩阵', level='important'):
        GMatrix = getGMatrix(gmatrix_file)

    # 填充好的数据
    imputed_data = []

    if ((2 * N_cross_snp) <= snp_num - 1) and ((2 * N_cross_ind) <= ind_num - 1):
        """
        新 数据预处理
        """
        # 为每个marker生成各自的训练数据和标签
        with timer(name='generate data', level='important'):
            for snp_index in range(0, snp_num):
                snp_train_seq = []
                snp_train_label_seq = []
                snp_test_seq = []
                for ind_index in range(0, ind_num):
                    if not matrix_info:  # 判断是否使用矩阵形式的信息
                        # 利用横纵两个方向的信息
                        p_seq = data_generator(geno_data_with_miss_array_snps, snp_index, ind_index, N_cross_snp, N_cross_ind, ind_num, snp_num, LDMatrix, GMatrix)

                    else:
                        # 利用矩阵形式的信息
                        p_seq = data_generator_using_matrix(geno_data_with_miss_array_snps, snp_index, ind_index, N_cross_snp, N_cross_ind, ind_num, snp_num)

                    if geno_data_with_miss_array_snps[snp_index, ind_index] == 5:
                        snp_test_seq.append(list(p_seq))
                    else:
                        snp_train_seq.append(list(p_seq))
                        snp_train_label_seq.append(geno_data_with_miss_array_snps[snp_index, ind_index])
                train_seq.append(snp_train_seq)
                train_label_seq.append(snp_train_label_seq)
                test_seq.append(snp_test_seq)

                if len(train_seq) > 100:
                    with open(data_dir + "sub_train_seq" + str(geno_file_num) + ".txt", "a") as f1:
                        lines = []
                        for t_s in train_seq:
                            line = ""
                            for item in t_s:
                                line += ' '.join(list(map(str, item))) + '\t'
                            line = line.strip("\t")
                            lines.append(line + '\n')
                        f1.writelines(lines)

                    with open(data_dir + "sub_train_label_seq" + str(geno_file_num) + ".txt", "a") as f2:
                        lines = []
                        for t_l_s in train_label_seq:
                            lines.append('\t'.join(list(map(str, t_l_s))) + '\n')
                        f2.writelines(lines)


                    with open(data_dir + "sub_test_seq" + str(geno_file_num) + ".txt", "a") as f3:
                        lines = []
                        for t_s in test_seq:
                            line = ""
                            for item in t_s:
                                line += ' '.join(list(map(str, item))) + '\t'
                            line = line.strip("\t")
                            lines.append(line + '\n')
                        f3.writelines(lines)

                    train_seq = []
                    train_label_seq = []
                    test_seq = []




        # # 写入文件 pickle版本
        # with timer(name='write train data to files', level='important'):
        #
        #     f1 = open(data_dir + "sub_train_seq" + str(geno_file_num) + ".txt", "wb")
        #     pickle.dump(train_seq, f1)
        #     f1.close()
        #
        #     f2 = open(data_dir + "sub_train_label_seq" + str(geno_file_num) + ".txt", "wb")
        #     pickle.dump(train_label_seq, f2)
        #     f2.close()
        #
        #     f3 = open(data_dir + "sub_test_seq" + str(geno_file_num) + ".txt", "wb")
        #     pickle.dump(test_seq, f3)
        #     f3.close()

        # 写入文件 手动写入
        with timer(name='write train data to files', level='important'):
            # train_seq文件
            with open(data_dir + "sub_train_seq" + str(geno_file_num) + ".txt", "a") as f1:
                lines = []
                for t_s in train_seq:
                    line = ""
                    for item in t_s:
                        line += ' '.join(list(map(str, item))) + '\t'
                    line = line.strip("\t")
                    lines.append(line + '\n')
                f1.writelines(lines)

            # train_label_seq文件
            with open(data_dir + "sub_train_label_seq" + str(geno_file_num) + ".txt", "a") as f2:
                lines = []
                for t_l_s in train_label_seq:
                    lines.append('\t'.join(list(map(str, t_l_s))) + '\n')
                f2.writelines(lines)

            # test_seq文件
            with open(data_dir + "sub_test_seq" + str(geno_file_num) + ".txt", "a") as f3:
                lines = []
                for t_s in test_seq:
                    line = ""
                    for item in t_s:
                        line += ' '.join(list(map(str, item))) + '\t'
                    line = line.strip("\t")
                    lines.append(line + '\n')
                f3.writelines(lines)

    else:
        print("请设置合适的N值")

    return


if __name__ == '__main__':
    data_dir = sys.argv[1]
    geno_file_num = sys.argv[2]
    # data_dir = './processed_data9/'
    N_c_s = int(sys.argv[3])
    N_c_i = int(sys.argv[4])
    main_func(data_dir=data_dir, geno_file_num=geno_file_num, matrix_info=False, N_cross_snp=N_c_s, N_cross_ind=N_c_i)












