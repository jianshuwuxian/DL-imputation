#!/usr/bin/env python
# coding=utf-8
# author: Tianfeng Shi<jianshuwuxian@163.com>
# 2019/11/4 21:16
#
import collections
import numpy as np
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


def get_gene_data(f):
    # VCF文件注释部分
    header_list = []

    line = f.readline()
    while line.startswith("##"):
        header_list.append(line)
        line = f.readline()
    # header_list.append(line)

    line = line.strip('#\n')
    title_list = line.split('\t')

    prev_data_list = []
    data_list = []

    for l in f:
        data_line = l.strip('\n').split('\t')
        # line_dict = collections.OrderedDict()

        # for index, title in enumerate(title_list[:9]):
        #     line_dict[title] = data_line[index]

        # line_dict['samples'] = data_line[9:]
        # data_list.append(line_dict)
        prev_data_list.append(data_line[:9])
        data_list.append(data_line[9:])

    return data_list, header_list, prev_data_list


def read_file(file):
    with open(file, 'r') as f:
        data, header, prev_data = get_gene_data(f)
    return data, header, prev_data


def get_raw_data(file):
    raw_gene_data, header, prev_data = read_file(file=file)
    gene_data = []
    for data in raw_gene_data:
        line = []
        for allele in data:
            final_allele = sum(list(map(int, allele.split("|"))))
            line.append(5 if final_allele >= 5 else final_allele)
        gene_data.append(line)

    return gene_data, header, prev_data


def imputed_file_writer(file_name, header, prev_data, geno_data):
    """
    将填充结果整理写入新文件
    :param file_name:
    :param header:
    :param prev_data:
    :param geno_data:
    :return:
    """
    with open(file_name, "w") as f:
        # 格式转化int到str

        # data = list(map(str, geno_data))

        # 写入文件注释和表头
        f.writelines(header)

        # 写入数据行
        new_data_field = []
        for index in range(0, len(geno_data)):

            data = list(map(str, geno_data[index]))

            new_line = prev_data[index] + data
            new_data_field.append('\t'.join(new_line)+'\n')

        f.writelines(new_data_field)


def multi_process_imputed_file_writer(file_name, prev_data, geno_data):
    with open(file_name, "w") as f:
        # 格式转化int到str

        # data = list(map(str, geno_data))

        # 写入数据行
        new_data_field = []
        for index in range(0, len(geno_data)):

            data = list(map(str, geno_data[index]))

            new_line = prev_data[index] + data
            new_data_field.append('\t'.join(new_line)+'\n')

        f.writelines(new_data_field)


def find_low_index(target, seq, N):
    """
    得到小训练数据
    :param target:
    :param seq:
    :param N:
    :return:
    """
    f_list = seq[np.where(seq < target)[0]].tolist()
    b_list = seq[np.where(seq > target)[0]].tolist()

    pred_list = []

    # # 左段
    # if len(f_list) <= N:
    #     pred_list += f_list
    # else:
    #     pred_list += f_list[(-1)*N:]
    #
    # # 右段
    # if len(b_list) <= N:
    #     pred_list += b_list
    # else:
    #     pred_list += b_list[:N]

    if len(f_list) <= N:
        pred_list += f_list + b_list[:2*N-len(f_list)]
    elif len(b_list) <= N:
        pred_list += f_list[(-1) * (2*N-len(b_list)):] + b_list
    else:
        pred_list += f_list[(-1) * N:] + b_list[:N]

    # # 右段
    # if len(b_list) <= N:
    #     pred_list += b_list
    # else:
    #     pred_list += b_list[:N]

    return pred_list





def generate_small_data(miss_data_array, unmiss_data_array, N):
    seq = []
    label_seq = []
    for i_index, i in enumerate(miss_data_array):
        ind = np.array(list(map(int, i)))

        # 构造miss_list和unmiss_list
        miss_list = np.where(ind == 5)[0]
        unmiss_list = np.where(ind != 5)[0]

        if (2 * N) > len(unmiss_list):
            return None, None
        # 对每一个缺失位点生成训练数据 window_size = 2 * N
        for u_index, u in enumerate(miss_list):
            p_list = find_low_index(u, unmiss_list, N)

            p_seq = list(ind[p_list])
            seq.append(p_seq)

            label = unmiss_data_array[i_index][miss_list[u_index]]
            label_seq.append(label)

    return seq, label_seq


def tt(target, sequ, N):
    seq = np.array(sequ)
    f_list = (np.where(seq < target)[0]).tolist()
    b_list = (np.where(seq > target)[0]).tolist()

    pred_list = []

    # # 左段
    # if len(f_list) <= N:
    #     pred_list += f_list
    # else:
    #     pred_list += f_list[(-1)*N:]
    #
    # # 右段
    # if len(b_list) <= N:
    #     pred_list += b_list
    # else:
    #     pred_list += b_list[:N]

    if len(f_list) <= N:
        pred_list += f_list + b_list[:2 * N - len(f_list)]

    elif len(b_list) <= N:
        pred_list += f_list[(-1) * (2 * N - len(b_list)):] + b_list
    else:
        pred_list += f_list[(-1) * N:] + b_list[:N]

    # # 右段
    # if len(b_list) <= N:
    #     pred_list += b_list
    # else:
    #     pred_list += b_list[:N]

    return pred_list


def generate_train_data_for_each_marker(single_snp, single_snp_index, geno_data_snps, geno_data_inds, N1, N2):
    """
    生成针对单个marker的训练数据和标签
    :param single_snp:
    :param single_snp_index:
    :param geno_data_snps:
    :param geno_data_inds:
    :param N:
    :return:
    """
    miss_indexes = np.where(single_snp == 5)[0]
    unmiss_indexes = np.where(single_snp != 5)[0]

    all_indexes_ind = list(range(0, len(single_snp)))

    train_seq = []

    train_label_seq = []

    test_seq = []
    for index in all_indexes_ind:
        # 生成横向序列:单个SNP
        snp_direction_index = tt(index, all_indexes_ind, N)
        snp_direction_seq = single_snp[snp_direction_index]

        # 生成纵向序列:单个ind
        single_ind = geno_data_inds[index]
        all_indexes_snp = list(range(0, len(single_ind)))

        ind_direction_index = tt(single_snp_index, all_indexes_snp, N)
        ind_direction_seq = single_ind[ind_direction_index]

        # 合并横向纵向序列
        p_seq = np.append(snp_direction_seq, ind_direction_seq)

        if index in unmiss_indexes:
            train_seq.append(list(p_seq))
            train_label_seq.append(single_snp[index])
        elif index in miss_indexes:
            test_seq.append(list(p_seq))


    # for unmiss_index in unmiss_indexes:
    #     # 生成横向序列:单个SNP
    #     snp_direction_index = tt(unmiss_index, all_indexes_ind, N)
    #     snp_direction_seq = single_snp[snp_direction_index]
    #
    #     # 生成纵向序列:单个ind
    #     single_ind = geno_data_inds[unmiss_index]
    #     all_indexes_snp = list(range(0, len(single_ind)))
    #
    #     ind_direction_index = tt(single_snp_index, all_indexes_snp, N)
    #     ind_direction_seq = single_ind[ind_direction_index]
    #
    #     # 合并横向纵向序列
    #     p_seq = np.append(snp_direction_seq, ind_direction_seq)
    #
    #     seq.append(list(p_seq))
    #     label_seq.append(single_snp[unmiss_index])

    return train_seq, train_label_seq, test_seq


def data_generator(geno_data, snp_index, ind_index, N1, N2, ind_num, snp_num, LDMatrix, GMatrix):
    """
    为每个位点生成序列 横纵两个方向信息
    :param geno_data:
    :param snp_index:
    :param ind_index:
    :param N1:
    :param N2:
    :param ind_num:
    :param snp_num:
    :return:
    """
    # 纵向 同一个snp上
    # 计算亲缘关系信息
    #if ind_index - N2 < 0:
    #    l_s = geno_data[snp_index, :ind_index]
    #    l_g = GMatrix[ind_index, :ind_index]

    #    r_s = geno_data[snp_index, ind_index + 1: 2 * N2 + 1]
    #    r_g = GMatrix[ind_index, ind_index+1: ind_index+1+N2+(N2-ind_index)]

    #elif ind_num - 1 - ind_index < N2:
    #    # l_s = geno_data[snp_index, ind_index - N2 - (N2 -(ind_num - 1 - ind_index)): ind_index]
    #    l_s = geno_data[snp_index, ind_num - 2 * N2 - 1: ind_index]
    #    l_g = GMatrix[ind_index, ind_num - 2 * N2 - 1: ind_index]

    #    r_s = geno_data[snp_index, ind_index + 1:]
    #    r_g = GMatrix[ind_index, ind_index + 1:]

    #else:
    #    l_s = geno_data[snp_index, ind_index - N2: ind_index]
    #    l_g = GMatrix[ind_index, ind_index - N2: ind_index]
    #    r_s = geno_data[snp_index, ind_index + 1: ind_index + 1 + N2]
    #    r_g = GMatrix[ind_index, ind_index + 1: ind_index + 1 + N2]

    ## 左序列+右序列+左G序列+右G序列
    #snp_direction_seq = np.append(np.append(np.append(l_s, r_s), l_g), r_g)
    ##snp_direction_seq = np.append(l_s, r_s)
    # 亲缘关系排序版本
    sorted_indexes = np.argsort(GMatrix[ind_index])[-1-(2*N2):-1][::-1]
    seq = geno_data[snp_index, sorted_indexes]
    g_seq = GMatrix[ind_index, sorted_indexes]
    snp_direction_seq = np.append(seq, g_seq)
    
    # 横向 同一个ind上
    # 计算LD信息
    if snp_index - N1 < 0:
        l_s = geno_data[:snp_index, ind_index]
        #l_ld = LDMatrix[:snp_index, snp_index]

        r_s = geno_data[snp_index + 1: 2 * N1 + 1, ind_index]
        #r_ld = LDMatrix[snp_index, snp_index + 1: 2 * N1 + 1]

    elif snp_num - 1 - snp_index < N1:
        # l_s = geno_data[snp_index - N1 - (N1 -(snp_num - 1 - snp_index)): snp_index, ind_index]
        l_s = geno_data[snp_num - 2 * N1 - 1: snp_index, ind_index]
        #l_ld = LDMatrix[snp_num - 2 * N1 - 1: snp_index, snp_index]

        r_s = geno_data[snp_index + 1:, ind_index]
        #r_ld = LDMatrix[snp_index, snp_index + 1:]

    else:
        l_s = geno_data[snp_index - N1: snp_index, ind_index]
        #l_ld = LDMatrix[snp_index - N1: snp_index, snp_index]

        r_s = geno_data[snp_index + 1: snp_index + 1 + N1, ind_index]
        #r_ld = LDMatrix[snp_index, snp_index + 1: snp_index + 1 + N1]

    # 左序列+右序列+左LD序列+右LD序列
    #ind_direction_seq = np.append(np.append(np.append(l_s, r_s), l_ld), r_ld)
    ind_direction_seq = np.append(l_s, r_s)
    # 合并横向纵向序列
    p_seq = np.append(ind_direction_seq, snp_direction_seq)

    return p_seq


def data_generator_using_matrix(geno_data, snp_index, ind_index, N1, N2, ind_num, snp_num):
    """
    为每个位点生成序列 以该位点为中心的矩阵
    :param geno_data:
    :param snp_index:
    :param ind_index:
    :param N1:
    :param N2:
    :param ind_num:
    :param snp_num:
    :return:
    """
    if snp_index - N1 < 0:
        if ind_index - N2 < 0:

            up = geno_data[:snp_index, :2*N2+1]

            mid_l = geno_data[snp_index, :ind_index]

            mid_r = geno_data[snp_index, ind_index + 1: 2 * N2 + 1]

            down = geno_data[snp_index+1: snp_index+1 + N1 + N1-snp_index, :2*N2+1]

        elif ind_num - 1 - ind_index < N2:
            up = geno_data[:snp_index, ind_num - 2*N2 - 1:]

            mid_l = geno_data[snp_index, ind_index - N2 - (N2 - (ind_num-1-ind_index)):ind_index]
            mid_r = geno_data[snp_index, ind_index+1:]

            down = geno_data[snp_index+1:snp_index+1+N1+N1-snp_index, ind_num - 2*N2 - 1:]

        else:
            up = geno_data[:snp_index, ind_index-N2:ind_index+N2+1]

            mid_l = geno_data[snp_index, ind_index - N2:ind_index]

            mid_r = geno_data[snp_index, ind_index+1: ind_index+1+N2]

            down = geno_data[snp_index+1:snp_index+1+N1+N1-snp_index, ind_index-N2:ind_index+N2+1]

    elif snp_num - 1 - snp_index < N1:
        if ind_index - N2 < 0:
            up = geno_data[snp_index - N1 - (N1 - (snp_num-1-snp_index)):snp_index, :2*N2+1]

            mid_l = geno_data[snp_index, :ind_index]

            mid_r = geno_data[snp_index, ind_index + 1: 2 * N2 + 1]

            down = geno_data[snp_index+1:, :2*N2+1]
        elif ind_num - 1 - ind_index < N2:
            up = geno_data[snp_index - N1 - (N1 - (snp_num-1-snp_index)):snp_index, ind_num - 2*N2 - 1:]

            mid_l = geno_data[snp_index, ind_index - N2 - (N2 - (ind_num-1-ind_index)):ind_index]

            mid_r = geno_data[snp_index, ind_index+1:]

            down = geno_data[snp_index+1:, ind_num - 2*N2 - 1:]
        else:
            up = geno_data[snp_index - N1 - (N1 - (snp_num-1-snp_index)):snp_index, ind_index-N2:ind_index+N2+1]

            mid_l = geno_data[snp_index, ind_index - N2:ind_index]

            mid_r = geno_data[snp_index, ind_index+1: ind_index+1+N2]

            down = geno_data[snp_index+1:, ind_index-N2:ind_index+N2+1]
    else:
        if ind_index - N2 < 0:
            up = geno_data[snp_index-N1:snp_index, :2*N2+1]

            mid_l = geno_data[snp_index, :ind_index]

            mid_r = geno_data[snp_index, ind_index + 1: 2 * N2 + 1]

            down = geno_data[snp_index+1:snp_index+1+N1, :2*N2+1]
        elif ind_num - 1 - ind_index < N2:
            up = geno_data[snp_index-N1:snp_index, ind_num - 2*N2 - 1:]

            mid_l = geno_data[snp_index, ind_index - N2 - (N2 - (ind_num-1-ind_index)):ind_index]

            mid_r = geno_data[snp_index, ind_index+1:]

            down = geno_data[snp_index+1:snp_index+1+N1, ind_num - 2*N2 - 1:]
        else:
            up = geno_data[snp_index-N1:snp_index, ind_index-N2:ind_index+N2+1]

            mid_l = geno_data[snp_index, ind_index - N2:ind_index]

            mid_r = geno_data[snp_index, ind_index+1:ind_index+1+N2]

            down = geno_data[snp_index+1:snp_index+1+N1, ind_index-N2:ind_index+N2+1]

    p_seq = np.append(np.append(np.append(up, mid_l),mid_r), down)

    return p_seq


def getLDMatrix(ld_file_name, snp_num=20):
    """
    返回LD矩阵
    :param file_name:
    :param snp_num:
    :return:
    """
    f = open(ld_file_name, 'r')

    lines = f.readlines()
    f.close()

    lines = lines[1:]

    first_line = lines[0].strip('\n').split()
    current = first_line[2]

    same_index = 0
    finished_snp_num = 0
    snp_list = []
    snp = np.zeros(snp_num)
    for line in lines:

        data = line.strip('\n').split()
        cur_snp = data[2]
        if current != cur_snp:
            finished_snp_num += 1
            snp_list.append(snp)
            current = cur_snp
            same_index = finished_snp_num
            snp = np.zeros(snp_num)
            snp[same_index + 1] = data[-1]
            same_index += 1
        else:
            snp[same_index+1] = data[-1]
            same_index += 1
    snp_list.append(snp)
    snp_list.append(np.zeros(snp_num))

    snp_list = np.array(snp_list)

    return snp_list


def getGMatrix(gmatrix_file_name):
    """
    返回亲缘关系矩阵
    :param gmatrix_file_name:
    :param ind_num:
    :return:
    """
    f = open(gmatrix_file_name, 'r')

    lines = f.readlines()
    f.close()

    ind_list = []
    for line in lines:

        values = np.array(line.strip('\n').split())
        ind_list.append(values)

    ind_list = np.array(ind_list)

    return ind_list


if __name__ == '__main__':
    getLDMatrix('./processed_data7/plink.ld', snp_num=4753)
    # getGMatrix(gmatrix_file_name='./processed_data7/result.txt')
