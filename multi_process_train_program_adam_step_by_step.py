#!/usr/bin/env python
# coding=utf-8
# author: Tianfeng Shi<jianshuwuxian@163.com>
# 2020/2/13 22:58
#

from data import get_raw_data, imputed_file_writer, data_generator, data_generator_using_matrix, getLDMatrix, getGMatrix, multi_process_imputed_file_writer
from keras import models
from keras import layers
import keras
import numpy as np
import time
import sys
from keras.preprocessing.text import Tokenizer
# from keras.utils import to_categorical
from miss_generator import make_missing_data
import matplotlib.pyplot as plt
import pickle
import os
from sklearn.metrics import classification_report
from keras.utils.np_utils import to_categorical
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
#N_cross_snp = 5  # the number of different snp alleles besides the loci
#N_cross_ind = 0  # the number of different ind alleles besides the loci


def build_model(lr):
    """
    构建网络
    :return:
    """
    model = models.Sequential()
    model.add(layers.Bidirectional(layers.LSTM(64, return_sequences=True)))
    #model.add(layers.Bidirectional(layers.LSTM(64, return_sequences=True)))
    # 分成0 1 2三类 多分类的激活函数使用softmax
    model.add(layers.Dense(3, activation='softmax'))

    # 多分类的损失函数使用categorical_crossentropy
    adam = keras.optimizers.Adam(learning_rate=lr)
    model.compile(optimizer=adam, loss='categorical_crossentropy', metrics=['accuracy'])
    return model


def load_data(train_data, train_label, test_data):
    m_d = np.reshape(train_data, (len(train_data), 1, len(train_data[0])))

    d = np.reshape(train_label, (len(train_label), 1, 1))
    x = np.array(test_data)
    test_d = np.reshape(x, (len(test_data), 1, len(test_data[0])))
    return m_d, d, test_d


def matrixToT(data):
    data_array = np.array(data)
    # 元素是每个ind
    return data_array.T


def train_test_file_decoder(single_line):
    samples = []
    s = single_line.strip('\n').split('\n')
    for item in s:
        samples += item.split('\t')
    #samples = single_line.strip('\n').split('\t')

    results = []
    for sample in samples:
        #results.append(list(map(float, sample.split())))
        results.append(sample.split())
    return results

def test_file_decoder(single_line):
    samples = []
    s = single_line.strip('\n').split(':')
    for item in s:
        temp = item.strip('\n').split('\t')
        if temp[0]:
            samples += temp

    results = []
    for sample in samples:
        results.append(sample.split())
    return results


def train_label_decoder(single_line):
    results = []
    s = single_line.strip('\n').split('\n')
    for item in s:
        results += item.split('\t')
    #results = single_line.strip('\n').split('\t')

    return results


def main_func(snp_num, data_dir, train_file_num, learning_rate):
    """
    数据预处理
    """

    # 读取数据
    raw_file_name = data_dir + "sub_geno_file" + str(train_file_num) + ".txt"

    with timer(name='get geno data', level='important'):
        # raw_geno_data_with_miss, header_data, prev_data = get_raw_data(file=file_name)  # 得到0125序列 5代表缺失
        with open(raw_file_name, 'r') as f:
            data_list = []
            prev_data_list = []
            for l in f:
                data_line = l.strip('\n').split('\t')
                prev_data_list.append(data_line[:9])
                data_list.append(data_line[9:])

        raw_geno_data_with_miss = []
        for data in data_list:
            line = []
            for allele in data:
                final_allele = sum(list(map(int, allele.split("|"))))
                line.append(5 if final_allele >= 5 else final_allele)
            raw_geno_data_with_miss.append(line)

    geno_data_with_miss_array_snps = np.array(raw_geno_data_with_miss)

    # 填充好的数据
    imputed_data = []

    # 训练数据文件
    train_data_file_name = data_dir + "sub_train_file" + str(train_file_num) + ".txt"

    # 训练标签文件
    train_label_file_name = data_dir + "sub_train_label_file" + str(train_file_num) + ".txt"

    # 测试数据文件
    test_data_file_name = data_dir + "sub_test_file" + str(train_file_num) + ".txt"

    '''
    读取训练数据
    '''
    with open(train_data_file_name, 'r') as f:
        train_seq = f.readlines()

    with open(train_label_file_name, 'r') as f:
        train_label_seq = f.readlines()

    with open(test_data_file_name, 'r') as f:
        test_seq = f.readlines()

    '''
    构建并训练模型
    '''
    # 针对每个marker
    for snp in range(0, int(snp_num)):
        print("----------------------------------------------------------------")
        print("MODEL", str(snp + 1), ":")
        with timer():
            '''
            读取该maker的数据
            '''
            # print("DATA LOADING......")
            #with timer("读取该marker的数据"):
            train_data_list = train_test_file_decoder(train_seq[snp])
            train_label_list = train_label_decoder(train_label_seq[snp])
            test_data_list = test_file_decoder(test_seq[snp])

            '''
            预处理成张量
            '''
            # print("DATA PREPROCESSING......")
            #with timer("预处理成张量"):
            train_data, train_labels, test_data = load_data(train_data_list, train_label_list, test_data_list)

            # one-hot编码
            train_labels = to_categorical(train_labels, num_classes=3)

            '''
            构建模型
            '''
            # print("MODEL BUILDING......")
            #with timer("构建模型"):
            model = build_model(learning_rate)

            '''
            训练模型
            '''
            #with timer("训练模型", "important"):
            history = model.fit(train_data, train_labels, epochs=num_epochs, batch_size=128, validation_split=0.15,
                                    verbose=False)

            '''
            输入缺失位点数据到模型，生成预测
            '''
            #with timer("生成预测"):
            pred = model.predict(test_data)

            pred_list = [np.argmax(item) for item in pred]
            
            '''
            将预测结果加入到原数据中
            '''
            raw_seq = geno_data_with_miss_array_snps[snp]
            
            raw_seq[np.where(raw_seq == 5)[0]] = pred_list

            imputed_seq = raw_seq

            imputed_data.append(imputed_seq)
    '''
    写入新文件
    '''
    #with timer(name='writing result file', level='important'):
    new_file_name = data_dir+"sub_result" + str(train_file_num) + ".txt"

    multi_process_imputed_file_writer(new_file_name, prev_data_list, imputed_data)

    return


if __name__ == '__main__':
    snp_num = sys.argv[1]
    data_dir = sys.argv[2]
    train_file_num = sys.argv[3]
    lr = float(sys.argv[4])
    print('This_is_file ' + train_file_num)
    main_func(snp_num=snp_num, data_dir=data_dir, train_file_num=train_file_num, learning_rate=lr)











