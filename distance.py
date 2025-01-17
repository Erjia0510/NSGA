"""
文件目的：所有与距离有关的函数


@author: Jia
"""

import globalv
import numpy as np
import math



# 加入K-means
# data = globalv.Tar_start
# a, b = k_means(data, globalv.U_num)
# globalv.UAV_start = b
# UAV_start = b

target_pos = globalv.Tar_start

row = globalv.U_num
column = globalv.T_num

# 无人机到各目标点的距离
def dis_UT():
    UAV_start = np.zeros((globalv.U_num, 2))
    UtoT = np.zeros(column)  #初始化矩阵信息
    for j in range(column):
        UtoT[j] = np.abs(target_pos[j][1]) + np.abs(target_pos[j][0])

    return UtoT  # 函数返回距离矩阵


# 各目标点之间的距离
def dis_TT():
    TtoT = np.zeros((column, column))

    for i in range(column):

        for j in range(column):
            d_x = target_pos[i][0] - target_pos[j][0]

            d_y = target_pos[i][1] - target_pos[j][1]

            TtoT[i][j] = np.abs(d_x) + np.abs(d_y)

    return TtoT


# 距离标准化
def normal():

    UtoT = dis_UT()

    TtoT = dis_TT()

    return UtoT, TtoT
