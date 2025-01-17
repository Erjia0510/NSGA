"""
文件目的：计算适应值

注：

@author: Jia
"""
import copy

import numpy as np

import globalv

import distance






class calculate_tcost(object):
    def __init__(self):
        self.UtoT, self.TtoT = distance.normal()
        # self.fitness2 = np.array([])



  ### 计算航迹代价
    def cost(self, data):
        disU = [[] for k in range(globalv.U_num)]  # 存储每个无人机飞过的每段距离
        timeU = [[] for k in range(globalv.U_num)]  # 存储每个无人机飞行所用时间
        time_task = [[] for k in range(globalv.T_num)]  # 存储每个任务的开始时间以及结束时间
        # fitness1 = np.array([])  # 距离
        # fitness2 = np.array([])  # 时间
        UandT = [[] for _ in range(globalv.U_num)]
        U_S = copy.deepcopy(globalv.U_S)

        for i in range(data.shape[1]):
            u = data[3][i]
            t = data[1][i]
            type = data[2][i]
            t_re = data[4][i]
            if len(UandT[u-1]) == 0:
                dis = self.UtoT[t-1]
            else:
                if(type == 2 and U_S[u-1] < t_re):  # 当资源不足时需要回仓库补充资源
                    dis = self.UtoT[t-1] + self.UtoT[UandT[u-1][-1] - 1]
                    timeU[u-1].append((globalv.U_S[u-1] - t_re) * globalv.U_T)
                else:
                    dis = self.TtoT[UandT[u - 1][-1] - 1][t - 1]
            if type == 2:
                U_S[u - 1] -= t_re
            UandT[u - 1].append(t)
            time_load = dis / globalv.U_V[u-1] + (t_re * globalv.T_T) if type == 2 else 30   # 飞行距离+任务时间
            timeU[u-1].append(time_load)
            disU[u-1].append(dis)
            # # 加入任务时序
            # if len(time_task[t-1]) == 0:
            #     time_task[t - 1].extend([time_load, time_load + globalv.time_T[type - 1]])
            # else:
            #     time1 = time_task[t-1][-1]
            #     time2 = sum(timeU[u-1])
            #     if time2 >= time1:     # 不需要等待
            #         time_task[t-1].extend([time2, time2 + globalv.time_T[type-1]])
            #     else:                  # 需要等待
            #         time_wait = time1 - time2
            #         timeU[u-1].append(time_wait)
            #         dis_wait = time_wait * globalv.U_V[u-1]
            #         disU[u-1].append(dis_wait)
            #         time_task[t - 1].extend([time1, time1 + globalv.time_T[type - 1]])

        # 算回仓库的距离
        for i in range(len(UandT)):
            if not UandT[i]:
                continue
            dis = self.UtoT[[UandT[i][-1]- 1]]
            disU[i].append(dis)
            timeU[i].append(dis / globalv.U_V[i])
            if sum(disU[i]) > globalv.max_flight:
                # fitness1 = np.append(fitness1, float('inf'))
                fitness1 = float('inf')
                fitness2 = float('inf')
            else:
                fitness1 = sum(sum(sublist) for sublist in disU)
                fitness2 = max([sum(sublist) for sublist in timeU])

        # fitness1 = np.append(fitness1, float('inf'))
        # cost = globalv.w1 * sum(sum(sublist) for sublist in disU) + globalv.w2 * max([sum(sublist) for sublist in disU])
        # fitness1 = np.append(fitness1, cost)


        return fitness1, fitness2














