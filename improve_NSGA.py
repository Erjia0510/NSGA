import math

import globalv
import copy
import random
from operator import itemgetter
import numpy as np
from fitness import calculate_tcost
from Chromo import Chromo
from K_means import k_means


def crowd_distance(pop, front):
    lenth = len(front)
    distance = [0 for i in range(lenth)]
    distance[0] = float('inf')
    distance[-1] = float('inf')
    if lenth > 2:
        fit1_values = [pop[individual]['fit1'] for individual in front]
        fit2_values = [pop[individual]['fit2'] for individual in front]
    for k in range(1, lenth - 1):
        # 拥挤度归一化
        distance[k] = distance[k] + abs(fit1_values[k + 1] - fit1_values[k - 1]) / (
                    max(fit1_values) - min(fit1_values)) + abs(fit2_values[k + 1] - fit2_values[k - 1]) / (
                                  max(fit2_values) - min(fit2_values))

    for k in range(len(front)):
        pop[front[k]]['distance'] = distance[k]
    return pop


def binary_tournament(rank, pop):
    if rank[0] != rank[1]:  # 如果两个个体有支配关系，即在两个不同的rank中，选择rank小的
        return pop[0] if rank[0] < rank[1] else pop[1]
    elif pop[0]['distance'] != pop[1]['distance']:  # 如果两个个体rank相同，比较拥挤度距离，选择拥挤读距离大的
        return pop[0] if pop[0]['distance'] > pop[1]['distance'] else pop[1]
    else:  # 如果rank和拥挤度都相同，返回任意一个都可以
        return pop[0]





class improve_NSGA(object):
    def __init__(self):
        self.population_size = globalv.population_size
        self.cross_p = globalv.cross_num  # 交叉概率
        self.mutation_p = globalv.mutation  # 变异概率
        self.max_step = globalv.time_max
        # 存储时间与与之对应的适应值
        self.time_mat = []
        self.fitness_mat = []
        a, b = k_means(globalv.Tar_start, len(globalv.A))
        self.cluster_U = [[x] for x in globalv.A]
        for i in range(len(globalv.S)):
            self.cluster_U[i % (len(globalv.A))].append(globalv.S[i])
        self.cluster_T = a
        pop = []
        ctc = calculate_tcost()
        # 初始化染色体数组
        for j in range(self.population_size):
            gene = np.zeros((6, globalv.T_type * globalv.T_num)).astype((int))
            # gene = gene.astype(int)
            # 第一行：生成自然数
            gene[0, :] = np.arange(1, gene.shape[1] + 1)
            # 第二行：任务重复三次
            temp = np.repeat(np.arange(1, globalv.T_num + 1), globalv.T_type)
            np.random.shuffle(temp)
            gene[1, :] = temp
            # 第三行：转换为目标序染色体，将任务类型C A V依次加入
            sorted_id = np.lexsort((gene[0, :], gene[1, :]))
            gene = gene[:, sorted_id]
            gene[2, :] = np.tile(np.arange(1, globalv.T_type + 1), globalv.T_num)
            # 第四行：根据第三行中的所需任务类型，从对应的无人机集合US/UA中随机选择一个有能力的无人机作为第四行中的分配无人机
            for i in range(len(self.cluster_T)):
                for ind in self.cluster_T[i]:
                    if random.random() < 0.7:
                        gene[3, (ind - 1) * 3] = random.choice(self.cluster_U[i][1:])
                        gene[3, (ind - 1) * 3 + 1] = self.cluster_U[i][0]
                        gene[3, (ind - 1) * 3 + 2] = random.choice(self.cluster_U[i][1:])
                    else:
                        gene[3, (ind - 1) * 3] = globalv.S[random.randint(0, len(globalv.S) - 1)]
                        gene[3, (ind - 1) * 3 + 1] = globalv.A[random.randint(0, len(globalv.A) - 1)]
                        gene[3, (ind - 1) * 3 + 2] = globalv.S[random.randint(0, len(globalv.S) - 1)]
                # 第五行：所需资源
                    gene[4, (ind - 1) * 3] = 0
                    gene[4, (ind - 1) * 3 + 2] = 0
                    gene[4, (ind - 1) * 3 + 1] = globalv.Tar_need[i - 1]
            # 转换回去
            # 第一行：自然数
            sorted_id = np.argsort(gene[0, :])
            gene = gene[:, sorted_id]
            fitness1, fitness2 = ctc.cost(gene)
            pop.append({'chromo': Chromo(data=gene), 'fit1': fitness1, 'fit2': fitness2})
        self.pop = pop
        self.first_pop = pop
        # self.bestindividual = self.selectBest(self.pop)  # 保存最好的个体数据{'Gene':Gene(), 'fitness':fitness}

    def dominates(self, fit11, fit12, fit21, fit22):
        value = False
        if (fit11 < fit21) or fit12 < fit22:
            if fit11 <= fit21 and fit12 <= fit22:
                value = True
        return value

    def fast_nondominated(self, pop):
        front = [[]]
        dominated_solutions = [[] for _ in range(len(pop))]
        dominated_count = [0 for _ in range(len(pop))]
        rank = [0 for _ in range(len(pop))]
        for ind in range(len(pop)):
            for p in range(len(pop)):
                if self.dominates(pop[ind]['fit1'], pop[ind]['fit2'], pop[p]['fit1'], pop[p]['fit2']):
                    dominated_solutions[ind].append(p)
                elif self.dominates(pop[p]['fit1'], pop[p]['fit2'], pop[ind]['fit1'], pop[ind]['fit2']):
                    dominated_count[ind] += 1
            # 如果一个解没有被任何其他解支配，则将其归为 Pareto 前沿的第一层
            if dominated_count[ind] == 0:
                rank[ind] = 0
                if ind not in front[0]:
                    front[0].append(ind)

        i = 0
        while front[i]:
            Q = []
            for p in front[i]:
                for q in dominated_solutions[p]:
                    dominated_count[q] -= 1
                    if dominated_count[q] == 0:
                        rank[q] = i + 1
                        if q not in Q:
                            Q.append(q)
            i += 1
            front.append(Q)
        del front[-1]
        return front, rank, dominated_solutions

    # 交叉函数：将选定的两条亲本染色体变换为目标序亲本染色体，随机选择两条交叉位点。
    def crossover(self, offspring):
        dim = offspring[0]['chromo'].data.shape[1]  # 获得数据维数，即基因位数

        chromoinfo1 = offspring[0]['chromo'].data  # 交叉的第一个数据
        chromoinfo2 = offspring[1]['chromo'].data  # 交叉的第二个数据

        # 将选定的两条亲本染色体变换为目标序亲本染色体
        sorted_id1 = np.lexsort((chromoinfo1[2, :], chromoinfo1[1, :]))
        sorted_id2 = np.lexsort((chromoinfo2[2, :], chromoinfo2[1, :]))

        # 交换任务点
        if dim == 1:
            pos1 = 1
            pos2 = 1
        else:
            pos1 = random.randrange(dim)
            pos2 = random.randrange(dim)
        newchromo1 = Chromo(data=[])
        newchromo2 = Chromo(data=[])

        # newchromo1改变 chromeinfo不改变

        newchromo1_temp = copy.deepcopy(chromoinfo1)
        newchromo2_temp = copy.deepcopy(chromoinfo2)

        for i in range(dim):
            if min(pos1, pos2) <= i < max(pos1, pos2):  # 交换的部分维度

                temp = newchromo1_temp[3][sorted_id1[i]]
                newchromo1_temp[3][sorted_id1[i]] = newchromo2_temp[3][sorted_id2[i]]
                newchromo2_temp[3][sorted_id2[i]] = temp

        newchromo1.data = newchromo1_temp
        newchromo2.data = newchromo2_temp

        return newchromo1, newchromo2

        # 变异:

    def mutate(self, crossoff):

        dim = crossoff.data.shape[1]
        # 选择单点变异的点
        if dim == 1:
            pos = 0
        else:
            pos = random.randrange(0, dim)

        temp1 = copy.deepcopy(globalv.S)
        temp2 = copy.deepcopy(globalv.A)

        if crossoff.data[2][pos] == 1 or crossoff.data[2][pos] == 3:
            if len(globalv.S) != 1:
                # temp1.remove(crossoff.data[3][pos])
                crossoff.data[3][pos] = temp1[random.randint(0, len(temp1) - 1)]
        else:
            if len(globalv.A) != 1:
                # temp2.remove(crossoff.data[3][pos])
                crossoff.data[3][pos] = temp2[random.randint(0, len(temp2) - 1)]
        return crossoff


    def evolve(self):
        print('Start of evolution')
        ctc = calculate_tcost()

        for now_step in range(self.max_step):
            print("############ Generation {} ##############".format(now_step))
            nextpop = []
            best_pop = []
            front, rank, dominated_solutions = self.fast_nondominated(self.pop)
            # crowd = []
            for i in range(len(front)):
                self.pop = crowd_distance(self.pop, front[i][:])
                # crowd.append(distance)
            next_p = self.pop[:]
            while len(next_p) != 2*globalv.population_size:
                idx = np.random.choice(np.arange(len(self.pop)), size=4, replace=False)
                parentpop1 = binary_tournament([rank[idx[0]], rank[idx[1]]], [self.pop[idx[0]], self.pop[idx[1]]])
                parentpop2 = binary_tournament([rank[idx[2]], rank[idx[3]]], [self.pop[idx[2]], self.pop[idx[3]]])
                if random.random() < self.cross_p:  # 交叉
                    crossoff1, crossoff2 = self.crossover([parentpop1, parentpop2])
                else:
                    crossoff1 = parentpop1['chromo']
                    crossoff2 = parentpop2['chromo']
                    # nextoff.extend(offspring)  # 直接追加两个后代
                if random.random() < self.mutation_p:  # 变异
                    muteoff1 = self.mutate(crossoff1)
                    muteoff2 = self.mutate(crossoff2)
                    fit_muteoff1 = ctc.cost(muteoff1.data)
                    fit_muteoff2 = ctc.cost(muteoff2.data)
                    next_p.append({'chromo': muteoff1, 'fit1': fit_muteoff1[0], 'fit2': fit_muteoff1[1]})
                    next_p.append({'chromo': muteoff2, 'fit1': fit_muteoff2[0], 'fit2': fit_muteoff2[1]})
                else:
                    fit_muteoff1 = ctc.cost(crossoff1.data)
                    fit_muteoff2 = ctc.cost(crossoff2.data)
                    next_p.append({'chromo': crossoff1, 'fit1': fit_muteoff1[0], 'fit2': fit_muteoff1[1]})
                    next_p.append({'chromo': crossoff2, 'fit1': fit_muteoff2[0], 'fit2': fit_muteoff2[1]})
            front, rank, dominated_solutions = self.fast_nondominated(next_p)
            for i in range(len(front)):
                self.pop = crowd_distance(next_p, front[i][:])

            i = 0
            while len(nextpop) + len(front[i]) <= globalv.population_size:
                for k in range(len(front[i])):
                    nextpop.append(next_p[front[i][k]])
                    # if i == 0:
                    #     best_pop.append(next_p[front[i][k]])
                i = i + 1
            if len(nextpop)<globalv.population_size:
                sorted_id = sorted([next_p[k] for k in front[i]],key = itemgetter("distance"), reverse=True)
                k = 0
                while len(nextpop) != globalv.population_size:
                    nextpop.append(sorted_id[k])
                    k = k + 1
            self.pop = nextpop
            best_pop = [next_p[k] for k in front[0]]
            # print("Best individual fit1 is {}".format([best_pop[i]['fit1'] for i in range(len(best_pop))]))
            # print("Best individual fit2 is {}".format([best_pop[i]['fit2'] for i in range(len(best_pop))]))
            # print("  Min fitness of current pop: {}".format(min(fits)))
        return best_pop, self.first_pop











