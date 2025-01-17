import numpy as np
import random
import globalv
import matplotlib.pyplot as plt

def cal_distance(node, centor):
    return np.sqrt(np.sum(np.square(node - centor)))


def random_centor(data, k):
    data = list(data)
    return random.sample(data, k)


def random_centor1(data, k):
    n = len(data[0])  # n维
    centor = np.array([[0] * n for _ in range(k)])  # 一定要将列表转换为数组
    for j in range(n):
        min_j = np.min(data[:, j])
        max_j = np.max(data[:, j])
        centor[:, j] = np.random.rand(k) * (max_j - min_j) + min_j
    return centor


def get_cluster(data, centor):
    cluster_dict = dict()
    k = len(centor)
    for node in data:
        cluster_class = -1
        min_distance = float('inf')
        for i in range(k):
            dist = cal_distance(node, centor[i])
            if dist < min_distance:
                cluster_class = i
                min_distance = dist
        if cluster_class not in cluster_dict.keys():
            cluster_dict[cluster_class] = []
        cluster_dict[cluster_class].append(node)
    return cluster_dict


def get_centor(cluster_dict, k):
    new_centor = []
    for i in range(k):
        centor = np.mean(cluster_dict[i], axis=0)
        new_centor.append(centor)
    return new_centor


def cal_varience(cluster_dict, centor):
    vsum = 0
    for i in range(len(centor)):
        cluster = cluster_dict[i]
        for j in cluster:
            vsum += cal_distance(j, centor[i])
    return vsum


def k_means(data, k):
    centor = random_centor(data, k)
    print(centor)
    cluster_dict = get_cluster(data, centor)
    new_varience = cal_varience(cluster_dict, centor)
    old_varience = 1
    while abs(old_varience - new_varience) > 0.1:
        centor = get_centor(cluster_dict, k)
        cluster_dict = get_cluster(data, centor)
        old_varience = new_varience
        new_varience = cal_varience(cluster_dict, centor)

    #取正
    centor_i = []
    for ind in centor:
        centor_i.append(ind.astype('int64'))

    # 确定每簇都有哪些任务点
    a = []
    # 生成随机颜色
    colors = ['#' + hex(random.randint(0, 0xFFFFFF))[2:].zfill(6) for _ in range(len(globalv.A))]
    for i in range(len(cluster_dict)):
        b = []
        for j in range(len(cluster_dict[i])):
            # plt.scatter(cluster_dict[i][j][0], cluster_dict[i][j][1], c='g', label='target')
            b.append((np.argwhere(globalv.Tar_start == cluster_dict[i][j])[0][0])+1)
            # plt.scatter(cluster_dict[i][j][0],cluster_dict[i][j][1], c = colors[i], label='target')

        a.append(b)
    # plt.show()
    print(centor)
    print(cluster_dict)
    return a, centor_i


# data = np.array([[1, 1, 1], [2, 2, 2], [1, 2, 1], [9, 8, 7], [7, 8, 9], [8, 9, 7]])

# print(globalv.UAV_start)
# print(a, b)