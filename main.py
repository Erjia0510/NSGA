from NSGA import NSGA
import matplotlib.pyplot as plt
from improve_NSGA import improve_NSGA

if __name__ == '__main__':
    nodominated = []
    improved = []
    for i in range(30):
        nsga = NSGA()
        best_pop1, first_pop1 = nsga.evolve()
        improve = improve_NSGA()
        best_pop2, first_pop2 = improve.evolve()
        nodominated.extend(best_pop1)
        improved.extend(best_pop2)


    plt.figure()
    plt.scatter([nodominated[i]['fit1'] for i in range(len(nodominated))], [nodominated[i]['fit2'] for i in range(len(nodominated))], color='blue', label='Nondominated Points')
    plt.scatter([improved[i]['fit1'] for i in range(len(improved))],
                [improved[i]['fit2'] for i in range(len(improved))], color='red', label='Nondominated Points')
    plt.title("2D Nondominated Points")
    plt.xlabel("fit1", fontsize=12)
    plt.ylabel("fit2", fontsize=12)
    plt.grid(True)
    plt.savefig(r'./fit.png')
    plt.show()
    # # 初始种群图
    # plt.figure()
    # plt.scatter([first_pop1[i]['fit1'] for i in range(len(first_pop1))],
    #             [first_pop1[i]['fit2'] for i in range(len(first_pop1))], color='blue', label='Nondominated Points')
    # plt.scatter([first_pop2[i]['fit1'] for i in range(len(first_pop2))],
    #             [first_pop2[i]['fit2'] for i in range(len(first_pop2))], color='red', label='Nondominated Points')
    # plt.title("first pop Points")
    # plt.show()

