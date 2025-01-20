# We use the package mpmath to calculate the distribution of the eigenstates
# We can change the boundary condition 


import mpmath as mp
import csv
import pandas
import numpy as np
import matplotlib.pyplot as plt

mp.dps=100

pi=mp.pi
omega=(mp.sqrt(5)-1)/2

def Eigensate_SFL():
 # we plot the typical eigenstate in the SFL region 
    L=100
    gamma=10000
    gamma=10^(-4)
    Lambda=0.5

    H=mp.zeros(L,L)
    for m in range(L-1):
        H[m+1,m]=1
        H[m,m]=2*Lambda*mp.cos(2*pi*omega*(m+1))
    H[L-1,L-1]=2*Lambda*mp.cos(2*pi*omega*(L))
    H[0,L-1]=gamma  # boudnary  condition 
    
    E,ER=mp.eig(H)

    # 将本征向量转为 NumPy 数组
    ER = np.array(ER.tolist(), dtype=complex)

    # 选择第 50、70 和 90 个本征态（注意索引从 0 开始）
    indices = [49, 69, 89]  # 第 50、70 和 90 个
    if max(indices) >= len(ER):
        print(f"矩阵规模为 {len(ER)}，无法访问超过范围的本征态！")
        return

    # 绘制本征态分布
    plt.figure(figsize=(10, 6))
    for idx in indices:
        plt.plot(range(L), np.abs(ER[:, idx])**2, label=f"Eigenstate {idx + 1}")

    # y轴坐标按照对数来画
    plt.yscale('log')
    # 设置图例和标签
    plt.xlabel("Site index")
    plt.ylabel("Probability density")
    plt.title("Eigenstate distributions")
    plt.legend()
    plt.grid()
    plt.show()
    
    print(H)


    return 


def Eigensate_AL():

    # 参数设置
    L = 100
    gamma = mp.mpf("10")**(-4)
    Lambda = mp.mpf("1.5")
    omega = (mp.sqrt(5)-1) / 2

    # 构建哈密顿量矩阵 H
    H = mp.zeros(L, L)
    for m in range(L - 1):
        H[m + 1, m] = 1
        H[m, m] = 2 * Lambda * mp.cos(2 * pi * omega * (m + 1))
    H[L - 1, L - 1] = 2 * Lambda * mp.cos(2 * pi * omega * L)
    H[0, L - 1] = gamma  # 边界条件

    # 计算本征值和本征向量
    E, ER = mp.eig(H)

    # 选择第 50、70 和 90 个本征态
    indices1 = [50, 70, 90]
    indices = []
    for x in indices1:
        target = 2 * Lambda * mp.cos(2 * omega * pi * (L - x))
        closest_index = min(range(len(E)), key=lambda i: abs(E[i] - target))
        indices.append(closest_index)
        print(f"Target energy: {target}, Closest eigenvalue: {E[closest_index]}")

    # 绘制本征态分布
    plt.figure(figsize=(10, 6))
    for idx in indices:
        eigenstate = [abs(ER[m, idx])**2 for m in range(L)]
        plt.plot(range(L), eigenstate, label=f"Eigenstate {idx + 1}")

    # 设置对数 y 轴
    plt.yscale('log')
    plt.xlabel("Site index")
    plt.ylabel("Probability density")
    plt.title("Eigenstate distributions")
    plt.legend()
    plt.grid()
    plt.show()



    return 


def find_closest_index(sequence, target):
    """
    在序列中找到最接近目标值的索引。

    参数:
        sequence (list or array): 数字序列。
        target (float): 目标值。

    返回:
        int: 最接近目标值的索引。
    """
    closest_index = min(range(len(sequence)), key=lambda i: abs(sequence[i] - target))
    return closest_index





if __name__=="__main__":
    Eigensate_AL()