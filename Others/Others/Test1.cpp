#include <iostream>
#include <mplapack_mpfr.h> // MPLAPACK 的 MPFR 支持
#include <mpreal.h>        // MPFR 的 C++ 接口
#include <vector>
#include <cmath>
#include <iomanip>         // 控制输出格式

using namespace std;
using mpreal = mpfr::mpreal;

int main() {
    // 设置高精度（256 位，约 77 位小数）
    mpreal::set_default_prec(256);

    // 参数设置
    const int L = 100;
    mpreal omega = (sqrt(mpreal(5)) - 1) / 2; // 高精度黄金分割比
    mpreal lambda = 1.5;

    // 构造 Hamiltonian 矩阵 H
    mpreal *H = new mpreal[L * L](); // 初始化为 0
    for (int i = 0; i < L; i++) {
        // 对角线元素
        H[i + i * L] = 2 * lambda * cos(2 * M_PI * omega * i);
        // 次对角线元素
        if (i < L - 1) {
            H[i + 1 + i * L] = 1; // 下次对角线
            H[i + (i + 1) * L] = 1; // 上次对角线
        }
    }
    H[L - 1] = pow(10, 4); // 设置 H(1, L) = 10^4

    // 求解特征值和特征向量
    mpreal *E = new mpreal[L]; // 存储特征值
    mpreal *Ev = new mpreal[L * L]; // 存储特征向量（列存储）

    // 使用 MPLAPACK 的对称特征值分解
    char jobz = 'V'; // 计算特征值和特征向量
    char uplo = 'U'; // 使用上三角部分
    mplapackint n = L, lda = L, info;
    mpreal *work = new mpreal[3 * L - 1];

    // 调用 MPLAPACK 特征值分解函数
    Rsyev(jobz, uplo, n, H, lda, E, work, 3 * L - 1, info);

    if (info != 0) {
        cerr << "Error: Eigenvalue computation failed with info = " << info << endl;
        return -1;
    }

    // 输出结果
    cout << fixed << setprecision(50); // 设置输出精度
    vector<int> m_all = {50, 70, 90};

    for (int m : m_all) {
        cout << "Eigenvalue E(" << m << ") = " << E[m - 1] << endl;

        // 输出特征向量的模平方
        cout << "Eigenvector |Ev|^2 for m = " << m << ":" << endl;
        for (int i = 0; i < L; i++) {
            mpreal amplitude = Ev[i + (m - 1) * L]; // 特征向量第 m 列
            cout << i << " " << amplitude * amplitude << endl;
        }
        cout << endl;
    }

    // 释放内存
    delete[] H;
    delete[] E;
    delete[] Ev;
    delete[] work;

    return 0;
}
