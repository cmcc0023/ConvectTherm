"""TDMA (Thomas 算法) 三对角矩阵直接求解器。

用于一维问题的线性方程组求解。式(37)。
"""

import numpy as np


def tdma(a_W: np.ndarray, a_P: np.ndarray,
         a_E: np.ndarray, b: np.ndarray) -> np.ndarray:
    """求解 FVM 三对角方程组: -a_W·T_{i-1} + a_P·T_i - a_E·T_{i+1} = b。

    注意：a_W, a_E 为 FVM 对流-扩散系数（均为正值），
    对应矩阵的负次对角线。内部转换为标准 Thomas 算法所需的符号。

    Parameters
    ----------
    a_W : 西侧对流-扩散系数 (长度 N，a_W[0] 无用)
    a_P : 主对角系数 (长度 N)
    a_E : 东侧对流-扩散系数 (长度 N，a_E[-1] 无用)
    b   : 右端项 (长度 N)

    Returns
    -------
    T : 解向量 (长度 N)
    """
    n = len(a_P)
    # 转换为标准三对角形式: c_sub·T_{i-1} + c_main·T_i + c_sup·T_{i+1} = b
    # c_sub = -a_W, c_main = a_P, c_sup = -a_E
    c_sub = -a_W  # 次对角（负值）
    c_sup = -a_E  # 上对角（负值）

    # 前向消元
    c = np.empty(n)
    d = np.empty(n)
    c[0] = c_sup[0] / a_P[0]
    d[0] = b[0] / a_P[0]
    for i in range(1, n):
        denom = a_P[i] - c_sub[i] * c[i - 1]
        c[i] = c_sup[i] / denom
        d[i] = (b[i] - c_sub[i] * d[i - 1]) / denom
    # 回代
    T = np.empty(n)
    T[-1] = d[-1]
    for i in range(n - 2, -1, -1):
        T[i] = d[i] - c[i] * T[i + 1]
    return T
