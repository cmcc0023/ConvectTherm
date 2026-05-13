"""相变物理模型：液相率、表观热容、物性插值。

对应公式编号参见《控制方程的统一形式》：
  - 液相率 f_l: 式(4)
  - 表观比热容 c_p_app: 式(5)
  - 物性插值 k, rho: 式(6)
  - 速度场 U = f_l U_w: 式(7)
"""

import numpy as np

from .constants import (
    RHO_W, RHO_I, K_W, K_I, CP_W, CP_I,
    T_MELT, DELTA_T_M, L_FUSION, EPS_T,
)


def liquid_fraction(T: np.ndarray) -> np.ndarray:
    """计算液相率 f_l，线性插值于 [T_s, T_l] 区间。式(4)。"""
    T_s = T_MELT - DELTA_T_M / 2.0
    T_l = T_MELT + DELTA_T_M / 2.0
    f_l = np.clip((T - T_s) / (T_l - T_s), 0.0, 1.0)
    return f_l


def apparent_cp(f_l: np.ndarray) -> np.ndarray:
    """表观比热容 c_p_app，相变区间内含潜热项。式(5)。"""
    T_s = T_MELT - DELTA_T_M / 2.0
    T_l = T_MELT + DELTA_T_M / 2.0
    solid = f_l < EPS_T
    liquid = f_l > 1.0 - EPS_T
    mushy = ~solid & ~liquid

    cp_app = np.where(solid, CP_I,
             np.where(liquid, CP_W, 0.0))
    cp_app[mushy] = (
        (1.0 - f_l[mushy]) * CP_I
        + f_l[mushy] * CP_W
        + L_FUSION / (T_l - T_s)
    )
    return cp_app


def interp_prop(f_l: np.ndarray, val_solid: float, val_liquid: float) -> np.ndarray:
    """按液相率线性插值物性参数。式(6)。"""
    return (1.0 - f_l) * val_solid + f_l * val_liquid


def effective_conductivity(k_P: np.ndarray, k_N: np.ndarray) -> np.ndarray:
    """界面调和平均导热系数。式(26)。"""
    return 2.0 * k_P * k_N / (k_P + k_N + 1e-30)
