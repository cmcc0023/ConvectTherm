"""一维隐式 FVM 求解器。

离散格式：全隐式 + 一阶迎风对流 + 中心差分扩散 + 调和平均界面导热。
线性求解：TDMA（追赶法）。
非线性处理：Picard 迭代。

对应公式：式(8)–(19)，边界条件式(31)–(36)。
"""

from dataclasses import dataclass
import numpy as np

from .constants import RHO_W, K_W, CP_W, CP_I, T_MELT
from .physics import liquid_fraction, apparent_cp, interp_prop, effective_conductivity
from .tdma import tdma


@dataclass
class BC1D:
    """一维边界条件。

    dirichlet: value = T_wall
    neumann:   value = q_in (进入计算域为正)
    robin:     value = (h, T_inf) 元组
    """
    left_type: str = "dirichlet"
    left_value: float | tuple = 0.0
    right_type: str = "dirichlet"
    right_value: float | tuple = 0.0


def solve_1d(
    T_init: np.ndarray,
    L: float,
    Nx: int,
    dt: float,
    Nt: int,
    bc: BC1D,
    u_water: float = 0.0,
    S_val: float = 0.0,
    tol_picard: float = 1e-4,
    max_picard: int = 100,
    alpha_picard: float = 0.5,
) -> dict:
    """一维隐式 FVM 求解器。

    Parameters
    ----------
    T_init      : 初始温度场 (K)，长度 Nx
    L           : 计算域长度 (m)
    Nx          : 网格数
    dt          : 时间步长 (s)
    Nt          : 总时间步数
    bc          : 边界条件
    u_water     : 水中参考速度 (m/s)
    S_val       : 体积热源 (W/m³)，与温度无关
    tol_picard  : Picard 迭代收敛容差 (K)
    max_picard  : 每时间步最大 Picard 迭代次数
    alpha_picard: Picard 欠松弛因子 (0,1]，同时作用于温度和物性

    Returns
    -------
    dict : 含 T, f_l, C_app, k, rho, U, n_picard, converged
    """
    dx = L / Nx
    x = np.linspace(dx / 2, L - dx / 2, Nx)

    # 初始化物性
    T = T_init.copy()
    f_l = liquid_fraction(T)
    rho = interp_prop(f_l, 917.0, RHO_W)
    k = interp_prop(f_l, 2.25, K_W)
    C_app = rho * apparent_cp(f_l)
    U = f_l * u_water

    S_p_line = 0.0   # 源项线性化系数 S_P'，须 ≤ 0
    S_p_const = S_val # 源项常数部分 S_P''

    total_picard = 0
    all_converged = True

    for n in range(Nt):
        T_old = T.copy()

        # Picard 迭代
        for m in range(max_picard):
            T_prev = T.copy()
            f_l_prev = f_l.copy()
            rho_prev = rho.copy()
            k_prev = k.copy()
            C_app_prev = C_app.copy()

            # 计算内部界面参数（cell i 与 cell i+1 之间，共 Nx-1 个）
            k_int = effective_conductivity(k[:-1], k[1:])  # 调和平均，式(26)
            cp = interp_prop(f_l, CP_I, CP_W)             # 比热容
            C = rho * cp                                    # 体积热容 C = ρ·cp
            C_int = 0.5 * (C[:-1] + C[1:])
            u_int = 0.5 * (U[:-1] + U[1:])
            F_int = C_int * u_int    # 对流通量系数，式(10)，F = C·u
            D_int = k_int / dx       # 扩散导热系数，式(11)

            # 组装三对角系数
            D_east = np.zeros(Nx)
            D_west = np.zeros(Nx)
            F_east = np.zeros(Nx)
            F_west = np.zeros(Nx)

            D_east[:-1] = D_int
            D_west[1:] = D_int
            F_east[:-1] = F_int
            F_west[1:] = F_int

            a_W = D_west + np.maximum(F_west, 0.0)   # 式(16a)
            a_E = D_east + np.maximum(-F_east, 0.0)   # 式(16b)

            # 时间项，式(16c)
            a_P0 = C_app * dx / dt
            b = a_P0 * T_old + S_p_const * dx   # 式(16e)

            # 主对角系数（不可压流简化），式(17)
            a_P = a_P0 + a_W + a_E + (F_east - F_west) - S_p_line * dx

            # --- 边界条件修正 ---
            if bc.left_type == "dirichlet":
                D_b = 2.0 * k[0] / dx
                a_P[0] += D_b
                b[0] += D_b * bc.left_value
            elif bc.left_type == "neumann":
                b[0] += bc.left_value
            elif bc.left_type == "robin":
                h, T_inf = bc.left_value
                h_eff = 1.0 / (dx / (2.0 * k[0]) + 1.0 / h)
                a_P[0] += h_eff
                b[0] += h_eff * T_inf

            if bc.right_type == "dirichlet":
                D_b = 2.0 * k[-1] / dx
                a_P[-1] += D_b
                b[-1] += D_b * bc.right_value
            elif bc.right_type == "neumann":
                b[-1] += bc.right_value
            elif bc.right_type == "robin":
                h, T_inf = bc.right_value
                h_eff = 1.0 / (dx / (2.0 * k[-1]) + 1.0 / h)
                a_P[-1] += h_eff
                b[-1] += h_eff * T_inf

            # 求解 + 欠松弛
            T_solved = tdma(a_W, a_P, a_E, b)
            T = (1.0 - alpha_picard) * T_prev + alpha_picard * T_solved

            # 更新物性并对物性做欠松弛（抑制相变区震荡）
            f_l_new = liquid_fraction(T)
            rho_new = interp_prop(f_l_new, 917.0, RHO_W)
            k_new = interp_prop(f_l_new, 2.25, K_W)
            C_app_new = rho_new * apparent_cp(f_l_new)
            f_l = (1.0 - alpha_picard) * f_l_prev + alpha_picard * f_l_new
            rho = (1.0 - alpha_picard) * rho_prev + alpha_picard * rho_new
            k = (1.0 - alpha_picard) * k_prev + alpha_picard * k_new
            C_app = (1.0 - alpha_picard) * C_app_prev + alpha_picard * C_app_new
            U = f_l * u_water

            # 收敛检查：比较相邻 Picard 迭代
            if np.max(np.abs(T - T_prev)) < tol_picard:
                total_picard += m + 1
                break
        else:
            total_picard += max_picard
            all_converged = False

        # 时间步结束，更新物性（确保返回值一致）
        f_l = liquid_fraction(T)
        rho = interp_prop(f_l, 917.0, RHO_W)
        k = interp_prop(f_l, 2.25, K_W)
        C_app = rho * apparent_cp(f_l)
        U = f_l * u_water

    return {
        "T": T, "f_l": f_l, "C_app": C_app,
        "k": k, "rho": rho, "U": U,
        "x": x, "dx": dx,
        "n_picard": total_picard, "converged": all_converged,
    }
