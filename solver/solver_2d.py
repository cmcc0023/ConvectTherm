"""二维隐式 FVM 求解器。

离散格式：全隐式 + 一阶迎风对流 + 中心差分扩散 + 调和平均界面导热。
线性求解：逐线 TDMA + SOR (SIP)。
非线性处理：Picard 迭代。

对应公式：式(20)–(30)，边界条件式(31)–(36)。
"""

from dataclasses import dataclass, field
import numpy as np

from .constants import RHO_W, K_W, CP_W, CP_I
from .physics import liquid_fraction, apparent_cp, interp_prop, effective_conductivity
from .tdma import tdma


@dataclass
class SideBC:
    """单侧边界条件。

    type: "dirichlet" | "neumann" | "robin"
    T_wall: 定温边界温度 (K)，dirichlet 时使用
    q_in: 热流密度 (W/m²)，neumann 时使用，进入计算域为正
    h: 对流换热系数 (W/(m²·K))，robin 时使用
    T_inf: 环境温度 (K)，robin 时使用
    """
    type: str = "dirichlet"
    T_wall: float = 293.15
    q_in: float = 0.0
    h: float = 0.0
    T_inf: float = 293.15


@dataclass
class BC2D:
    """二维四边边界条件。"""
    left: SideBC = field(default_factory=SideBC)
    right: SideBC = field(default_factory=SideBC)
    bottom: SideBC = field(default_factory=SideBC)
    top: SideBC = field(default_factory=SideBC)


def solve_2d(
    T_init: np.ndarray,
    Lx: float,
    Ly: float,
    Nx: int,
    Ny: int,
    dt: float,
    Nt: int,
    bc: BC2D | None = None,
    T_boundary: float | None = None,
    u_water: float = 0.0,
    v_water: float = 0.0,
    S_val: float = 0.0,
    tol_picard: float = 1e-4,
    max_picard: int = 100,
    alpha_picard: float = 0.5,
    sor_omega: float = 1.5,
    sor_tol: float = 1e-6,
    sor_maxiter: int = 5000,
) -> dict:
    """二维隐式 FVM 求解器。

    Parameters
    ----------
    T_init      : 初始温度场 (K)，形状 (Ny, Nx)
    Lx, Ly      : 计算域尺寸 (m)
    Nx, Ny      : x, y 方向网格数
    dt          : 时间步长 (s)
    Nt          : 总时间步数
    bc          : 四边边界条件 (BC2D)，优先使用
    T_boundary  : 四边定温温度 (K)，当 bc=None 时使用（兼容旧接口）
    u_water     : 水中 x 方向参考速度 (m/s)
    v_water     : 水中 y 方向参考速度 (m/s)
    S_val       : 体积热源 (W/m³)
    tol_picard  : Picard 迭代收敛容差 (K)
    max_picard  : 每时间步最大 Picard 迭代次数
    alpha_picard: Picard 欠松弛因子 (0,1]，1 为无松弛
    sor_omega   : SOR 松弛因子
    sor_tol     : SOR 收敛容差
    sor_maxiter : SOR 最大迭代次数

    Returns
    -------
    dict
    """
    # 兼容旧接口：若传入 T_boundary 而非 bc，构建四边 Dirichlet 的 BC2D
    if bc is None:
        T_val = T_boundary if T_boundary is not None else 293.15
        bc = BC2D(
            left=SideBC("dirichlet", T_wall=T_val),
            right=SideBC("dirichlet", T_wall=T_val),
            bottom=SideBC("dirichlet", T_wall=T_val),
            top=SideBC("dirichlet", T_wall=T_val),
        )

    dx = Lx / Nx
    dy = Ly / Ny
    x = np.linspace(dx / 2, Lx - dx / 2, Nx)
    y = np.linspace(dy / 2, Ly - dy / 2, Ny)

    T = T_init.copy()
    f_l = liquid_fraction(T)
    rho = interp_prop(f_l, 917.0, RHO_W)
    k = interp_prop(f_l, 2.25, K_W)
    C_app = rho * apparent_cp(f_l)
    U = f_l * u_water
    V = f_l * v_water

    S_p_line = 0.0
    S_p_const = S_val

    total_picard = 0
    all_converged = True
    sor_all_converged = True
    max_sor_iters_used = 0

    for n in range(Nt):
        T_old = T.copy()

        # Picard 迭代
        for m in range(max_picard):
            T_prev = T.copy()
            f_l_prev = f_l.copy()
            rho_prev = rho.copy()
            k_prev = k.copy()
            C_app_prev = C_app.copy()

            # --- 内部界面参数 ---
            cp = interp_prop(f_l, CP_I, CP_W)
            C = rho * cp

            k_xint = effective_conductivity(k[:, :-1], k[:, 1:])
            C_xint = 0.5 * (C[:, :-1] + C[:, 1:])
            u_xint = 0.5 * (U[:, :-1] + U[:, 1:])
            F_xint = C_xint * u_xint * dy
            D_xint = k_xint * dy / dx

            k_yint = effective_conductivity(k[:-1, :], k[1:, :])
            C_yint = 0.5 * (C[:-1, :] + C[1:, :])
            v_yint = 0.5 * (V[:-1, :] + V[1:, :])
            F_yint = C_yint * v_yint * dx
            D_yint = k_yint * dx / dy

            D_east = np.zeros((Ny, Nx))
            D_west = np.zeros((Ny, Nx))
            F_east = np.zeros((Ny, Nx))
            F_west = np.zeros((Ny, Nx))
            D_east[:, :-1] = D_xint
            D_west[:, 1:] = D_xint
            F_east[:, :-1] = F_xint
            F_west[:, 1:] = F_xint

            D_north = np.zeros((Ny, Nx))
            D_south = np.zeros((Ny, Nx))
            F_north = np.zeros((Ny, Nx))
            F_south = np.zeros((Ny, Nx))
            D_north[:-1, :] = D_yint
            D_south[1:, :] = D_yint
            F_north[:-1, :] = F_yint
            F_south[1:, :] = F_yint

            # --- 逐线 TDMA + SOR ---
            T_new = T.copy()
            for sor_iter in range(sor_maxiter):
                T_sor_prev = T_new.copy()

                for j in range(Ny):
                    a_W = D_west[j, :] + np.maximum(F_west[j, :], 0.0)
                    a_E = D_east[j, :] + np.maximum(-F_east[j, :], 0.0)
                    a_S = D_south[j, :] + np.maximum(F_south[j, :], 0.0)
                    a_N = D_north[j, :] + np.maximum(-F_north[j, :], 0.0)

                    a_P0 = C_app[j, :] * dx * dy / dt
                    rhs = a_P0 * T_old[j, :] + S_p_const * dx * dy

                    if j > 0:
                        rhs += a_S * T_new[j - 1, :]
                    if j < Ny - 1:
                        rhs += a_N * T_new[j + 1, :]

                    a_P = a_P0 + a_W + a_E + a_S + a_N + (F_east[j, :] - F_west[j, :]) + (F_north[j, :] - F_south[j, :]) - S_p_line * dx * dy

                    # === 左边界 (i=0) ===
                    if bc.left.type == "dirichlet":
                        D_bl = 2.0 * k[j, 0] * dy / dx
                        a_P[0] += D_bl
                        rhs[0] += D_bl * bc.left.T_wall
                    elif bc.left.type == "neumann":
                        rhs[0] += bc.left.q_in * dy
                    elif bc.left.type == "robin":
                        h_eff = 1.0 / (dx / (2.0 * k[j, 0]) + 1.0 / bc.left.h)
                        a_P[0] += h_eff * dy
                        rhs[0] += h_eff * dy * bc.left.T_inf

                    # === 右边界 (i=Nx-1) ===
                    if bc.right.type == "dirichlet":
                        D_br = 2.0 * k[j, -1] * dy / dx
                        a_P[-1] += D_br
                        rhs[-1] += D_br * bc.right.T_wall
                    elif bc.right.type == "neumann":
                        rhs[-1] += bc.right.q_in * dy
                    elif bc.right.type == "robin":
                        h_eff = 1.0 / (dx / (2.0 * k[j, -1]) + 1.0 / bc.right.h)
                        a_P[-1] += h_eff * dy
                        rhs[-1] += h_eff * dy * bc.right.T_inf

                    # === 下边界 (j=0) ===
                    if j == 0:
                        if bc.bottom.type == "dirichlet":
                            D_bs = 2.0 * k[0, :] * dx / dy
                            a_P += D_bs
                            rhs += D_bs * bc.bottom.T_wall
                        elif bc.bottom.type == "neumann":
                            rhs += bc.bottom.q_in * dx
                        elif bc.bottom.type == "robin":
                            h_eff = 1.0 / (dy / (2.0 * k[0, :]) + 1.0 / bc.bottom.h)
                            a_P += h_eff * dx
                            rhs += h_eff * dx * bc.bottom.T_inf

                    # === 上边界 (j=Ny-1) ===
                    if j == Ny - 1:
                        if bc.top.type == "dirichlet":
                            D_bn = 2.0 * k[-1, :] * dx / dy
                            a_P += D_bn
                            rhs += D_bn * bc.top.T_wall
                        elif bc.top.type == "neumann":
                            rhs += bc.top.q_in * dx
                        elif bc.top.type == "robin":
                            h_eff = 1.0 / (dy / (2.0 * k[-1, :]) + 1.0 / bc.top.h)
                            a_P += h_eff * dx
                            rhs += h_eff * dx * bc.top.T_inf

                    T_line = tdma(a_W, a_P, a_E, rhs)
                    T_new[j, :] = (1.0 - sor_omega) * T_new[j, :] + sor_omega * T_line

                max_change = np.max(np.abs(T_new - T_sor_prev))
                if max_change < sor_tol:
                    break

            # SOR 收敛追踪
            sor_iters_this = sor_iter + 1 if max_change < sor_tol else sor_maxiter
            max_sor_iters_used = max(max_sor_iters_used, sor_iters_this)
            if sor_iters_this == sor_maxiter and max_change >= sor_tol:
                sor_all_converged = False

            # 欠松弛
            T = (1.0 - alpha_picard) * T_prev + alpha_picard * T_new

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
            V = f_l * v_water

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
        V = f_l * v_water

    return {
        "T": T, "f_l": f_l, "C_app": C_app,
        "k": k, "rho": rho, "U": U, "V": V,
        "x": x, "y": y, "dx": dx, "dy": dy,
        "n_picard": total_picard,
        "converged": all_converged and sor_all_converged,
        "sor_converged": sor_all_converged,
        "max_sor_iters": max_sor_iters_used,
    }
