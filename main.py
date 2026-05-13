"""主入口：6 个算例覆盖 3 种边界条件 × 2 种维度。

算例编号与边界条件：
  Case 1  1D-Dirichlet  →  左右定温 293.15 K
  Case 2  1D-Neumann    →  左侧热流 1000 W/m²，右侧定温 293.15 K
  Case 3  1D-Robin      →  左侧对流 h=200, T∞=293.15 K，右侧定温 293.15 K
  Case 4  2D-Dirichlet  →  四边定温 293.15 K
  Case 5  2D-Neumann    →  左右热流 1000 W/m²，上下定温 293.15 K
  Case 6  2D-Robin      →  左右对流 h=200, 上下定温 293.15 K
"""

import os
import sys
import time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from solver.constants import T_MELT
from solver.solver_1d import solve_1d, BC1D
from solver.solver_2d import solve_2d, BC2D, SideBC
from solver.visualize import (
    plot_temperature_history,
    plot_temperature_field,
    plot_liquid_fraction,
    save_data,
    OUT_DIR,
)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def create_ice_mask_2d(x, y, ice_x0, ice_x1, ice_y0, ice_y1):
    """创建二维冰区域布尔掩码。"""
    X, Y = np.meshgrid(x, y)
    return (X >= ice_x0) & (X <= ice_x1) & (Y >= ice_y0) & (Y <= ice_y1)


def find_ice_center_index(x, ice_mask):
    """找到冰区域中心的网格索引。"""
    ice_idx = np.where(ice_mask)[0]
    return ice_idx[len(ice_idx) // 2]


# ────────────────────────────────────────────
# 1D 算例
# ────────────────────────────────────────────

def run_1d_dirichlet():
    """Case 1: 一维 Dirichlet — 左右定温加热。"""
    print("=" * 60)
    print("Case 1: 一维 Dirichlet（左右定温 293.15 K）")
    print("=" * 60)

    L, Nx = 0.1, 100
    dt, Nt_per_call, n_calls = 10.0, 100, 6
    T_h = 293.15

    dx = L / Nx
    x = np.linspace(dx / 2, L - dx / 2, Nx)
    ice_mask = (x >= 0.03) & (x <= 0.07)
    T = np.where(ice_mask, 268.15, 275.15)
    ice_c = find_ice_center_index(x, ice_mask)

    bc = BC1D("dirichlet", T_h, "dirichlet", T_h)

    total_n_picard = 0
    all_converged = True
    t_hist, Tc_hist = [0.0], [T[ice_c]]
    for ci in range(n_calls):
        result = solve_1d(T_init=T, L=L, Nx=Nx,
                          dt=dt, Nt=Nt_per_call, bc=bc,
                          tol_picard=1e-3, max_picard=50)
        total_n_picard += result['n_picard']
        all_converged &= result['converged']
        T = result["T"]
        tc = (ci + 1) * dt * Nt_per_call
        t_hist.append(tc)
        Tc_hist.append(T[ice_c])
        if (ci + 1) % max(1, n_calls // 3) == 0:
            print(f"  t={tc/60:.0f}min  Tc={T[ice_c]:.2f}K  fl={result['f_l'][ice_c]:.3f}")

    result['n_picard'] = total_n_picard
    result['converged'] = all_converged
    t_h = np.array(t_hist) / 3600
    plot_temperature_history(t_h, [Tc_hist], ["定温"], "case1_1d_dirichlet_Tc.png",
                             title="Case 1: 一维定温边界冰块中心温度")
    _plot_1d_profile(x, T, ice_mask, "Case 1: 一维定温边界最终温度分布", "case1_1d_profile.png")
    return result, t_hist, Tc_hist


def run_1d_neumann():
    """Case 2: 一维 Neumann — 左侧给定热流，右侧定温。"""
    print("=" * 60)
    print("Case 2: 一维 Neumann（左 q=1000 W/m²，右 T=293.15 K）")
    print("=" * 60)

    L, Nx = 0.1, 100
    dt, Nt_per_call, n_calls = 10.0, 100, 6
    T_h = 293.15
    q_in = 1000.0

    dx = L / Nx
    x = np.linspace(dx / 2, L - dx / 2, Nx)
    ice_mask = (x >= 0.03) & (x <= 0.07)
    T = np.where(ice_mask, 268.15, 275.15)
    ice_c = find_ice_center_index(x, ice_mask)

    bc = BC1D("neumann", q_in, "dirichlet", T_h)

    total_n_picard = 0
    all_converged = True
    t_hist, Tc_hist = [0.0], [T[ice_c]]
    for ci in range(n_calls):
        result = solve_1d(T_init=T, L=L, Nx=Nx,
                          dt=dt, Nt=Nt_per_call, bc=bc,
                          tol_picard=1e-3, max_picard=50)
        total_n_picard += result['n_picard']
        all_converged &= result['converged']
        T = result["T"]
        tc = (ci + 1) * dt * Nt_per_call
        t_hist.append(tc)
        Tc_hist.append(T[ice_c])
        if (ci + 1) % max(1, n_calls // 3) == 0:
            print(f"  t={tc/60:.0f}min  Tc={T[ice_c]:.2f}K  fl={result['f_l'][ice_c]:.3f}")

    result['n_picard'] = total_n_picard
    result['converged'] = all_converged
    t_h = np.array(t_hist) / 3600
    plot_temperature_history(t_h, [Tc_hist], ["热流"], "case2_1d_neumann_Tc.png",
                             title="Case 2: 一维热流边界冰块中心温度")
    _plot_1d_profile(x, T, ice_mask, "Case 2: 一维热流边界最终温度分布", "case2_1d_profile.png")
    return result, t_hist, Tc_hist


def run_1d_robin():
    """Case 3: 一维 Robin — 左侧对流换热，右侧定温。"""
    print("=" * 60)
    print("Case 3: 一维 Robin（左 h=200 T∞=293.15K，右 T=293.15K）")
    print("=" * 60)

    L, Nx = 0.1, 100
    dt, Nt_per_call, n_calls = 10.0, 100, 6
    T_h = 293.15
    h_conv = 200.0
    T_inf = 293.15

    dx = L / Nx
    x = np.linspace(dx / 2, L - dx / 2, Nx)
    ice_mask = (x >= 0.03) & (x <= 0.07)
    T = np.where(ice_mask, 268.15, 275.15)
    ice_c = find_ice_center_index(x, ice_mask)

    bc = BC1D("robin", (h_conv, T_inf), "dirichlet", T_h)

    total_n_picard = 0
    all_converged = True
    t_hist, Tc_hist = [0.0], [T[ice_c]]
    for ci in range(n_calls):
        result = solve_1d(T_init=T, L=L, Nx=Nx,
                          dt=dt, Nt=Nt_per_call, bc=bc,
                          tol_picard=1e-3, max_picard=50)
        total_n_picard += result['n_picard']
        all_converged &= result['converged']
        T = result["T"]
        tc = (ci + 1) * dt * Nt_per_call
        t_hist.append(tc)
        Tc_hist.append(T[ice_c])
        if (ci + 1) % max(1, n_calls // 3) == 0:
            print(f"  t={tc/60:.0f}min  Tc={T[ice_c]:.2f}K  fl={result['f_l'][ice_c]:.3f}")

    result['n_picard'] = total_n_picard
    result['converged'] = all_converged
    t_h = np.array(t_hist) / 3600
    plot_temperature_history(t_h, [Tc_hist], ["对流"], "case3_1d_robin_Tc.png",
                             title="Case 3: 一维对流边界冰块中心温度")
    _plot_1d_profile(x, T, ice_mask, "Case 3: 一维对流边界最终温度分布", "case3_1d_profile.png")
    return result, t_hist, Tc_hist


# ────────────────────────────────────────────
# 2D 算例
# ────────────────────────────────────────────

def run_2d_dirichlet():
    """Case 4: 二维 Dirichlet — 四边定温。"""
    print("=" * 60)
    print("Case 4: 二维 Dirichlet（四边定温 293.15 K）")
    print("=" * 60)

    Lx, Ly, Nx, Ny = 0.1, 0.1, 50, 50
    dt, Nt_per_call, n_calls = 10.0, 50, 12
    T_h = 293.15

    dx, dy = Lx / Nx, Ly / Ny
    x = np.linspace(dx / 2, Lx - dx / 2, Nx)
    y = np.linspace(dy / 2, Ly - dy / 2, Ny)
    ice_mask = create_ice_mask_2d(x, y, 0.03, 0.07, 0.03, 0.07)
    T = np.where(ice_mask, 268.15, 275.15)
    ice_j = find_ice_center_index(y, ice_mask[:, Nx // 2])
    ice_i = find_ice_center_index(x, ice_mask[Ny // 2, :])

    bc = BC2D(
        left=SideBC("dirichlet", T_wall=T_h),
        right=SideBC("dirichlet", T_wall=T_h),
        bottom=SideBC("dirichlet", T_wall=T_h),
        top=SideBC("dirichlet", T_wall=T_h),
    )

    total_n_picard = 0
    all_converged = True
    all_sor_converged = True
    overall_max_sor_iters = 0
    t_hist, Tc_hist = [0.0], [T[ice_j, ice_i]]
    for ci in range(n_calls):
        result = solve_2d(T_init=T, Lx=Lx, Ly=Ly,
                          Nx=Nx, Ny=Ny, dt=dt, Nt=Nt_per_call, bc=bc,
                          tol_picard=1e-3, max_picard=50,
                          sor_omega=1.5, sor_tol=1e-5, sor_maxiter=2000)
        total_n_picard += result['n_picard']
        all_converged &= result['converged']
        all_sor_converged &= result.get('sor_converged', True)
        overall_max_sor_iters = max(overall_max_sor_iters, result.get('max_sor_iters', 0))
        T = result["T"]
        tc = (ci + 1) * dt * Nt_per_call
        t_hist.append(tc)
        Tc_hist.append(T[ice_j, ice_i])
        if (ci + 1) % max(1, n_calls // 3) == 0:
            print(f"  t={tc/60:.0f}min  Tc={T[ice_j,ice_i]:.2f}K  fl={result['f_l'][ice_j,ice_i]:.3f}")

    result['n_picard'] = total_n_picard
    result['converged'] = all_converged
    result['sor_converged'] = all_sor_converged
    result['max_sor_iters'] = overall_max_sor_iters
    t_h = np.array(t_hist) / 3600
    plot_temperature_history(t_h, [Tc_hist], ["定温"], "case4_2d_dirichlet_Tc.png",
                             title="Case 4: 二维定温边界冰块中心温度")
    plot_temperature_field(T, x, y, "Case 4: 二维定温边界温度场 (K)", "case4_2d_T.png")
    plot_liquid_fraction(result["f_l"], x, y, "Case 4: 二维定温边界液相率分布", "case4_2d_fl.png")
    return result, t_hist, Tc_hist


def run_2d_neumann():
    """Case 5: 二维 Neumann — 左右热流，上下定温。"""
    print("=" * 60)
    print("Case 5: 二维 Neumann（左右 q=1000，上下 T=293.15 K）")
    print("=" * 60)

    Lx, Ly, Nx, Ny = 0.1, 0.1, 50, 50
    dt, Nt_per_call, n_calls = 10.0, 50, 12
    T_h = 293.15
    q_in = 1000.0

    dx, dy = Lx / Nx, Ly / Ny
    x = np.linspace(dx / 2, Lx - dx / 2, Nx)
    y = np.linspace(dy / 2, Ly - dy / 2, Ny)
    ice_mask = create_ice_mask_2d(x, y, 0.03, 0.07, 0.03, 0.07)
    T = np.where(ice_mask, 268.15, 275.15)
    ice_j = find_ice_center_index(y, ice_mask[:, Nx // 2])
    ice_i = find_ice_center_index(x, ice_mask[Ny // 2, :])

    bc = BC2D(
        left=SideBC("neumann", q_in=q_in),
        right=SideBC("neumann", q_in=q_in),
        bottom=SideBC("dirichlet", T_wall=T_h),
        top=SideBC("dirichlet", T_wall=T_h),
    )

    total_n_picard = 0
    all_converged = True
    all_sor_converged = True
    overall_max_sor_iters = 0
    t_hist, Tc_hist = [0.0], [T[ice_j, ice_i]]
    for ci in range(n_calls):
        result = solve_2d(T_init=T, Lx=Lx, Ly=Ly,
                          Nx=Nx, Ny=Ny, dt=dt, Nt=Nt_per_call, bc=bc,
                          tol_picard=1e-3, max_picard=50,
                          sor_omega=1.5, sor_tol=1e-5, sor_maxiter=2000)
        total_n_picard += result['n_picard']
        all_converged &= result['converged']
        all_sor_converged &= result.get('sor_converged', True)
        overall_max_sor_iters = max(overall_max_sor_iters, result.get('max_sor_iters', 0))
        T = result["T"]
        tc = (ci + 1) * dt * Nt_per_call
        t_hist.append(tc)
        Tc_hist.append(T[ice_j, ice_i])
        if (ci + 1) % max(1, n_calls // 3) == 0:
            print(f"  t={tc/60:.0f}min  Tc={T[ice_j,ice_i]:.2f}K  fl={result['f_l'][ice_j,ice_i]:.3f}")

    result['n_picard'] = total_n_picard
    result['converged'] = all_converged
    result['sor_converged'] = all_sor_converged
    result['max_sor_iters'] = overall_max_sor_iters
    t_h = np.array(t_hist) / 3600
    plot_temperature_history(t_h, [Tc_hist], ["热流"], "case5_2d_neumann_Tc.png",
                             title="Case 5: 二维热流边界冰块中心温度")
    plot_temperature_field(T, x, y, "Case 5: 二维热流边界温度场 (K)", "case5_2d_T.png")
    plot_liquid_fraction(result["f_l"], x, y, "Case 5: 二维热流边界液相率分布", "case5_2d_fl.png")
    return result, t_hist, Tc_hist


def run_2d_robin():
    """Case 6: 二维 Robin — 左右对流换热，上下定温。"""
    print("=" * 60)
    print("Case 6: 二维 Robin（左右 h=200 T∞=293.15K，上下 T=293.15K）")
    print("=" * 60)

    Lx, Ly, Nx, Ny = 0.1, 0.1, 50, 50
    dt, Nt_per_call, n_calls = 10.0, 50, 12
    T_h = 293.15
    h_conv = 200.0
    T_inf = 293.15

    dx, dy = Lx / Nx, Ly / Ny
    x = np.linspace(dx / 2, Lx - dx / 2, Nx)
    y = np.linspace(dy / 2, Ly - dy / 2, Ny)
    ice_mask = create_ice_mask_2d(x, y, 0.03, 0.07, 0.03, 0.07)
    T = np.where(ice_mask, 268.15, 275.15)
    ice_j = find_ice_center_index(y, ice_mask[:, Nx // 2])
    ice_i = find_ice_center_index(x, ice_mask[Ny // 2, :])

    bc = BC2D(
        left=SideBC("robin", h=h_conv, T_inf=T_inf),
        right=SideBC("robin", h=h_conv, T_inf=T_inf),
        bottom=SideBC("dirichlet", T_wall=T_h),
        top=SideBC("dirichlet", T_wall=T_h),
    )

    total_n_picard = 0
    all_converged = True
    all_sor_converged = True
    overall_max_sor_iters = 0
    t_hist, Tc_hist = [0.0], [T[ice_j, ice_i]]
    for ci in range(n_calls):
        result = solve_2d(T_init=T, Lx=Lx, Ly=Ly,
                          Nx=Nx, Ny=Ny, dt=dt, Nt=Nt_per_call, bc=bc,
                          tol_picard=1e-3, max_picard=50,
                          sor_omega=1.5, sor_tol=1e-5, sor_maxiter=2000)
        total_n_picard += result['n_picard']
        all_converged &= result['converged']
        all_sor_converged &= result.get('sor_converged', True)
        overall_max_sor_iters = max(overall_max_sor_iters, result.get('max_sor_iters', 0))
        T = result["T"]
        tc = (ci + 1) * dt * Nt_per_call
        t_hist.append(tc)
        Tc_hist.append(T[ice_j, ice_i])
        if (ci + 1) % max(1, n_calls // 3) == 0:
            print(f"  t={tc/60:.0f}min  Tc={T[ice_j,ice_i]:.2f}K  fl={result['f_l'][ice_j,ice_i]:.3f}")

    result['n_picard'] = total_n_picard
    result['converged'] = all_converged
    result['sor_converged'] = all_sor_converged
    result['max_sor_iters'] = overall_max_sor_iters
    t_h = np.array(t_hist) / 3600
    plot_temperature_history(t_h, [Tc_hist], ["对流"], "case6_2d_robin_Tc.png",
                             title="Case 6: 二维对流边界冰块中心温度")
    plot_temperature_field(T, x, y, "Case 6: 二维对流边界温度场 (K)", "case6_2d_T.png")
    plot_liquid_fraction(result["f_l"], x, y, "Case 6: 二维对流边界液相率分布", "case6_2d_fl.png")
    return result, t_hist, Tc_hist


# ────────────────────────────────────────────
# 对比图
# ────────────────────────────────────────────

def plot_comparison(results_1d, results_2d):
    """绘制 6 个算例的冰块中心温度对比曲线。"""
    labels = ["定温", "热流", "对流"]

    # 1D 对比
    fig, ax = plt.subplots(figsize=(8, 5))
    for res, label in zip(results_1d, labels):
        ax.plot(np.array(res["t_hist"]) / 3600, res["Tc_hist"],
                linewidth=1.5, label=label)
    ax.set_xlabel("时间 (h)")
    ax.set_ylabel("冰块中心温度 (K)")
    ax.set_title("一维算例：三种边界条件冰块中心温度对比")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "comparison_1d.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("  -> output/comparison_1d.png")

    # 2D 对比
    fig, ax = plt.subplots(figsize=(8, 5))
    for res, label in zip(results_2d, labels):
        ax.plot(np.array(res["t_hist"]) / 3600, res["Tc_hist"],
                linewidth=1.5, label=label)
    ax.set_xlabel("时间 (h)")
    ax.set_ylabel("冰块中心温度 (K)")
    ax.set_title("二维算例：三种边界条件冰块中心温度对比")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "comparison_2d.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("  -> output/comparison_2d.png")


def _plot_1d_profile(x, T, ice_mask, title, filename):
    """一维最终温度分布折线图。"""
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(x * 100, T, "b-", linewidth=1.5)
    ax.axhline(y=T_MELT, color="gray", linestyle="--", alpha=0.5, label=f"熔点 T$_m$={T_MELT} K")
    ice_x = x[ice_mask] * 100
    if len(ice_x) > 0:
        ax.axvspan(ice_x[0], ice_x[-1], alpha=0.1, color="blue", label="冰区域")
    ax.set_xlabel("x (cm)")
    ax.set_ylabel("温度 T (K)")
    ax.set_title(title)
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, filename), dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  -> output/{filename}")


# ────────────────────────────────────────────
# Main
# ────────────────────────────────────────────

def main():
    print("对流传热相变数值求解器")
    print("物理场景：矩形域水-冰系统，6 个算例覆盖 3 种边界条件 × 2 种维度\n")

    t_start = time.time()

    # 1D 算例
    r1, t1, Tc1 = run_1d_dirichlet()
    r2, t2, Tc2 = run_1d_neumann()
    r3, t3, Tc3 = run_1d_robin()

    # 2D 算例
    r4, t4, Tc4 = run_2d_dirichlet()
    r5, t5, Tc5 = run_2d_neumann()
    r6, t6, Tc6 = run_2d_robin()

    # 对比图
    results_1d = [{"t_hist": t1, "Tc_hist": Tc1}, {"t_hist": t2, "Tc_hist": Tc2}, {"t_hist": t3, "Tc_hist": Tc3}]
    results_2d = [{"t_hist": t4, "Tc_hist": Tc4}, {"t_hist": t5, "Tc_hist": Tc5}, {"t_hist": t6, "Tc_hist": Tc6}]
    plot_comparison(results_1d, results_2d)

    # 保存温度场历史数据
    all_results = [
        ("case1_1d_dirichlet", r1, t1, Tc1),
        ("case2_1d_neumann",   r2, t2, Tc2),
        ("case3_1d_robin",     r3, t3, Tc3),
        ("case4_2d_dirichlet", r4, t4, Tc4),
        ("case5_2d_neumann",   r5, t5, Tc5),
        ("case6_2d_robin",     r6, t6, Tc6),
    ]
    for name, r, th, tch in all_results:
        save_data(
            f"{name}.npz",
            t=np.array(th),
            Tc=np.array(tch),
            T_final=r["T"],
            f_l_final=r["f_l"],
            x=r["x"],
            y=r.get("y", np.array([])),
            dx=r["dx"],
            dy=r.get("dy", np.array([])),
        )

    elapsed = time.time() - t_start
    print(f"\n全部算例完成，总耗时 {elapsed:.1f}s")

    # 打印汇总表
    cases = [
        ("1D-Dirichlet", r1), ("1D-Neumann", r2), ("1D-Robin", r3),
        ("2D-Dirichlet", r4), ("2D-Neumann", r5), ("2D-Robin", r6),
    ]
    print(f"\n{'Case':<18} {'T_min(K)':>8} {'T_max(K)':>8} {'Tc(K)':>8} {'f_l_c':>7} {'Picard':>7} {'Conv':>5} {'SOR':>5}")
    print("-" * 73)
    for name, r in cases:
        T = r["T"]
        fl = r["f_l"]
        # 中心温度和中心液相率
        if T.ndim == 1:
            Tc = T[len(T) // 2]
            flc = fl[len(fl) // 2]
        else:
            Tc = T[T.shape[0] // 2, T.shape[1] // 2]
            flc = fl[fl.shape[0] // 2, fl.shape[1] // 2]
        sor_info = f"{r['sor_converged']!s:>5}" if 'sor_converged' in r else "    -"
        print(f"{name:<18} {T.min():8.2f} {T.max():8.2f} {Tc:8.2f} {flc:7.3f} {r['n_picard']:7d} {r['converged']!s:>5} {sor_info}")


if __name__ == "__main__":
    main()
