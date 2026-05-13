"""可视化工具：温度曲线、温度场云图、液相率分布。

所有图形自动保存为 PNG 文件，输出至 output/ 目录。
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

# 尝试加载中文字体
try:
    plt.rcParams["font.sans-serif"] = ["Arial Unicode MS", "SimHei", "DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False
except Exception:
    pass

OUT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "output")


def _ensure_out_dir():
    os.makedirs(OUT_DIR, exist_ok=True)


def _out(filename: str) -> str:
    _ensure_out_dir()
    return os.path.join(OUT_DIR, filename)


def plot_temperature_history(t_hours, T_history, labels, filename="temperature_history.png", title=None):
    """绘制冰块中心温度随时间变化曲线。"""
    fig, ax = plt.subplots(figsize=(8, 5))
    for T_arr, label in zip(T_history, labels):
        ax.plot(t_hours, T_arr, linewidth=1.5, label=label)
    ax.set_xlabel("时间 (h)")
    ax.set_ylabel("温度 (K)")
    ax.set_title(title or "冰块中心温度随时间变化")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(_out(filename), dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  -> output/{filename}")


def plot_temperature_field(T, x, y, title, filename, T_melt=273.15):
    """绘制二维温度场云图，标注相变等温线。"""
    fig, ax = plt.subplots(figsize=(8, 6))
    X, Y = np.meshgrid(x, y)
    vmin = max(T.min(), T_melt - 5)
    vmax = min(T.max(), T_melt + 5)
    norm = TwoSlopeNorm(vmin=vmin, vcenter=T_melt, vmax=vmax)
    pcm = ax.pcolormesh(X * 100, Y * 100, T, cmap="RdBu_r", norm=norm, shading="auto")
    fig.colorbar(pcm, ax=ax, label="温度 (K)")
    cs = ax.contour(X * 100, Y * 100, T, levels=[T_melt], colors="black", linewidths=1.5)
    ax.clabel(cs, fmt="%.1f K", fontsize=8)
    ax.set_xlabel("x (cm)")
    ax.set_ylabel("y (cm)")
    ax.set_title(title)
    ax.set_aspect("equal")
    fig.tight_layout()
    fig.savefig(_out(filename), dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  -> output/{filename}")


def plot_liquid_fraction(f_l, x, y, title, filename):
    """绘制液相率分布图。"""
    fig, ax = plt.subplots(figsize=(8, 6))
    X, Y = np.meshgrid(x, y)
    pcm = ax.pcolormesh(X * 100, Y * 100, f_l, cmap="Blues",
                        vmin=0, vmax=1, shading="auto")
    fig.colorbar(pcm, ax=ax, label="液相率 $f_l$")
    ax.set_xlabel("x (cm)")
    ax.set_ylabel("y (cm)")
    ax.set_title(title)
    ax.set_aspect("equal")
    fig.tight_layout()
    fig.savefig(_out(filename), dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  -> output/{filename}")


def save_data(filename: str, **kwargs):
    """保存温度场历史数据为 .npz 文件。

    Parameters
    ----------
    filename : 输出文件名（不含路径）
    **kwargs : 要保存的数组/标量，如 t, T_history, f_l_history, x, y, bc_type 等
    """
    _ensure_out_dir()
    path = os.path.join(OUT_DIR, filename)
    np.savez_compressed(path, **kwargs)
    print(f"  -> output/{filename}")
