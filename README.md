# ConvectTherm — 对流传热相变数值模拟

## 项目概述

本项目针对 **含相变的非稳态对流-扩散传热问题**，采用有限体积法（FVM）编制数值求解程序。

<img src="images/题目方程与条件.png" alt="题目方程与条件" style="zoom:33%;" />

**物理场景**：矩形域内装水，水中含矩形冰块，通过不同边界条件加热，求解瞬态温度场分布，分析冰块融化过程。



<img src="images/算例.png" alt="算例" style="zoom:33%;" />

## 公式推导

### 1. 控制方程

含相变的对流-扩散传热统一形式：

<img src="https://latex.codecogs.com/png.image?C_%5Cmathrm%7Bapp%7D%5Cfrac%7B%5Cpartial%20T%7D%7B%5Cpartial%20t%7D%20%2B%20%5Cnabla%5Ccdot%28C%5Cmathbf%7BU%7D%20T%29%20%3D%20%5Cnabla%5Ccdot%28k%5Cnabla%20T%29%20%2B%20S" alt="formula" />

其中 <img src="https://latex.codecogs.com/png.image?C%20%3D%20%5Crho%20c_p" alt="formula" /> 为体积热容，<img src="https://latex.codecogs.com/png.image?C_%5Cmathrm%7Bapp%7D" alt="formula" /> 为表观热容（含相变潜热），<img src="https://latex.codecogs.com/png.image?%5Cmathbf%7BU%7D%20%3D%20f_l%20%5Cmathbf%7BU%7D_w" alt="formula" /> 为液相速度场。

### 2. 表观热容法处理相变

在熔点 <img src="https://latex.codecogs.com/png.image?T_m" alt="formula" /> 两侧引入相变区间 <img src="https://latex.codecogs.com/png.image?%5BT_s%2C%20T_l%5D" alt="formula" />，避免显式追踪冰-水界面：

<img src="https://latex.codecogs.com/png.image?T_s%20%3D%20T_m%20-%20%5Cfrac%7B%5CDelta%20T_m%7D%7B2%7D%2C%20%5Cquad%20T_l%20%3D%20T_m%20%2B%20%5Cfrac%7B%5CDelta%20T_m%7D%7B2%7D" alt="formula" />

**液相率** 在区间内线性插值：

<img src="https://latex.codecogs.com/png.image?f_l%20%3D%20%5Cbegin%7Bcases%7D%200%2C%20%26%20T%20%5Cle%20T_s%20%5C%5C%20%5Cfrac%7BT%20-%20T_s%7D%7BT_l%20-%20T_s%7D%2C%20%26%20T_s%20%3C%20T%20%3C%20T_l%20%5C%5C%201%2C%20%26%20T%20%5Cge%20T_l%20%5Cend%7Bcases%7D" alt="formula" />

**表观比热容** 在相变区间内叠加潜热项：

<img src="https://latex.codecogs.com/png.image?c_%7Bp%2C%5Cmathrm%7Bapp%7D%7D%20%3D%20c_p%20%2B%20%5Cfrac%7BL%7D%7BT_l%20-%20T_s%7D" alt="formula" />

**物性插值**：<img src="https://latex.codecogs.com/png.image?k%20%3D%20f_l%20k_w%20%2B%20%281-f_l%29%20k_i" alt="formula" />，<img src="https://latex.codecogs.com/png.image?%5Crho%20%3D%20f_l%20%5Crho_w%20%2B%20%281-f_l%29%20%5Crho_i" alt="formula" />

### 3. 有限体积法离散

对控制方程在控制体 <img src="https://latex.codecogs.com/png.image?V_P" alt="formula" /> 上积分，时间项全隐式，对流项一阶迎风，扩散项中心差分：

<img src="https://latex.codecogs.com/png.image?C_%7B%5Cmathrm%7Bapp%7D%2CP%7D%20V_P%20%5Cfrac%7BT_P%5E%7Bn%2B1%7D-T_P%5En%7D%7B%5CDelta%20t%7D%20%2B%20%28F_e%20T_e%20-%20F_w%20T_w%29%20%3D%20D_e%28T_E%20-%20T_P%29%20-%20D_w%28T_P%20-%20T_W%29%20%2B%20S_P%20V_P" alt="formula" />

整理为标准代数方程：

<img src="https://latex.codecogs.com/png.image?a_P%20T_P%5E%7Bn%2B1%7D%20%3D%20a_W%20T_W%5E%7Bn%2B1%7D%20%2B%20a_E%20T_E%5E%7Bn%2B1%7D%20%2B%20b_P" alt="formula" />

其中 <img src="https://latex.codecogs.com/png.image?a_W%20%3D%20D_w%20%2B%20%5Cmax%28F_w%2C%200%29" alt="formula" />，<img src="https://latex.codecogs.com/png.image?a_E%20%3D%20D_e%20%2B%20%5Cmax%28-F_e%2C%200%29" alt="formula" />，<img src="https://latex.codecogs.com/png.image?a_P%20%3D%20a_P%5E0%20%2B%20a_W%20%2B%20a_E" alt="formula" />，<img src="https://latex.codecogs.com/png.image?a_P%5E0%20%3D%20C_%7B%5Cmathrm%7Bapp%7D%2CP%7D%20V_P%20%2F%20%5CDelta%20t" alt="formula" />。

界面导热系数采用调和平均：<img src="https://latex.codecogs.com/png.image?D_e%20%3D%20%5Cfrac%7B2%20k_P%20k_E%7D%7Bk_P%20%2B%20k_E%7D%20%5Cfrac%7BA_e%7D%7B%5Cdelta_e%7D" alt="formula" />

### 4. 边界条件

| 类型 | a_P 修正 | b 修正 |
|------|----------|--------|
| 定温 (Dirichlet) | <img src="https://latex.codecogs.com/png.image?a_P%20%5Cmathrel%7B%2B%7D%3D%20D_b" alt="formula" /> | <img src="https://latex.codecogs.com/png.image?b%20%5Cmathrel%7B%2B%7D%3D%20D_b%20T_%5Cmathrm%7Bwall%7D" alt="formula" /> |
| 热流 (Neumann) | 不变 | <img src="https://latex.codecogs.com/png.image?b%20%5Cmathrel%7B%2B%7D%3D%20q_%5Cmathrm%7Bin%7D%20A_b" alt="formula" /> |
| 对流 (Robin) | <img src="https://latex.codecogs.com/png.image?a_P%20%5Cmathrel%7B%2B%7D%3D%20h_%5Cmathrm%7Beff%7D%20A_b" alt="formula" /> | <img src="https://latex.codecogs.com/png.image?b%20%5Cmathrel%7B%2B%7D%3D%20h_%5Cmathrm%7Beff%7D%20A_b%20T_%5Cinfty" alt="formula" /> |

其中 <img src="https://latex.codecogs.com/png.image?h_%5Cmathrm%7Beff%7D%20%3D%201%2F%28d_b%2Fk_b%20%2B%201%2Fh%29" alt="formula" />。定温边界是 <img src="https://latex.codecogs.com/png.image?h%20%5Cto%20%5Cinfty" alt="formula" /> 的极限。

### 5. 非线性求解

由于 <img src="https://latex.codecogs.com/png.image?C_%5Cmathrm%7Bapp%7D" alt="formula" />、<img src="https://latex.codecogs.com/png.image?k" alt="formula" />、<img src="https://latex.codecogs.com/png.image?rho" alt="formula" /> 均随温度变化，每一时间步需 Picard 迭代 + 欠松弛：

1. 取 <img src="https://latex.codecogs.com/png.image?T%5E%7Bn%2B1%2C0%7D%20%3D%20T%5En" alt="formula" />
2. 根据当前温度更新液相率及物性
3. 组装系数矩阵，TDMA（1D）/ 逐线 TDMA+SOR（2D）求解
4. 温度和物性欠松弛：<img src="https://latex.codecogs.com/png.image?T%20%3D%20%281-%5Calpha%29T_%5Cmathrm%7Bprev%7D%20%2B%20%5Calpha%20T_%5Cmathrm%7Bsolved%7D" alt="formula" />
5. 检验收敛：<img src="https://latex.codecogs.com/png.image?%5Cmax%7CT%5E%7Bm%2B1%7D%20-%20T%5Em%7C%20%3C%20%5Cvarepsilon_T" alt="formula" />

---

## 数值格式

| 项目 | 方案 | 说明 |
|------|------|------|
| 空间离散 | 有限体积法（FVM） | 守恒性好，适合不均匀物性 |
| 时间推进 | 全隐式（一阶） | 无条件稳定，适合相变问题 |
| 对流项 | 一阶迎风 | 稳定，可升级为 QUICK |
| 扩散项 | 中心差分 | 二阶精度 |
| 界面导热系数 | 调和平均 | 保证冰-水界面热流连续 |
| 相变处理 | 表观热容法 | 无需追踪相界面 |
| 线性求解 | TDMA (1D) / 逐线 TDMA+SOR (2D) | |
| 非线性处理 | Picard 迭代 + 欠松弛 | 抑制相变区震荡 |

---

## 模块结构

```
ConvectTherm/
├── main.py            主程序：6 个算例 + 可视化 + 数据存储
├── solver/
│   ├── constants.py   物理常数：水/冰物性，相变参数
│   ├── physics.py     相变模型：液相率、表观热容、物性插值
│   ├── tdma.py        TDMA 三对角求解器
│   ├── solver_1d.py   一维隐式 FVM 求解器
│   ├── solver_2d.py   二维隐式 FVM 求解器
│   └── visualize.py   可视化 + 数据存储
└── output/            输出目录（图片 + .npz 数据）
```

---

## 算例设计

### 物理场景

- 计算域：10 cm × 10 cm 矩形
- 冰块：位于域中央 3–7 cm，初始 4 cm × 4 cm
- 初始温度：冰块 268.15 K (−5°C)，外围水 275.15 K (2°C)
- 物理时间：1 小时
- 纯导热（无对流），<img src="https://latex.codecogs.com/png.image?u%20%3D%20v%20%3D%200" alt="formula" />

### 算例矩阵

| Case | 维度 | 左/右边界 | 上/下边界 | 物理意义 |
|------|------|-----------|-----------|----------|
| 1 | 1D | 定温 293.15 K (两侧) | — | 恒温水浴 |
| 2 | 1D | 热流 q=1000 W/m² (左) + 定温 (右) | — | 单侧加热 |
| 3 | 1D | 对流 h=200 (左) + 定温 (右) | — | 对流+恒温 |
| 4 | 2D | 定温 293.15 K | 定温 293.15 K | 四面恒温 |
| 5 | 2D | 热流 q=1000 W/m² | 定温 293.15 K | 侧面加热 |
| 6 | 2D | 对流 h=200 | 定温 293.15 K | 侧面对流 |

---

## 计算结果

### 一维算例

Case 1冰块中心温度随时间变化：

<img src="output/case1_1d_dirichlet_Tc.png" alt="case1_1d_dirichlet_Tc" style="zoom:33%;" />

Case 1 最终温度分布：

<img src="output/case1_1d_profile.png" alt="case1_1d_profile" style="zoom: 33%;" />

Case 2 冰块中心温度随时间变化：

<img src="output/case2_1d_neumann_Tc.png" alt="case2_1d_neumann_Tc" style="zoom:33%;" />

Case 2 最终温度分布：

<img src="output/case2_1d_profile.png" alt="case2_1d_profile" style="zoom:33%;" />

Case 3 冰块中心温度随时间变化：

<img src="output/case3_1d_robin_Tc.png" alt="case3_1d_robin_Tc" style="zoom:33%;" />

Case 3 最终温度分布：

<img src="output/case3_1d_profile.png" alt="case3_1d_profile" style="zoom:33%;" />

### 二维算例

Case 4 温度场云图：

<img src="output/case4_2d_T.png" alt="case4_2d_T" style="zoom:33%;" />



Case 4 液相率分布：

<img src="output/case4_2d_fl.png" alt="case4_2d_fl" style="zoom:33%;" />

Case 5 温度场云图：

<img src="output/case5_2d_T.png" alt="case5_2d_T" style="zoom:33%;" />

Case 5 液相率分布：

<img src="output/case5_2d_fl.png" alt="case5_2d_fl" style="zoom:33%;" />

Case 6 温度场云图：



<img src="output/case6_2d_T.png" alt="case6_2d_T" style="zoom:33%;" />

Case 6 液相率分布：

<img src="output/case6_2d_fl.png" alt="case6_2d_fl" style="zoom:33%;" />

---

## 运行方式

```bash
python main.py
```

运行约 10 分钟，结果保存在 `output/` 目录（PNG 图片 + .npz 数据文件）。

### 环境依赖

- Python 3.10+
- NumPy, SciPy, Matplotlib

---

## 物理参数

| 参数 | 符号 | 数值 | 单位 |
|------|------|------|------|
| 水密度 | <img src="https://latex.codecogs.com/png.image?%5Crho_w" alt="formula" /> | 999.8 | kg/m³ |
| 水导热系数 | <img src="https://latex.codecogs.com/png.image?k_w" alt="formula" /> | 0.569 | W/(m·K) |
| 水比热容 | <img src="https://latex.codecogs.com/png.image?c_%7Bp%2Cw%7D" alt="formula" /> | 4217 | J/(kg·K) |
| 冰密度 | <img src="https://latex.codecogs.com/png.image?%5Crho_i" alt="formula" /> | 917.0 | kg/m³ |
| 冰导热系数 | <img src="https://latex.codecogs.com/png.image?k_i" alt="formula" /> | 2.25 | W/(m·K) |
| 冰比热容 | <img src="https://latex.codecogs.com/png.image?c_%7Bp%2Ci%7D" alt="formula" /> | 2090 | J/(kg·K) |
| 熔点 | <img src="https://latex.codecogs.com/png.image?T_m" alt="formula" /> | 273.15 | K |
| 相变区间半宽 | <img src="https://latex.codecogs.com/png.image?%5CDelta%20T_m%2F2" alt="formula" /> | 0.25 | K |
| 相变潜热 | <img src="https://latex.codecogs.com/png.image?L" alt="formula" /> | 334000 | J/kg |
