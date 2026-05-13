"""物理常数与默认参数。"""

# 水 (液相)
RHO_W = 999.8       # kg/m³
K_W = 0.569         # W/(m·K)
CP_W = 4217.0       # J/(kg·K)

# 冰 (固相)
RHO_I = 917.0       # kg/m³
K_I = 2.25          # W/(m·K)
CP_I = 2090.0       # J/(kg·K)

# 相变参数
T_MELT = 273.15     # K
DELTA_T_M = 0.5     # K，相变区间宽度
L_FUSION = 334000.0 # J/kg

EPS_T = 1e-3        # K，液相率计算容差
