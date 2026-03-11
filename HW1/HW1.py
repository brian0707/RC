# -*- coding: utf-8 -*-
"""
HW1: 鋼筋混凝土柱斷面分析 (OpenSeesPy)
斷面: 60cm x 60cm
主筋: 16支 #8
箍筋: #4
單軸向力行為分析

圍束效應依 Mander et al. (1988) "Theoretical Stress-Strain Model for Confined Concrete"
J. Struct. Eng., ASCE, 114(8), 1804-1826.
"""

import openseespy.opensees as ops
import numpy as np
import os
from pathlib import Path

# 輸出目錄
RESULTS_DIR = Path(__file__).resolve().parent / 'Results'
RESULTS_DIR.mkdir(exist_ok=True)

# ==================== 設計參數 ====================
# 單位: N, mm, MPa

# 斷面幾何
B = 600          # 斷面寬度 (mm)
H = 600          # 斷面深度 (mm)
cover = 50       # 保護層厚度: 混凝土表面至箍筋外緣 (mm)

# 鋼筋規格
# #8: 直徑 25.4mm, 面積 506.7 mm²
# #4: 直徑 12.7mm, 面積 126.7 mm²
d_bar_main = 25.4
A_bar_main = np.pi * (d_bar_main / 2) ** 2   # 506.7 mm²
d_bar_stirrup = 12.7
A_bar_stirrup = np.pi * (d_bar_stirrup / 2) ** 2  # 126.7 mm²
s_stirrup = 100   # 箍筋間距 (mm)

# 主筋中心至混凝土表面距離
# = 保護層 + 箍筋半徑 + 主筋半徑
dist_to_bar_center = cover + (d_bar_stirrup / 2) + (d_bar_main / 2)  # 69 mm

# 主筋配置: 4x4 網格, 最外排中心距斷面邊緣 69mm
# 有效範圍: 600 - 2*69 = 462 mm, 間距 462/3 = 154 mm
# 主筋座標 (以斷面中心為原點): ±77, ±231 mm
bar_positions = []
for y in [-231, -77, 77, 231]:
    for z in [-231, -77, 77, 231]:
        bar_positions.append((y, z))

# 纖維細分數
n_fib_y = 12
n_fib_z = 12
n_fib_cover = 2   # 保護層纖維數 (較薄)

# 材料參數
fpc = -28       # 混凝土抗壓強度 (MPa, 負值)
epsc0 = 0.002   # 達最大強度時應變
fpcu = -5.6     # 壓碎強度 (MPa, 負值)
epsu = 0.005    # 壓碎時應變

Fy = 420        # 鋼筋降伏強度 (MPa)
fyh = 420       # 箍筋降伏強度 (MPa)
E0 = 200000     # 鋼筋彈性模數 (MPa)
b_steel = 0.01  # 應變硬化比

# ==================== 圍束效應 (Mander et al. 1988) ====================
# 核心區尺寸: 箍筋中心線至中心線 (confined region = inside hoops)
cover_to_stirrup_center = cover + (d_bar_stirrup / 2)  # 56.35 mm
bc = B - 2 * cover_to_stirrup_center   # 核心寬度 (mm)
dc = H - 2 * cover_to_stirrup_center   # 核心深度 (mm)
core_bound = (bc + dc) / 4   # 核心邊界距斷面中心 = 243.65 mm

# 箍筋體積比 ρs (Mander Eq. 矩形斷面)
# Ash = 單向箍筋總斷面積 (每側箍筋腿數 × 單支面積)
# 矩形箍筋: 每向 2 支 #4 抵抗側向力
Ash = 2 * A_bar_stirrup
rho_s = (2 * Ash * (bc + dc)) / (s_stirrup * bc * dc)

# 圍束有效性係數 ke (Mander Eq. 考慮拱效應與箍筋間距)
# wi' = 主筋間淨距; s' = 箍筋間距
bar_span = 2 * (H/2 - dist_to_bar_center)   # 主筋有效跨度 462 mm
bar_spacing = bar_span / 3                   # 主筋中心距 154 mm
wi_prime = bar_spacing - d_bar_main          # 主筋淨距 128.6 mm
s_prime = s_stirrup
# ke = (1 - Σwi'²/(6bc·dc)) × (1 - s'/(2bc)) × (1 - s'/(2dc))
# Σwi'²: 每邊 3 個間距, 4 邊共 12 個
sum_wi_sq = 12 * (wi_prime ** 2)
ke = (1 - sum_wi_sq / (6 * bc * dc)) * (1 - s_prime / (2 * bc)) * (1 - s_prime / (2 * dc))
ke = max(ke, 0.3)  # 下限避免不合理值

# 有效側向圍束應力 fl (MPa) - Mander Eq.
fc_abs = abs(fpc)
fl = 0.5 * ke * rho_s * fyh

# 圍束混凝土強度 f'cc - Mander 極限強度準則 (Popovics 型)
# f'cc/f'c = -1.254 + 2.254√(1 + 7.94·fl/f'c) - 2·fl/f'c
fcc_fc_ratio = -1.254 + 2.254 * np.sqrt(1 + 7.94 * fl / fc_abs) - 2 * fl / fc_abs
fcc_fc_ratio = max(fcc_fc_ratio, 1.0)
fpc_confined = fpc * fcc_fc_ratio

# 圍束混凝土峰值應變 εcc - Mander Eq.
# εcc = εc·[1 + 5(f'cc/f'c - 1)]
epsc0_confined = epsc0 * (1 + 5 * (fcc_fc_ratio - 1))

# 圍束混凝土極限應變 εcu - 受箍筋斷裂限制 (Mander: 0.01~0.06)
# NZSEE C5 修正: εcu = 0.004 + 0.6·ρs·fyh·εsu/f'cc, εsu≈0.09
eps_su = 0.09   # 鋼筋極限拉應變
epsu_confined = 0.004 + 0.6 * rho_s * fyh * eps_su / abs(fpc_confined)
epsu_confined = min(max(epsu_confined, 0.01), 0.06)  # 限制於 0.01~0.06

# 殘餘強度 fpcu (壓碎時應力)
fpcu_confined = 0.2 * fpc_confined

# 輸出圍束參數
print("=== Mander et al. (1988) 圍束參數 ===")
print(f"  ke = {ke:.4f}, ρs = {rho_s:.4f}, fl = {fl:.2f} MPa")
print(f"  f'cc/f'c = {fcc_fc_ratio:.3f}, εcc = {epsc0_confined:.4f}, εcu = {epsu_confined:.4f}")

# ==================== 建立模型 ====================
ops.wipe()
ops.model('basic', '-ndm', 2, '-ndf', 3)

# 材料定義
# 1: 未圍束混凝土 (保護層) - Concrete01
# 2: 圍束混凝土 (核心) - Concrete02 較佳表示 Mander 壓碎後行為
# 3: 鋼筋
ops.uniaxialMaterial('Concrete01', 1, fpc, epsc0, fpcu, epsu)
# Concrete02: lambda=0.1 表示卸載勁度較低 (圍束混凝土韌性較高)
ops.uniaxialMaterial('Concrete02', 2, fpc_confined, epsc0_confined, fpcu_confined, epsu_confined, 0.1, 0.0, 0.0)
ops.uniaxialMaterial('Steel01', 3, Fy, E0, b_steel)

# 纖維斷面定義 (分區: 保護層 + 核心)
# 頂點 I: 左下, 頂點 J: 右上 (y, z 以斷面中心為原點)
# 核心邊界 = 箍筋中心線 (core_bound 已於 Mander 區塊計算)

ops.section('Fiber', 1, '-GJ', 1e10)
# 核心區 (圍束混凝土, 箍筋內側)
ops.patch('rect', 2, n_fib_y, n_fib_z, -core_bound, -core_bound, core_bound, core_bound)
# 保護層 (箍筋外側至混凝土表面)
ops.patch('rect', 1, n_fib_cover, n_fib_z, -H/2, -B/2, -core_bound, B/2)   # 左
ops.patch('rect', 1, n_fib_cover, n_fib_z, core_bound, -B/2, H/2, B/2)     # 右
ops.patch('rect', 1, n_fib_y, n_fib_cover, -H/2, -B/2, H/2, -core_bound)   # 下
ops.patch('rect', 1, n_fib_y, n_fib_cover, -H/2, core_bound, H/2, B/2)    # 上

# 主筋 (16 支 #8)
for y, z in bar_positions:
    ops.fiber(y, z, A_bar_main, 3)

# ==================== 節點與元素 ====================
ops.node(1, 0.0, 0.0)
ops.node(2, 0.0, 0.0)

# zeroLengthSection: 兩節點同位置, 以斷面連接
# -orient: 使斷面軸向對齊全域 Y (柱軸向)
ops.element('zeroLengthSection', 1, 1, 2, 1, '-orient', 0, 1, -1, 0)

# ==================== 邊界條件 ====================
ops.fix(1, 1, 1, 1)   # 節點1 完全固定
ops.fix(2, 1, 0, 1)   # 節點2 僅 Y 向自由 (軸向)

# ==================== 分析設定 ====================
ops.timeSeries('Linear', 1)
ops.pattern('Plain', 1, 1)
ops.load(2, 0, -1, 0)  # 參考軸力 (向下為壓)

ops.constraints('Plain')
ops.numberer('RCM')
ops.system('BandGen')
ops.test('NormDispIncr', 1e-8, 50)
ops.algorithm('Newton')
ops.integrator('LoadControl', 0.0)
ops.analysis('Static')

# 施加初始軸力 (可選, 此處設為 0 開始)
ops.analyze(1)

# ==================== 位移控制分析 (單軸向力) ====================
ops.wipeAnalysis()
ops.loadConst('-time', 0.0)

# 位移控制: 軸向壓縮
ops.timeSeries('Linear', 2)
ops.pattern('Plain', 2, 2)
ops.load(2, 0, -1, 0)

ops.constraints('Plain')
ops.numberer('RCM')
ops.system('BandGen')
ops.test('NormDispIncr', 1e-8, 100)
ops.algorithm('Newton')
# 位移增量: 較小步長以捕捉峰值後軟化行為 (Mander 圍束混凝土)
disp_incr = -0.005   # mm/step
ops.integrator('DisplacementControl', 2, 2, disp_incr)
ops.analysis('Static')

# 記錄: 時間 | 軸力 P | 彎矩 Mz (切換至 Results 目錄確保 OpenSeesPy 正確寫入)
output_txt = RESULTS_DIR / 'axial_force_disp.txt'
orig_cwd = os.getcwd()
os.chdir(RESULTS_DIR)
ops.recorder('Element', '-file', 'axial_force_disp.txt', '-time', '-ele', 1, 'section', 'force')

# 執行分析: 目標位移依 εcu 估算 (柱高假設 3m → 應變 0.06 對應 180mm)
# 保守取 -8mm 以涵蓋彈性、降伏、軟化階段
max_steps = 2000
target_disp = -8.0   # 目標軸向位移 (mm)

print("開始單軸向力分析 (位移控制)...")
for i in range(max_steps):
    ok = ops.analyze(1)
    if ok != 0:
        # 收斂失敗時嘗試 ModifiedNewton
        ops.algorithm('ModifiedNewton')
        ok = ops.analyze(1)
        ops.algorithm('Newton')
        if ok != 0:
            print(f"分析收斂失敗於步驟 {i}, 位移 = {ops.nodeDisp(2, 2):.4f} mm")
            break
    disp = ops.nodeDisp(2, 2)
    if disp < target_disp:
        break

os.chdir(orig_cwd)  # 還原工作目錄
print(f"分析完成。結果輸出至 {output_txt}")
print("格式: 時間(累積位移) 軸力P 彎矩Mz")

# 繪製軸力-位移曲線
try:
    import matplotlib.pyplot as plt
    data = np.loadtxt(output_txt)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    disp = -data[:, 0]  # 位移取正值 (壓縮)
    force = -data[:, 1] if data.shape[1] > 1 else -data[:, 0]  # 軸力 (壓力為正)
    plt.figure(figsize=(8, 5))
    plt.plot(disp, force, 'b-', linewidth=1.5)
    plt.xlabel('Disp (mm)')
    plt.ylabel('P (N)')
    plt.title('Axial behavior of column')
    plt.grid(True)
    plt.tight_layout()
    output_png = RESULTS_DIR / 'axial_force_disp.png'
    plt.savefig(str(output_png), dpi=150)
    plt.show()
    print(f"圖表已儲存至 {output_png}")
except Exception as e:
    print(f"繪圖跳過: {e}")
