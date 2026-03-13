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
# 主筋配置: 16支全部配置於邊緣 (每邊 5 支含角隅)
# 有效範圍: 600 - 2*69 = 462 mm, 每邊 5 支間距 462/4 = 115.5 mm
# 主筋座標 (以斷面中心為原點): 僅邊緣 ±231, ±115.5, 0
bar_span = 2 * (H/2 - dist_to_bar_center)   # 462 mm
bar_spacing_edge = bar_span / 4   # 115.5 mm (5 bars per edge, 4 gaps)
outer_edge = H/2 - dist_to_bar_center   # 231 mm
edge_coords = [-outer_edge, -outer_edge + bar_spacing_edge, 0, outer_edge - bar_spacing_edge, outer_edge]
mid_coords = [-outer_edge + bar_spacing_edge, 0, outer_edge - bar_spacing_edge]  # excl. corners
bar_positions = []
# Top edge (5 bars)
for y in edge_coords:
    bar_positions.append((y, outer_edge))
# Right edge (3 bars, excl. corners)
for z in mid_coords:
    bar_positions.append((outer_edge, z))
# Bottom edge (5 bars, incl. corners)
for y in edge_coords:
    bar_positions.append((y, -outer_edge))
# Left edge (3 bars, excl. corners)
for z in mid_coords:
    bar_positions.append((-outer_edge, z))
# 纖維細分數
n_fib_y = 12
n_fib_z = 12
n_fib_cover = 2   # 保護層纖維數 (較薄)
# 分析組合: (fpc MPa, fy MPa)
COMBINATIONS = [
    (28, 420),   # 1. 一般混凝土 + 一般鋼筋
    (70, 420),   # 2. 高強度混凝土 + 一般鋼筋
    (28, 550),   # 3. 一般混凝土 + 高強度鋼筋
    (70, 550),   # 4. 高強度混凝土 + 高強度鋼筋
]
E0 = 200000     # 鋼筋彈性模數 (MPa)
b_steel = 0.01  # 應變硬化比
epsu = 0.005    # 未圍束混凝土壓碎應變
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
# 圍束有效性係數 ke (Mander Eq.)
bar_spacing = bar_spacing_edge
wi_prime = bar_spacing - d_bar_main
s_prime = s_stirrup
sum_wi_sq = 16 * (wi_prime ** 2)
ke = (1 - sum_wi_sq / (6 * bc * dc)) * (1 - s_prime / (2 * bc)) * (1 - s_prime / (2 * dc))
ke = max(ke, 0.3)
# 分析參數
disp_incr = -0.0000005   # mm/step
max_steps = 30000
target_disp = -8.0  # 目標軸向位移 (mm)
def run_analysis(fpc_val, fy_val, confined=True):
    """執行單一組合分析，回傳 (disp, force) 陣列。confined=True 考慮圍束效應，False 則全斷面用未圍束混凝土。"""
    fpc = -fpc_val
    fyh = fy_val
    epsc0 = 0.002   # 未圍束混凝土峰值應變 (Mander: 與強度無關)
    fpcu = 0.2 * fpc
    # 建立模型
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)
    ops.uniaxialMaterial('Concrete01', 1, fpc, epsc0, fpcu, epsu)
    if confined:
        fc_abs = abs(fpc)
        fl = 0.5 * ke * rho_s * fyh
        fcc_fc_ratio = -1.254 + 2.254 * np.sqrt(1 + 7.94 * fl / fc_abs) - 2 * fl / fc_abs
        fcc_fc_ratio = max(fcc_fc_ratio, 1.0)
        fpc_confined = fpc * fcc_fc_ratio
        epsc0_confined = epsc0 * (1 + 5 * (fcc_fc_ratio - 1))
        eps_su = 0.09
        epsu_confined = 0.004 + 0.6 * rho_s * fyh * eps_su / abs(fpc_confined)
        epsu_confined = min(max(epsu_confined, 0.01), 0.06)
        fpcu_confined = 0.2 * fpc_confined
        ops.uniaxialMaterial('Concrete02', 2, fpc_confined, epsc0_confined, fpcu_confined, epsu_confined, 0.1, 0.0, 0.0)
    ops.uniaxialMaterial('Steel01', 3, fy_val, E0, b_steel)
    ops.section('Fiber', 1, '-GJ', 1e10)
    core_mat = 2 if confined else 1
    ops.patch('rect', core_mat, n_fib_y, n_fib_z, -core_bound, -core_bound, core_bound, core_bound)
    ops.patch('rect', 1, n_fib_cover, n_fib_z, -H/2, -B/2, -core_bound, B/2)
    ops.patch('rect', 1, n_fib_cover, n_fib_z, core_bound, -B/2, H/2, B/2)
    ops.patch('rect', 1, n_fib_y, n_fib_cover, -H/2, -B/2, H/2, -core_bound)
    ops.patch('rect', 1, n_fib_y, n_fib_cover, -H/2, core_bound, H/2, B/2)
    for y, z in bar_positions:
        ops.fiber(y, z, A_bar_main, 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, 0.0, 0.0)
    ops.element('zeroLengthSection', 1, 1, 2, 1)
    ops.fix(1, 1, 1, 1)
    ops.fix(2, 0, 1, 1)
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    ops.load(2, -1, 0, 0)
    ops.constraints('Plain')
    ops.numberer('RCM')
    ops.system('BandGen')
    ops.test('NormDispIncr', 1e-8, 50)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 0.0)
    ops.analysis('Static')
    ops.analyze(1)
    ops.wipeAnalysis()
    ops.loadConst('-time', 0.0)
    ops.timeSeries('Linear', 2)
    ops.pattern('Plain', 2, 2)
    ops.load(2, -1, 0, 0)
    ops.constraints('Plain')
    ops.numberer('RCM')
    ops.system('BandGen')
    ops.test('NormDispIncr', 1e-8, 100)
    ops.algorithm('Newton')
    ops.integrator('DisplacementControl', 2, 1, disp_incr)
    ops.analysis('Static')
    suffix = 'confined' if confined else 'unconfined'
    rec_file = f'axial_force_disp_fpc{fpc_val}_fy{fy_val}_{suffix}.txt'
    orig_cwd = os.getcwd()
    os.chdir(RESULTS_DIR)
    ops.recorder('Element', '-file', rec_file, '-time', '-ele', 1, 'section', 'force')
    for _ in range(max_steps):
        ok = ops.analyze(1)
        if ok != 0:
            ops.algorithm('ModifiedNewton')
            ok = ops.analyze(1)
            ops.algorithm('Newton')
            if ok != 0:
                break
        if ops.nodeDisp(2, 1) < target_disp:
            break
    os.chdir(orig_cwd)
    data = np.loadtxt(RESULTS_DIR / rec_file)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    n_steps = len(data)
    disp = (np.arange(n_steps) + 1) * abs(disp_incr)
    force = -data[:, 1] if data.shape[1] > 1 else -data[:, 0]
    return disp, force
# 執行四種組合分析 (含圍束與未圍束)
results = []
for fpc_val, fy_val in COMBINATIONS:
    print(f"分析中: fpc={fpc_val} MPa, fy={fy_val} MPa (圍束)")
    disp_c, force_c = run_analysis(fpc_val, fy_val, confined=True)
    print(f"分析中: fpc={fpc_val} MPa, fy={fy_val} MPa (未圍束)")
    disp_u, force_u = run_analysis(fpc_val, fy_val, confined=False)
    label = 'f\'c={} MPa, fy={} MPa'.format(fpc_val, fy_val)
    results.append((label, disp_c, force_c, disp_u, force_u, fpc_val, fy_val))
# ==================== 繪圖 ====================
try:
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle, Circle
    from matplotlib.lines import Line2D
    # 配色與圖例
    colors = ['#2E86AB', '#E94F37', '#44AF69', '#F4A261']
    # --- 軸力-位移曲線 (三張圖：僅圍束、僅未圍束、兩者合併) ---
    legend_handles = [Line2D([0], [0], color=colors[i], lw=2, label=label) for i, (label, *_) in enumerate(results)]
    legend_handles_full = legend_handles + [
        Line2D([0], [0], color='gray', lw=2, linestyle='-', label='confined'),
        Line2D([0], [0], color='gray', lw=2, linestyle='--', label='unconfined')]
    def plot_common(ax):
        ax.set_xlabel('Axial Displacement (mm)')
        ax.set_ylabel('Axial Force (kN)')
        ax.grid(True)
        ax.set_xlim(0, 0.015)
    # 1. 僅圍束 (confined only)
    fig1, ax1 = plt.subplots(1, 1, figsize=(12, 7))
    for i, (label, disp_c, force_c, disp_u, force_u, fpc_val, fy_val) in enumerate(results):
        ax1.plot(disp_c, force_c / 1000, color=colors[i], linewidth=1.5, linestyle='-')
    ax1.legend(handles=legend_handles, loc='upper center', bbox_to_anchor=(0.5, -0.08), ncol=2, fontsize=9, frameon=True)
    ax1.set_title('Force-Displacement Curve (Confined Only)')
    plot_common(ax1)
    plt.tight_layout(rect=[0, 0.1, 1, 0.98])
    fig1.savefig(str(RESULTS_DIR / 'axial_force_disp_confined.png'), dpi=150, bbox_inches='tight')
    plt.show()
    print(f"Chart saved to {RESULTS_DIR / 'axial_force_disp_confined.png'}")
    # 2. 僅未圍束 (unconfined only)
    fig2, ax2 = plt.subplots(1, 1, figsize=(12, 7))
    for i, (label, disp_c, force_c, disp_u, force_u, fpc_val, fy_val) in enumerate(results):
        ax2.plot(disp_u, force_u / 1000, color=colors[i], linewidth=1.5, linestyle='--')
    ax2.legend(handles=legend_handles, loc='upper center', bbox_to_anchor=(0.5, -0.08), ncol=2, fontsize=9, frameon=True)
    ax2.set_title('Force-Displacement Curve (Unconfined Only)')
    plot_common(ax2)
    plt.tight_layout(rect=[0, 0.1, 1, 0.98])
    fig2.savefig(str(RESULTS_DIR / 'axial_force_disp_unconfined.png'), dpi=150, bbox_inches='tight')
    plt.show()
    print(f"Chart saved to {RESULTS_DIR / 'axial_force_disp_unconfined.png'}")
    # 3. 兩者合併 (confined + unconfined)
    fig3, ax3 = plt.subplots(1, 1, figsize=(12, 7))
    for i, (label, disp_c, force_c, disp_u, force_u, fpc_val, fy_val) in enumerate(results):
        ax3.plot(disp_c, force_c / 1000, color=colors[i], linewidth=1.5, linestyle='-')
        ax3.plot(disp_u, force_u / 1000, color=colors[i], linewidth=1.5, linestyle='--')
    ax3.legend(handles=legend_handles_full, loc='upper center', bbox_to_anchor=(0.5, -0.08), ncol=3, fontsize=9, frameon=True)
    ax3.set_title('Force-Displacement Curve (Confined + Unconfined)')
    plot_common(ax3)
    plt.tight_layout(rect=[0, 0.1, 1, 0.98])
    fig3.savefig(str(RESULTS_DIR / 'axial_force_disp_combined.png'), dpi=150, bbox_inches='tight')
    plt.show()
    print(f"Chart saved to {RESULTS_DIR / 'axial_force_disp_combined.png'}")
    # --- 匯出各組合 CSV (圍束與未圍束) ---
    for label, disp_c, force_c, disp_u, force_u, fpc_val, fy_val in results:
        for suffix, disp, force in [('confined', disp_c, force_c), ('unconfined', disp_u, force_u)]:
            output_csv = RESULTS_DIR / f'axial_force_disp_fpc{fpc_val}_fy{fy_val}_{suffix}.csv'
            with open(output_csv, 'w') as f:
                f.write('Displacement_mm,Axial_Force_N\n')
                for d, ff in zip(disp, force):
                    f.write(f'{d:.6e},{ff:.6e}\n')
            print(f"P-Delta data saved to {output_csv}")
    # --- 單獨斷面圖 (較大) ---
    fig_sec, ax_sec = plt.subplots(1, 1, figsize=(6, 6))
    ax_sec.set_aspect('equal')
    ax_sec.set_xlim(-320, 320)
    ax_sec.set_ylim(-320, 320)
    concrete2 = Rectangle((-B/2, -H/2), B, H, fill=True, facecolor='#E8E8E8', edgecolor='#333', linewidth=2)
    ax_sec.add_patch(concrete2)
    for y, z in bar_positions:
        ax_sec.add_patch(Circle((y, z), d_bar_main/2, fill=True, facecolor='#2C3E50', edgecolor='#1a252f', linewidth=1))
    ax_sec.plot([], [], 'o', markersize=8, color='#2C3E50', markeredgecolor='#1a252f', label='Longitudinal bars #8 (16 nos.)')
    ax_sec.set_xlabel('y (mm)')
    ax_sec.set_ylabel('z (mm)')
    ax_sec.set_title('Section Layout (60x60 cm)')
    ax_sec.legend(loc='upper right')
    ax_sec.grid(True, alpha=0.3)
    plt.tight_layout()
    output_section = RESULTS_DIR / 'section.png'
    fig_sec.savefig(str(output_section), dpi=150)
    plt.show()
    print(f"Section plot saved to {output_section}")
    # --- 混凝土應力-應變曲線 (fpc=28 與 fpc=70，實線=confined、虛線=unconfined) ---
    def concrete01_stress_strain(eps, fpc, epsc0, fpcu, epsu):
        """Concrete01/KSP: ascending parabolic, descending linear."""
        sig = np.zeros_like(eps)
        for i, e in enumerate(eps):
            if e <= 0:
                sig[i] = 0
            elif e <= epsc0:
                sig[i] = fpc * (2*e/epsc0 - (e/epsc0)**2)
            elif e <= epsu:
                sig[i] = fpc - (fpc - fpcu) * (e - epsc0) / (epsu - epsc0)
            else:
                sig[i] = fpcu
        return sig
    def mander_confined_params(fpc_val, fy_val):
        fpc = -fpc_val
        epsc0 = 0.002
        fl = 0.5 * ke * rho_s * fy_val
        fcc_ratio = max(-1.254 + 2.254 * np.sqrt(1 + 7.94 * fl / fpc_val) - 2 * fl / fpc_val, 1.0)
        fpc_conf = fpc * fcc_ratio
        epsc0_conf = epsc0 * (1 + 5 * (fcc_ratio - 1))
        epsu_conf = 0.004 + 0.6 * rho_s * fy_val * 0.09 / abs(fpc_conf)
        epsu_conf = min(max(epsu_conf, 0.01), 0.06)
        fpcu_conf = 0.2 * fpc_conf
        return fpc_conf, epsc0_conf, fpcu_conf, epsu_conf
    conc_colors = ['#2E86AB', '#E94F37']
    conc_fpc_list = [28, 70]
    fy_for_conf = 420
    n_pts = 200
    eps_max = 0.02   # 應變最大值 2%
    strain = np.linspace(0, eps_max, n_pts)
    strain_pct = strain * 100
    fig3, ax3 = plt.subplots(1, 1, figsize=(10, 6))
    conc_data = {}
    for i, fpc_val in enumerate(conc_fpc_list):
        epsc0 = 0.002
        fpcu = 0.2 * (-fpc_val)
        sig_un = concrete01_stress_strain(strain, -fpc_val, epsc0, fpcu, epsu)
        fpc_c, epsc0_c, fpcu_c, epsu_c = mander_confined_params(fpc_val, fy_for_conf)
        sig_conf = concrete01_stress_strain(strain, fpc_c, epsc0_c, fpcu_c, epsu_c)
        conc_data[fpc_val] = (sig_un, sig_conf)
        ax3.plot(strain_pct, -np.array(sig_un), color=conc_colors[i], linewidth=1.5, linestyle='--')
        ax3.plot(strain_pct, -np.array(sig_conf), color=conc_colors[i], linewidth=1.5, linestyle='-')
    legend_handles_conc = [Line2D([0], [0], color=conc_colors[i], lw=2, label=f"f'c={fpc_val} MPa") for i, fpc_val in enumerate(conc_fpc_list)]
    legend_handles_conc.append(Line2D([0], [0], color='gray', lw=2, linestyle='-', label='confined'))
    legend_handles_conc.append(Line2D([0], [0], color='gray', lw=2, linestyle='--', label='unconfined'))
    ax3.legend(handles=legend_handles_conc, loc='upper right', ncol=1, fontsize=9, frameon=True)
    ax3.set_xlabel('Strain (%)')
    ax3.set_ylabel('Stress (MPa)')
    ax3.set_title('Concrete Stress-Strain Curves')
    ax3.grid(True)
    plt.tight_layout()
    output_conc = RESULTS_DIR / 'concrete_stress_strain.png'
    fig3.savefig(str(output_conc), dpi=150, bbox_inches='tight')
    plt.show()
    print(f"Concrete stress-strain plot saved to {output_conc}")
    # Export to CSV (fpc=28 與 fpc=70)
    output_csv = RESULTS_DIR / 'concrete_stress_strain.csv'
    with open(output_csv, 'w') as f:
        f.write('Strain_pct,Stress_fpc28_Unconfined,Stress_fpc28_Confined,Stress_fpc70_Unconfined,Stress_fpc70_Confined\n')
        s28u, s28c = conc_data[28][0], conc_data[28][1]
        s70u, s70c = conc_data[70][0], conc_data[70][1]
        for j in range(n_pts):
            f.write(f'{strain[j]*100:.6e},{s28u[j]:.6e},{s28c[j]:.6e},{s70u[j]:.6e},{s70c[j]:.6e}\n')
    print(f"Concrete stress-strain data saved to {output_csv}")
    # 印出混凝土參數 (fpc=28 與 fpc=70 為例)
    for fpc_val in [28, 70]:
        epsc0 = 0.002
        fpcu = 0.2 * (-fpc_val)
        fpc_c, epsc0_c, fpcu_c, epsu_c = mander_confined_params(fpc_val, 420)
        print(f"--- fpc={fpc_val} MPa ---")
        print(f"  epsc0(unconfined)={epsc0}, epsc0_confined={epsc0_c}")
        print(f"  epsu(unconfined)={epsu}, epsu_confined={epsu_c}")
        print(f"  fpcu(unconfined)={fpcu}, fpcu_confined={fpcu_c}")
        print(f"  fpc(unconfined)={-fpc_val}, fpc_confined={fpc_c}")
except Exception as e:
    print(f"Plotting skipped: {e}")