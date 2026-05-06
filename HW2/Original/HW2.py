# Introduction ------------------------------------------------------------------------------------------------------------------------------------------------#


# i import DLC ------------------------------------------------------------------------------------------------------------------------------------------------#
import math
import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import openseespy.opensees as ops

# 與腳本同 repo 之 HW1/… 內含 unit_SI、rebar（以 unit_SI.py 實際路徑為準，避免資料夾名或 cwd 差異）
def _ensure_vendor_path() -> None:
    here = Path(__file__).resolve().parent
    repo = here.parent  # HW2_yuhong 之上一層 = RCCHung
    candidates = [
        repo / "HW1" / "宥宏",
        here / "宥宏",
    ]
    for d in candidates:
        if (d / "unit_SI.py").is_file():
            sys.path.insert(0, str(d))
            return
    hw1 = repo / "HW1"
    if hw1.is_dir():
        for found in hw1.rglob("unit_SI.py"):
            if found.name == "unit_SI.py":
                sys.path.insert(0, str(found.parent))
                return
    raise ImportError(
        "找不到 unit_SI.py。預期位置："
        f"「{repo / 'HW1' / '宥宏'}」或腳本同層「宥宏」。"
    )


_ensure_vendor_path()

# 與腳本同層「data results」：所有匯出之數據、圖檔
_DATA_DIR = Path(__file__).resolve().parent
DATA_RESULTS_DIR = _DATA_DIR / "data results"

# ii User-define function import ------------------------------------------------------------------------------------------------------------------------------#

# Unit, Rebar
# 這章節程式碼定義了使用毫米（unit_SI.unit_SI.mm）、千牛頓（unit_SI.kN）和秒（sec）作為單位的系統。
# 包含了長度、面積、力和應力的轉換關係，並設置了一些常數，如圓周率（PI）和重力加速度（g）、極大值（Ubig）極小值（Usmall）等。
# 可用單位：
    # 長度：unit_SI.unit_SI.mm, inch, ft, cm, m
    # 面積：inch2, unit_SI.unit_SI.mm2, cm2, m2
    # 慣性矩：inch4
    # 力：kip, N, MN, tf
    # 應力：unit_SI.MPa, GPa, ksi, psi, lbf, pcf
import unit_SI
unit_SI.define_units()

# 可使用鋼筋號數：
    # D6, D10, D13, D16, D19, D22, D25, D29, D32, D36, D38, D43, D50, D57, D64
# 可使用鋼筋面積：
    # A6, A10, A13, A16, A19, A22, A25, A29, A32, A36, A38, A43, A50, A57, A64
import rebar
rebar.define_rebar()



# iii User-define function ------------------------------------------------------------------------------------------------------------------------------------#

# 計算圍束混凝土強度
def func_get_mander_params(fco, fyh, BCol, HCol, clear_cover, n_long_y, n_long_z, d_long, d_tie, s_spacing):
    """
    計算 Mander 受限混凝土參數
    fco: 未圍束混凝土強度 (MPa)
    fyh: 箍筋屈服強度 (MPa)
    B, H: 斷面寬、高 (mm)
    cover: 淨保護層厚度 (mm)
    n_long_y, n_long_z: Y軸與 Z軸向的縱筋支數
    d_long: 縱筋直徑 (mm)
    d_tie: 箍筋直徑 (mm)
    s_spacing: 箍筋垂直間距 (mm)
    """
    # 1. 計算核心尺寸 (至箍筋中心線) [cite: 132, 233]
    bc = BCol - 2 * clear_cover - d_tie
    dc = HCol - 2 * clear_cover - d_tie
    Acc = bc * dc # 核心總面積 [cite: 132, 234]
    
    # 2. 計算縱筋淨距之平方和 (w'i^2) 以求 ke [cite: 226, 231, 236]
    # 假設鋼筋均勻分佈，計算相鄰鋼筋間的淨距離
    wy = (bc-(n_long_y*d_long)) / (n_long_y-1)
    wz = (dc-(n_long_z*d_long)) / (n_long_z-1)
    sum_wi_sq = (n_long_y - 1) * 2 * (wy**2) + (n_long_z - 1) * 2 * (wz**2)
    
    # 3. 計算圍束有效係數 ke [cite: 137, 236, 577]
    s_prime = s_spacing - d_tie # 箍筋間淨距 [cite: 171, 726]
    rho_cc = ( (n_long_y*2 + (n_long_z-2)*2) * (0.25 * math.pi * d_long**2) ) / Acc # 縱筋比 [cite: 142, 143]
    
    ke = (1 - sum_wi_sq / (6 * bc * dc)) * (1 - s_prime / (2 * bc)) * (1 - s_prime / (2 * dc)) / (1 - rho_cc)
    
    # 4. 計算側向圍束應力 fl' [cite: 134, 219, 255]
    # 假設為 4 肢箍筋配置 (你的圖面配置)
    Asx = 4 * (0.25 * math.pi * d_tie**2) # X向總箍筋面積
    Asy = 4 * (0.25 * math.pi * d_tie**2) # Y向總箍筋面積
    flx = (Asx / (s_spacing * dc)) * fyh # 全域側壓力 [cite: 248, 249]
    fly = (Asy / (s_spacing * bc)) * fyh
    
    flx_prime = ke * flx # 有效側壓力 [cite: 255, 256]
    fly_prime = ke * fly
    
    # 5. 計算受限強度提升 fcc' [cite: 91, 266, 267]
    # 使用 Mander 針對不同兩向側壓力的圖解簡化公式或威廉-瓦克準則
    # 此處取兩向平均值進行強度計算 (對稱斷面常用)
    fl_avg = (flx_prime + fly_prime) / 2
    fcc = fco * (-1.254 + 2.254 * math.sqrt(1 + 7.94 * fl_avg / fco) - 2 * fl_avg / fco)
    
    # 6. 計算對應應變 ecc [cite: 51, 94, 579]
    eco = 0.002 # 未圍束峰值應變 [cite: 117]
    ecc = eco * (1 + 5 * (fcc / fco - 1))
    
    # 7. 計算形狀參數 r 與 Ec [cite: 118, 120, 581]
    fco_MPa = fco / unit_SI.MPa # 轉換為 MPa 單位
    Ec = 5000 * math.sqrt(fco_MPa)*unit_SI.MPa # MPa [cite: 120]
    Esec = fcc / ecc
    r = Ec / (Ec - Esec)
    
    return fcc, ecc, Ec, r

# 4. 定義材料
def func_MatDef(BCol, HCol, D_Longi, D_StirUp):
    # 材料定義參數
    # unconfined concrete (28MPa)
    ID28MPaUC = 4001
    fpc28MPaUC = -28.0 * unit_SI.MPa
    epsc028MPaUC = -0.002
    fpcu28MPaUC = 0.0
    epsU28MPaUC = -0.004
    ops.uniaxialMaterial("Concrete01", ID28MPaUC, fpc28MPaUC, epsc028MPaUC, fpcu28MPaUC, epsU28MPaUC)


    # unconfined concrete (70MPa)
    ID70MPaUC = 4002
    fpc70MPaUC = -70 * unit_SI.MPa
    epsc070MPaUC = -0.002
    fpcu70MPaUC = 0.0
    epsU70MPaUC = -0.004
    ops.uniaxialMaterial("Concrete01", ID70MPaUC, fpc70MPaUC, epsc070MPaUC, fpcu70MPaUC, epsU70MPaUC)


    # confined concrete (28MPa)
    ID28MPaCC = 4101
    fyh = 420.0 * unit_SI.MPa
    clear_cover = 40.0 * unit_SI.mm
    n_long_y = 5
    n_long_z = 5
    s_spacing = 100.0 * unit_SI.mm
    fcc_28, ecc_28, Ec_28, r_28 = func_get_mander_params(-fpc28MPaUC, fyh, BCol, HCol, clear_cover, n_long_y, n_long_z, D_Longi, D_StirUp, s_spacing)
    # Concrete04: fpc, epsc, fpcu, epsu, ec, r
    # ops.uniaxialMaterial('Concrete04', ID28MPaCC, -fcc_28, -ecc_28, -0.2*fcc_28, -0.02, Ec_28, r_28)
    ops.uniaxialMaterial('Concrete01', ID28MPaCC, -fcc_28, -ecc_28, -0.2*fcc_28, -0.02)
    
    # confined concrete (70MPa)
    ID70MPaCC = 4102
    fyh = 420.0 * unit_SI.MPa
    clear_cover = 40.0 * unit_SI.mm
    n_long_y = 5
    n_long_z = 5
    s_spacing = 100.0 * unit_SI.mm
    fcc_70, ecc_70, Ec_70, r_70 = func_get_mander_params(-fpc70MPaUC, fyh, BCol, HCol, clear_cover, n_long_y, n_long_z, D_Longi, D_StirUp, s_spacing)
    # ops.uniaxialMaterial('Concrete04', ID70MPaCC, -fcc_70, -ecc_70, -0.2*fcc_70, -0.015, Ec_70, r_70)
    ops.uniaxialMaterial('Concrete01', ID70MPaCC, -fcc_70, -ecc_70, -0.2*fcc_70, -0.015)

    # rebar (SD420)
    IDSD420 = 4201
    fySD420 = 420 * unit_SI.MPa
    fuSD420 = fySD420*1.25
    EsSD420 = 200 * unit_SI.GPa
    EshSD420 = EsSD420/30
    bSD420 = 0.018
    R0SD420 = 0.13
    ops.uniaxialMaterial('ReinforcingSteel', IDSD420, fySD420, fuSD420, EsSD420, EshSD420, bSD420, R0SD420)

    # rebar (SD550)
    IDSD550 = 4202
    fySD550 = 550 * unit_SI.MPa
    fuSD550 = fySD550*1.25
    EsSD550 = 200 * unit_SI.GPa
    EshSD550 = EsSD550/30
    bSD550 = 0.018
    R0SD550 = 0.13
    ops.uniaxialMaterial('ReinforcingSteel', IDSD550, fySD550, fuSD550, EsSD550, EshSD550, bSD550, R0SD550)

    return fcc_28, ecc_28, fcc_70, ecc_70

# 斷面定義參數
def func_sectionDef(IDSec, IDUC, IDCC, IDMat_Longi, BCol, HCol, DLongi, ALongi, DStirUp):
    # Y軸為橫軸，Z軸為縱軸

    # Define Section
    ops.section('Fiber', IDSec)


    # cover
    cover = 50*unit_SI.mm + DStirUp + DLongi/2

    # edge
    yE1 = -BCol/2
    yE2 = BCol/2
    zE1 = -HCol/2
    zE2 = HCol/2
    # core
    yC1 = yE1 + cover - DLongi/2
    yC2 = -yC1
    zC1 = zE1 + cover - DLongi/2
    zC2 = -zC1

    # Longitudinal
    spacing = (BCol-2*cover)/4
    yL1 = yE1 + cover + spacing*0
    yL2 = yE1 + cover + spacing*1
    yL3 = yE1 + cover + spacing*2
    yL4 = yE1 + cover + spacing*3
    yL5 = yE1 + cover + spacing*4

    zL1 = zE1 + cover + spacing*0
    zL2 = zE1 + cover + spacing*1
    zL3 = zE1 + cover + spacing*2
    zL4 = zE1 + cover + spacing*3
    zL5 = zE1 + cover + spacing*4

    # 核心混凝土
    ops.patch('rect', IDCC, 20, 20, yC1, zC1, yC2, zC2)

    # 保護層混凝土 (上, 下, 左, 下)
    ops.patch('rect', IDUC, 24, 2, yE1, zC2, yE2, zE2)
    ops.patch('rect', IDUC, 24, 2, yE1, zE1, yE2, zC1)
    ops.patch('rect', IDUC, 2, 24,  yE1, zC1, yC1, zC2)
    ops.patch('rect', IDUC, 2, 24,  yC2, zC1, yE2, zC2)

    # Longitudinal bars (左至右直欄：5, 2, 2, 2, 5)
    ops.layer('straight', IDMat_Longi, 5, ALongi, yL1, zL1, yL1, zL5)
    ops.layer('straight', IDMat_Longi, 2, ALongi, yL2, zL1, yL2, zL5)
    ops.layer('straight', IDMat_Longi, 2, ALongi, yL3, zL1, yL3, zL5)
    ops.layer('straight', IDMat_Longi, 2, ALongi, yL4, zL1, yL4, zL5)
    ops.layer('straight', IDMat_Longi, 5, ALongi, yL5, zL1, yL5, zL5)


    print("Define Fiber Section")
    print("-------------------------------------------------------------------------------------")

# 執行單調軸壓分析的獨立函數
def func_monotonic_axial_analysis(LCol, handles, caseNum):
    """
    執行單調軸壓分析的獨立函數
    """
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    ops.load(3003, 0.0, -1.0, 0.0) # 施加單位軸壓

    ops.system('BandGeneral')
    ops.numberer('Plain')
    ops.constraints('Plain')
    ops.integrator('DisplacementControl', 3003, 2, -0.04) # 每次推進 -0.04
    ops.algorithm('Newton')
    ops.analysis('Static')

    plot_data = []
    
    for i in range(1000):
        ok = ops.analyze(1)
        if ok == 0:
            # 讀取節點位移與元素力量
            force_ele1 = ops.eleResponse(6001, 'globalForce')
            force_ele2 = ops.eleResponse(6002, 'globalForce')
            
            d = [
                ops.nodeDisp(3001, 2), force_ele1[1],
                ops.nodeDisp(3002, 2), force_ele1[4],
                ops.nodeDisp(3003, 2), force_ele2[4]
            ]
            
            # 手動寫入檔案
            for j in range(6):
                handles[j].write(f"{d[j]:12.6e}\n")
            
            # 紀錄繪圖數據 [軸向應變, 軸向力量(kN)]
            plot_data.append([abs(d[4]) / LCol, abs(d[5])])
        else:
            print(f"Case {caseNum} 分析停止於第 {i} 步 (可能是壓碎或不收斂)")
            break
            
    return plot_data


def _section_moment_bending(ele_tag=1, sec_num=1):
    """2D fiber section：eleResponse 回傳 [軸力 N, 彎矩 M]（與模型單位一致：N、N·mm）。"""
    sf = ops.eleResponse(ele_tag, "section", sec_num, "force")
    return float(sf[1])


def func_moment_curvature_export(
    case_label,
    case_num,
    IDSec,
    IDUC,
    IDCC,
    IDMat_Longi,
    BCol,
    HCol,
    DLongi,
    ALongi,
    DStirUp,
    axial_load_N,
    mu=15.0,
    num_incr=100,
):
    """
    與 MomentCurvature.py 相同概念：獨立 wipe 後建立 zeroLengthSection 纖維斷面，
    先施加常軸力，再以位移控制逐步加曲率，匯出曲率–彎矩數據。

    曲率：
    - OpenSees 內部以 mm 為長度單位時，節點轉角對參考長度 1 即為 κ（1/mm）。
    - **匯出 CSV 第一欄為 κ（1/m）**：κ(1/m) = κ(1/mm) × 1000；第二欄保留 1/mm 供核對。
    - 回傳給繪圖之曲率欄位同為 **1/m**。
    彎矩：斷面力第二分量（N·mm）。

    須與柱分析分開呼叫（內部會 ops.wipe），不可與 forceBeamColumn 模型並存。
    """
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)

    func_MatDef(BCol, HCol, DLongi, DStirUp)
    func_sectionDef(IDSec, IDUC, IDCC, IDMat_Longi, BCol, HCol, DLongi, ALongi, DStirUp)

    ops.node(1, 0.0, 0.0)
    ops.node(2, 0.0, 0.0)
    ops.fix(1, 1, 1, 1)
    ops.fix(2, 0, 1, 0)
    ops.element("zeroLengthSection", 1, 1, 2, IDSec)

    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, axial_load_N, 0.0, 0.0)

    ops.integrator("LoadControl", 0.0)
    ops.system("SparseGeneral", "-piv")
    # SI 單位下軸力極大，殘差常落在 ~1e-8；1e-9 過嚴易誤報不收斂
    ops.test("NormUnbalance", 1e-7, 25)
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.algorithm("Newton")
    ops.analysis("Static")
    ops.analyze(1)

    kappa_hist = [ops.nodeDisp(2, 3)]
    moment_hist = [_section_moment_bending()]

    cover = 50 * unit_SI.mm + DStirUp + DLongi / 2
    d_eff = float(HCol - cover)
    if IDMat_Longi == 4201:
        fy_long = 420.0 * unit_SI.MPa
    else:
        fy_long = 550.0 * unit_SI.MPa
    Es_long = 200.0 * unit_SI.GPa
    epsy = fy_long / Es_long
    Ky = epsy / (0.7 * d_eff)
    max_kappa = Ky * mu
    dK = max_kappa / num_incr

    ops.timeSeries("Linear", 2)
    ops.pattern("Plain", 2, 2)
    ops.load(2, 0.0, 0.0, 1.0)
    ops.integrator("DisplacementControl", 2, 3, dK, 1, dK, dK)

    for _ in range(num_incr):
        ok = ops.analyze(1)
        if ok != 0:
            break
        kappa_hist.append(ops.nodeDisp(2, 3))
        moment_hist.append(_section_moment_bending())

    data_dir = str(DATA_RESULTS_DIR / f"HW2_case{case_num}_Data")
    os.makedirs(data_dir, exist_ok=True)
    out_path = os.path.join(data_dir, f"MomentCurvature_{case_label}.csv")
    kappa_arr = np.asarray(kappa_hist, dtype=float)
    kappa_per_m = kappa_arr * 1000.0
    header = "curvature_1_per_m,curvature_1_per_mm,moment_Nmm"
    np.savetxt(
        out_path,
        np.column_stack([kappa_per_m, kappa_arr, moment_hist]),
        delimiter=",",
        header=header,
        comments="",
    )
    print(f"已匯出 moment–curvature：{out_path}")
    return np.column_stack([kappa_per_m, moment_hist]), out_path


def _moment_nmm_to_knm(m_nmm):
    """斷面彎矩 N·mm → kN·m"""
    return np.asarray(m_nmm, dtype=float) / 1e6


def plot_moment_curvature_single(curve, label, out_png):
    """curve: (N,2) 曲率(1/m)、彎矩(N·mm)；每案一張圖。"""
    k = curve[:, 0]
    m_knm = _moment_nmm_to_knm(curve[:, 1])
    x_plot = k
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    ax.plot(x_plot, m_knm, "-", linewidth=1.4, color="#1f77b4")
    ax.set_xlabel(r"Curvature $\phi$ (m$^{-1}$)")
    ax.set_ylabel("Moment (kN·m)")
    ax.set_title(f"Moment–curvature: {label}")
    ax.grid(True, linestyle=":", alpha=0.7)
    ax.xaxis.set_major_locator(mticker.MaxNLocator(nbins=7, prune="both"))
    ax.xaxis.set_major_formatter(mticker.FormatStrFormatter("%.3f"))
    fig.tight_layout()
    fig.savefig(out_png, dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"已匯出 moment–curvature 圖：{out_png}")


def plot_moment_curvature_combined(series, out_png):
    """series: [(label, curve_ndarray), ...]，曲率欄位為 1/m；四案合圖。"""
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]
    fig, ax = plt.subplots(figsize=(9, 5.5))
    for i, (label, curve) in enumerate(series):
        if curve is None or len(curve) == 0:
            continue
        k = curve[:, 0]
        m_knm = _moment_nmm_to_knm(curve[:, 1])
        x_plot = k
        c = colors[i % len(colors)]
        ax.plot(x_plot, m_knm, "-", linewidth=1.5, color=c, label=label)
    ax.set_xlabel(r"Curvature $\phi$ (m$^{-1}$)")
    ax.set_ylabel("Moment (kN·m)")
    ax.set_title("Moment–curvature (four cases)")
    ax.grid(True, linestyle=":", alpha=0.7)
    ax.xaxis.set_major_locator(mticker.MaxNLocator(nbins=8, prune="both"))
    ax.xaxis.set_major_formatter(mticker.FormatStrFormatter("%.3f"))
    ax.legend(loc="best", fontsize=9)
    fig.tight_layout()
    fig.savefig(out_png, dpi=240, bbox_inches="tight")
    plt.close(fig)
    print(f"已匯出四案合圖：{out_png}")


# 主函數：執行柱子分析
def func_run_column_analysis(caseNum, IDSec, BCol, HCol, IDUC, IDCC, IDMat_Longi, DLongi, ALongi,  DStirUp):

    # 1. initialize system setting --------------------------------------------------------------------------------------------------------------------------------#
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)
    
    dataDir = str(DATA_RESULTS_DIR / f"HW2_case{caseNum}_Data")
    os.makedirs(dataDir, exist_ok=True)


    # 2. Geometry -------------------------------------------------------------------------------------------------------------------------------------------------#
    LCol = 3200 *unit_SI.mm
    BCol =  600 * unit_SI.mm
    HCol =  600 * unit_SI.mm

    # 3. Node -----------------------------------------------------------------------------------------------------------------------------------------------------#
    ops.node(3001, 0, 0)
    ops.node(3002, 0, LCol/2)
    ops.node(3003, 0, LCol)
    ops.fix(3001, 1, 1, 1)

    # 4. Materials ------------------------------------------------------------------------------------------------------------------------------------------------#
    fcc_28, ecc_28, fcc_70, ecc_70 = func_MatDef(BCol, HCol, DLongi, DStirUp)
        
    # 5. Section --------------------------------------------------------------------------------------------------------------------------------------------------#
    func_sectionDef(IDSec, IDUC, IDCC, IDMat_Longi, BCol, HCol, DLongi, ALongi, DStirUp)
    # 6. Elements -------------------------------------------------------------------------------------------------------------------------------------------------#
    
    ColTransfTag = 6901
    ColTransType = 'Linear'
    ops.geomTransf(ColTransType, ColTransfTag)

    numIntPts = 5
    IDIntegrationTag = 1
    maxIter = 50
    tol = 1e-7
    ops.beamIntegration('Lobatto', IDIntegrationTag, IDSec, numIntPts)

    # element('forceBeamColumn', eleTag, iNode, jNode, transfTag, integrationTag)
    ops.element('forceBeamColumn', 6001, 3001, 3002, ColTransfTag, IDIntegrationTag, 'iter', maxIter, tol)
    ops.element('forceBeamColumn', 6002, 3002, 3003, ColTransfTag, IDIntegrationTag, 'iter', maxIter, tol)

    # 7. recorder -------------------------------------------------------------------------------------------------------------------------------------------------#
    file_names = ['Disp1.out', 'Force1.out', 'Disp2.out', 'Force2.out', 'Disp3.out', 'Force3.out']
    handles = [open(os.path.join(dataDir, name), "w") for name in file_names]

    # 8. gravity analysis -----------------------------------------------------------------------------------------------------------------------------------------#

    
    # 9. Monotonic Axial Analysis----------------------------------------------------------------------------------------------------------------------------------#

    plot_data = func_monotonic_axial_analysis(LCol, handles, caseNum)
    
    for h in handles: h.close()
    print (f"Case {caseNum} analysis completed. Data saved in {dataDir}.")
    return np.array(plot_data), fcc_28, ecc_28, fcc_70, ecc_70
    



# HW2 main program --------------------------------------------------------------------------------------------------------------------------------------------#


DLongi = rebar.D25
ALongi = rebar.A25
DStirUp = rebar.D13

BCol = 600 * unit_SI.mm
HCol = 600 * unit_SI.mm
LCol = 3200 * unit_SI.mm

case1 = 1 ; case2 = 2 ; case3 = 3 ; case4 = 4
ID28MPaUC = 4001 ; ID70MPaUC = 4002
ID28MPaCC = 4101 ; ID70MPaCC = 4102
IDSD420   = 4201 ; IDSD550   = 4202
IDsec1 = 5001 ; IDsec2 = 5002 ; IDsec3 = 5003 ; IDsec4 = 5004

# --- 執行四種組合 ---
DATA_RESULTS_DIR.mkdir(parents=True, exist_ok=True)
print(f"輸出目錄：{DATA_RESULTS_DIR}")
print("Running Column Analyses for All Cases...")
PDeltaCurve1, fcc_28c1, ecc_28c1, fcc_70c1, ecc_70c1 = func_run_column_analysis(case1, IDsec1, BCol, HCol,  ID28MPaUC, ID28MPaCC, IDSD420, DLongi, ALongi, DStirUp)
PDeltaCurve2, fcc_28c2, ecc_28c2, fcc_70c2, ecc_70c2 = func_run_column_analysis(case2, IDsec2, BCol, HCol,  ID28MPaUC, ID28MPaCC, IDSD550, DLongi, ALongi, DStirUp)
PDeltaCurve3, fcc_28c3, ecc_28c3, fcc_70c3, ecc_70c3 = func_run_column_analysis(case3, IDsec3, BCol, HCol,  ID70MPaUC, ID70MPaCC, IDSD420, DLongi, ALongi, DStirUp)
PDeltaCurve4, fcc_28c4, ecc_28c4, fcc_70c4, ecc_70c4 = func_run_column_analysis(case4, IDsec4, BCol, HCol,  ID70MPaUC, ID70MPaCC, IDSD550, DLongi, ALongi, DStirUp)

# --- Moment–curvature（zeroLengthSection 斷面分析，概念同 MomentCurvature.py；與柱模型分開 wipe）---
# 軸力取該案圍束混凝土強度之 20%×斷面積（可依設計軸力修改）
Ag = BCol * HCol
_mphi_series = []

curve1, _ = func_moment_curvature_export(
    "Case1_28MPa_SD420", case1, IDsec1, ID28MPaUC, ID28MPaCC, IDSD420,
    BCol, HCol, DLongi, ALongi, DStirUp, -0.2 * abs(fcc_28c1) * Ag,
)
plot_moment_curvature_single(
    curve1,
    "Case1_28MPa_SD420",
    DATA_RESULTS_DIR / "HW2_case1_Data" / "MomentCurvature_Case1_28MPa_SD420.png",
)
_mphi_series.append(("Case1_28MPa_SD420", curve1))

curve2, _ = func_moment_curvature_export(
    "Case2_28MPa_SD550", case2, IDsec2, ID28MPaUC, ID28MPaCC, IDSD550,
    BCol, HCol, DLongi, ALongi, DStirUp, -0.2 * abs(fcc_28c2) * Ag,
)
plot_moment_curvature_single(
    curve2,
    "Case2_28MPa_SD550",
    DATA_RESULTS_DIR / "HW2_case2_Data" / "MomentCurvature_Case2_28MPa_SD550.png",
)
_mphi_series.append(("Case2_28MPa_SD550", curve2))

curve3, _ = func_moment_curvature_export(
    "Case3_70MPa_SD420", case3, IDsec3, ID70MPaUC, ID70MPaCC, IDSD420,
    BCol, HCol, DLongi, ALongi, DStirUp, -0.2 * abs(fcc_70c3) * Ag,
)
plot_moment_curvature_single(
    curve3,
    "Case3_70MPa_SD420",
    DATA_RESULTS_DIR / "HW2_case3_Data" / "MomentCurvature_Case3_70MPa_SD420.png",
)
_mphi_series.append(("Case3_70MPa_SD420", curve3))

curve4, _ = func_moment_curvature_export(
    "Case4_70MPa_SD550", case4, IDsec4, ID70MPaUC, ID70MPaCC, IDSD550,
    BCol, HCol, DLongi, ALongi, DStirUp, -0.2 * abs(fcc_70c4) * Ag,
)
plot_moment_curvature_single(
    curve4,
    "Case4_70MPa_SD550",
    DATA_RESULTS_DIR / "HW2_case4_Data" / "MomentCurvature_Case4_70MPa_SD550.png",
)
_mphi_series.append(("Case4_70MPa_SD550", curve4))

plot_moment_curvature_combined(_mphi_series, DATA_RESULTS_DIR / "MomentCurvature_all_cases.png")

results = {
        'Case1_28MPa_SD420': (PDeltaCurve1, fcc_28c1, ecc_28c1, fcc_70c1, ecc_70c1),
        'Case2_28MPa_SD550': (PDeltaCurve2, fcc_28c2, ecc_28c2, fcc_70c2, ecc_70c2),
        'Case3_70MPa_SD420': (PDeltaCurve3, fcc_28c3, ecc_28c3, fcc_70c3, ecc_70c3),
        'Case4_70MPa_SD550': (PDeltaCurve4, fcc_28c4, ecc_28c4, fcc_70c4, ecc_70c4)
    }
 
# --- 輸出四個 Case 的曲線數據到 TXT 檔 ---
print("正在輸出數據檔案...")
for label, (data, fcc28, ecc28, fcc70, ecc70) in results.items():
    if len(data) > 0:

        fcc_28MPa = fcc28 / unit_SI.MPa
        fcc_70MPa = fcc70 / unit_SI.MPa


        # 儲存為 txt 檔，第一欄是 Strain，第二欄是 Axial Load (kN)
        filename = str(DATA_RESULTS_DIR / f"PDeltaData_{label}.txt")
        header = f"f'cc28MPa = {fcc_28MPa:.2f} MPa, ecc28MPa = {ecc28:.4f}\nf'cc70MPa = {fcc_70MPa:.2f} MPa, ecc70MPa = {ecc70:.4f}\nAxial_Strain Axial_Load(kN)"
        np.savetxt(filename, data, fmt='%12.6e', header=header, comments='')
        print("************************************")
        print(f"已儲存: {filename}")
        print(f"f'cc = {fcc_28MPa:.2f} MPa, ecc = {ecc28:.4f}")
        print(f"f'cc = {fcc_70MPa:.2f} MPa, ecc = {ecc70:.4f}")

        # 儲存為 txt 檔，第一欄是 Strain，第二欄是 Axial Load (kN)
        filename = str(DATA_RESULTS_DIR / f"PData_{label}.txt")
        header = f"f'cc28MPa = {fcc_28MPa:.2f} MPa, ecc28MPa = {ecc28:.4f}\nf'cc70MPa = {fcc_70MPa:.2f} MPa, ecc70MPa = {ecc70:.4f}\nAxial_Strain Axial_Load(kN)"
        np.savetxt(filename, data, fmt='%12.6e', header=header, comments='')
        print("************************************")
        print(f"已儲存: {filename}")
        print(f"f'cc = {fcc_28MPa:.2f} MPa, ecc = {ecc28:.4f}")
        print(f"f'cc = {fcc_70MPa:.2f} MPa, ecc = {ecc70:.4f}")




# --- 繪圖 ---
plt.figure(figsize=(10,6))
for label, val in results.items():
    data = val[0]
    if len(data) > 0:
        plt.plot(data[:, 0], data[:, 1], label=label, linewidth=2)
        
plt.xlabel('Axial Strain', fontsize=12)
plt.ylabel('Axial Load (kN)', fontsize=12)
plt.title('Monotonic Axial Response of RC Columns', fontsize=14)
plt.legend()
plt.grid(True)

_pdelta_png = DATA_RESULTS_DIR / "PDelta_Curve_Result.png"
plt.savefig(_pdelta_png, dpi=300, bbox_inches="tight")
print(f"已成功輸出圖片：{_pdelta_png}")

plt.show()



