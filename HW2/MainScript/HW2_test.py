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

# unit_SI、rebar：以 unit_SI.py 實際路徑為準
def _ensure_vendor_path() -> None:
    here = Path(__file__).resolve().parent
    candidates = [
        here,
        here / "宥宏",
        here.parent / "宥宏",
        here.parent / "HW1" / "宥宏",
        here.parent.parent / "HW1" / "宥宏",
    ]
    for d in candidates:
        if (d / "unit_SI.py").is_file():
            sys.path.insert(0, str(d))
            return
    for root in (here.parent, here.parent.parent):
        hw1 = root / "HW1"
        if hw1.is_dir():
            for found in hw1.rglob("unit_SI.py"):
                if found.name == "unit_SI.py":
                    sys.path.insert(0, str(found.parent))
                    return
    raise ImportError(
        "找不到 unit_SI.py。預期位置：腳本同層、同層「宥宏」、或專案內「HW1/宥宏」。"
    )


_ensure_vendor_path()

# 與腳本同層（MainScript）「results_test」：所有匯出之數據、圖檔
_DATA_DIR = Path(__file__).resolve().parent
DATA_RESULTS_DIR = _DATA_DIR / "results_test"
# PData／PDeltaData 文字檔：同屬 MainScript/results_test 子資料夾
PDATA_EXPORT_DIR = DATA_RESULTS_DIR / "PData"
PDELTADATA_EXPORT_DIR = DATA_RESULTS_DIR / "PDeltaData"
# fc／fy 兩兩比較圖（各 2 條 M–φ 曲線）
COMPA_DIR = _DATA_DIR / "compa"

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

# 4. 定義材料（四種柱混凝土：28／35／42／70 MPa；主程式八案：四強度 × SD420／SD550）
def func_MatDef(BCol, HCol, D_Longi, D_StirUp):
    # 材料定義參數
    fyh = 420.0 * unit_SI.MPa
    clear_cover = 40.0 * unit_SI.mm
    n_long_y = 5
    n_long_z = 5
    s_spacing = 100.0 * unit_SI.mm

    # unconfined concrete (28MPa)
    ID28MPaUC = 4001
    fpc28MPaUC = -28.0 * unit_SI.MPa
    epsc028MPaUC = -0.002
    fpcu28MPaUC = 0.0
    epsU28MPaUC = -0.004
    ops.uniaxialMaterial("Concrete01", ID28MPaUC, fpc28MPaUC, epsc028MPaUC, fpcu28MPaUC, epsU28MPaUC)

    # unconfined concrete (35MPa)
    ID35MPaUC = 4002
    fpc35MPaUC = -35.0 * unit_SI.MPa
    epsc035MPaUC = -0.002
    fpcu35MPaUC = 0.0
    epsU35MPaUC = -0.004
    ops.uniaxialMaterial("Concrete01", ID35MPaUC, fpc35MPaUC, epsc035MPaUC, fpcu35MPaUC, epsU35MPaUC)

    # unconfined concrete (42MPa)
    ID42MPaUC = 4003
    fpc42MPaUC = -42.0 * unit_SI.MPa
    epsc042MPaUC = -0.002
    fpcu42MPaUC = 0.0
    epsU42MPaUC = -0.004
    ops.uniaxialMaterial("Concrete01", ID42MPaUC, fpc42MPaUC, epsc042MPaUC, fpcu42MPaUC, epsU42MPaUC)

    # unconfined concrete (70MPa)
    ID70MPaUC = 4004
    fpc70MPaUC = -70.0 * unit_SI.MPa
    epsc070MPaUC = -0.002
    fpcu70MPaUC = 0.0
    epsU70MPaUC = -0.004
    ops.uniaxialMaterial("Concrete01", ID70MPaUC, fpc70MPaUC, epsc070MPaUC, fpcu70MPaUC, epsU70MPaUC)

    # confined concrete (28MPa)
    ID28MPaCC = 4101
    fcc_28, ecc_28, Ec_28, r_28 = func_get_mander_params(-fpc28MPaUC, fyh, BCol, HCol, clear_cover, n_long_y, n_long_z, D_Longi, D_StirUp, s_spacing)
    ops.uniaxialMaterial('Concrete01', ID28MPaCC, -fcc_28, -ecc_28, -0.2 * fcc_28, -0.02)

    # confined concrete (35MPa)
    ID35MPaCC = 4102
    fcc_35, ecc_35, Ec_35, r_35 = func_get_mander_params(-fpc35MPaUC, fyh, BCol, HCol, clear_cover, n_long_y, n_long_z, D_Longi, D_StirUp, s_spacing)
    ops.uniaxialMaterial('Concrete01', ID35MPaCC, -fcc_35, -ecc_35, -0.2 * fcc_35, -0.02)

    # confined concrete (42MPa)
    ID42MPaCC = 4103
    fcc_42, ecc_42, Ec_42, r_42 = func_get_mander_params(-fpc42MPaUC, fyh, BCol, HCol, clear_cover, n_long_y, n_long_z, D_Longi, D_StirUp, s_spacing)
    ops.uniaxialMaterial('Concrete01', ID42MPaCC, -fcc_42, -ecc_42, -0.2 * fcc_42, -0.015)

    # confined concrete (70MPa)
    ID70MPaCC = 4104
    fcc_70, ecc_70, Ec_70, r_70 = func_get_mander_params(-fpc70MPaUC, fyh, BCol, HCol, clear_cover, n_long_y, n_long_z, D_Longi, D_StirUp, s_spacing)
    ops.uniaxialMaterial('Concrete01', ID70MPaCC, -fcc_70, -ecc_70, -0.2 * fcc_70, -0.015)

    # rebar (SD420)
    IDSD420 = 4201
    fySD420 = 420 * unit_SI.MPa
    fuSD420 = fySD420 * 1.25
    EsSD420 = 200 * unit_SI.GPa
    EshSD420 = EsSD420 / 30
    bSD420 = 0.018
    R0SD420 = 0.13
    ops.uniaxialMaterial("ReinforcingSteel", IDSD420, fySD420, fuSD420, EsSD420, EshSD420, bSD420, R0SD420)

    # rebar (SD550)
    IDSD550 = 4202
    fySD550 = 550 * unit_SI.MPa
    fuSD550 = fySD550 * 1.25
    EsSD550 = 200 * unit_SI.GPa
    EshSD550 = EsSD550 / 30
    bSD550 = 0.018
    R0SD550 = 0.13
    ops.uniaxialMaterial("ReinforcingSteel", IDSD550, fySD550, fuSD550, EsSD550, EshSD550, bSD550, R0SD550)

    return fcc_28, ecc_28, fcc_35, ecc_35, fcc_42, ecc_42, fcc_70, ecc_70

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


# Foroughi & Yuksel (2020)「曲率延性」：μ_φ=φ_u/φ_y；φ_y=受拉鋼筋達 fy；
# φ_u 採 **ε_cu=0.004**（與本腳本 func_MatDef 未圍束 Concrete01 之 epsU=-0.004 一致）。
_MPHI_PHI_U_EPS = -0.004


def _mphi_cover_sample_yz(BCol, HCol, DLongi, DStirUp):
    """
    與 func_sectionDef 相同幾何，於保護層**最外緣表面**取查詢點 (y,z)（mm），供 fiber 應變讀取。
    先前若取帶中點 (zC±zE)/2，彎曲時該處 |ε| 小於外緣纖維，會高估達 εU 所需之 κ。
    """
    yE1 = -BCol / 2
    yE2 = BCol / 2
    zE1 = -HCol / 2
    zE2 = HCol / 2
    return [
        (0.0, zE2),
        (0.0, zE1),
        (yE1, zE2),
        (yE2, zE2),
        (yE1, zE1),
        (yE2, zE1),
    ]


def _fiber_strain_at_mat(ele_tag, sec_num, y, z, mat_tag):
    """以 eleResponse 讀最接近 (y,z) 且材質為 mat_tag 之纖維應變（壓為負、拉為正）。"""
    yf, zf = float(y), float(z)
    mt = int(mat_tag)
    arg_lists = (
        (ele_tag, "section", sec_num, "fiber", yf, zf, mt, "stressStrain"),
        (ele_tag, "section", sec_num, "fiber", yf, zf, mt, "strain"),
        (ele_tag, "section", "fiber", yf, zf, mt, "stressStrain"),
        (ele_tag, "section", "fiber", yf, zf, mt, "strain"),
        (ele_tag, "section", sec_num, "fiber", str(yf), str(zf), str(mt), "stressStrain"),
    )
    for args in arg_lists:
        try:
            r = ops.eleResponse(*args)
            if isinstance(r, (list, tuple)) and len(r) >= 2:
                return float(r[1])
            return float(r)
        except BaseException:
            continue
    raise RuntimeError("eleResponse fiber strain 失敗")


def _fiber_strain_at_uc(ele_tag, sec_num, y, z, id_uc):
    return _fiber_strain_at_mat(ele_tag, sec_num, y, z, id_uc)


def _min_cover_uc_strain(ele_tag, sec_num, id_uc, sample_yz):
    """多點取最小應變（最大壓縮）。任一點失敗則略過該點。"""
    strains = []
    for y, z in sample_yz:
        try:
            strains.append(_fiber_strain_at_uc(ele_tag, sec_num, y, z, id_uc))
        except BaseException:
            continue
    if not strains:
        return None
    return min(strains)


def _mphi_longitudinal_bar_yz(BCol, HCol, DLongi, ALongi, DStirUp):
    """與 func_sectionDef 相同之縱筋位置 (y,z)（mm），共 5×5 點。"""
    cover = 50 * unit_SI.mm + DStirUp + DLongi / 2
    yE1 = -BCol / 2
    spacing = (BCol - 2 * cover) / 4
    yL1 = yE1 + cover + spacing * 0
    yL2 = yE1 + cover + spacing * 1
    yL3 = yE1 + cover + spacing * 2
    yL4 = yE1 + cover + spacing * 3
    yL5 = yE1 + cover + spacing * 4
    zE1 = -HCol / 2
    zL1 = zE1 + cover + spacing * 0
    zL2 = zE1 + cover + spacing * 1
    zL3 = zE1 + cover + spacing * 2
    zL4 = zE1 + cover + spacing * 3
    zL5 = zE1 + cover + spacing * 4
    pts = []
    for y in (yL1, yL2, yL3, yL4, yL5):
        for z in (zL1, zL2, zL3, zL4, zL5):
            pts.append((y, z))
    return pts


def _max_steel_strain(ele_tag, sec_num, id_steel, bar_pts):
    """縱筋多點取 max(ε)（受拉側最先達 εy 之代表，文獻 φ_y 定義）。"""
    strains = []
    for y, z in bar_pts:
        try:
            strains.append(_fiber_strain_at_mat(ele_tag, sec_num, y, z, id_steel))
        except BaseException:
            continue
    if not strains:
        return None
    return max(strains)


def _interp_kappa_moment(k_prev, k_curr, m_prev_knm, m_curr_knm, s_prev, s_curr, s_target):
    """線性內插：應變由 s_prev→s_curr 時，跨過 s_target 之 κ(1/m) 與 M(kN·m)。"""
    denom = s_curr - s_prev
    if abs(denom) <= 1e-16:
        return float(k_curr), float(m_curr_knm)
    t = (s_target - s_prev) / denom
    kappa = k_prev + (k_curr - k_prev) * t
    moment = m_prev_knm + (m_curr_knm - m_prev_knm) * t
    return float(kappa), float(moment)


def _section_moment_bending(ele_tag=1, sec_num=1):
    """
    2D fiber section：eleResponse 回傳 [軸力, 彎矩]。
    本腳本材料以 unit_SI（mm、kN 基底）定義，斷面力與 **彎矩單位為 kN·mm**（非 N·mm）。
    """
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
    - **匯出 CSV 第一欄為 κ（1/m）**：κ(1/m) = κ(1/mm) × 1000；第二欄保留 1/mm 供核對；第三欄為彎矩 **moment_kNm（kN·m）**。
    - 回傳給繪圖之曲率欄位同為 **1/m**，彎矩欄位為 **kN·m**。
    彎矩：與 unit_SI 一致時，斷面第二分量為 **kN·mm**；**匯出與回傳之 kN·m = kN·mm ÷ 1000**（勿誤用 N·mm 之 ÷10⁶）。

    須與柱分析分開呼叫（內部會 ops.wipe），不可與 forceBeamColumn 模型並存。

    額外狀態（供延性比 μ_φ=φ_u/φ_y；擷取對齊 Foroughi & Yuksel 2020 式 (1) 說明）：
    - **φ_y、M_y**：**受拉鋼筋**纖維 **max(ε)** 首次達 **εy=fy/Es**（降伏曲率）。
    - **φ_u、M_u**：最外緣**未圍束**混凝土 **min(ε)** 首次達 **ε_cu=0.004**（模型 ε=-0.004，同 Concrete01 epsU）。
    相鄰兩步線性內插；加曲率前已達則 κ=0；未達則 None。

    回傳第三項為 dict：phi_y, My_knm, phi_u, Mu_knm（鍵名見 return）。
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

    _sec_num = 1
    _ele_tag = 1
    _cover_pts = _mphi_cover_sample_yz(BCol, HCol, DLongi, DStirUp)
    _bar_pts = _mphi_longitudinal_bar_yz(BCol, HCol, DLongi, ALongi, DStirUp)

    if IDMat_Longi == 4201:
        fy_long = 420.0 * unit_SI.MPa
    else:
        fy_long = 550.0 * unit_SI.MPa
    Es_long = 200.0 * unit_SI.GPa
    epsy = fy_long / Es_long

    kappa_yield_m = None
    moment_yield_knm = None
    kappa_spall_m = None
    moment_spall_knm = None

    def _cover_strain_min():
        return _min_cover_uc_strain(_ele_tag, _sec_num, IDUC, _cover_pts)

    def _steel_strain_max():
        return _max_steel_strain(_ele_tag, _sec_num, IDMat_Longi, _bar_pts)

    k_prev = float(ops.nodeDisp(2, 3)) * 1000.0
    m_prev_knm = float(moment_hist[0]) / 1000.0
    s_cov_prev = _cover_strain_min()
    s_steel_prev = _steel_strain_max()

    if s_steel_prev is not None and s_steel_prev >= epsy - 1e-12:
        kappa_yield_m = k_prev
        moment_yield_knm = m_prev_knm
    if s_cov_prev is not None and s_cov_prev <= _MPHI_PHI_U_EPS + 1e-12:
        kappa_spall_m = k_prev
        moment_spall_knm = m_prev_knm

    cover = 50 * unit_SI.mm + DStirUp + DLongi / 2
    d_eff = float(HCol - cover)
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
        moment_knmm = _section_moment_bending()
        moment_hist.append(moment_knmm)
        k_curr = float(ops.nodeDisp(2, 3)) * 1000.0
        m_curr_knm = float(moment_knmm) / 1000.0
        s_cov_curr = _cover_strain_min()
        s_steel_curr = _steel_strain_max()

        if kappa_yield_m is None and s_steel_prev is not None and s_steel_curr is not None:
            if s_steel_curr >= epsy - 1e-12:
                if s_steel_prev < epsy - 1e-12:
                    kappa_yield_m, moment_yield_knm = _interp_kappa_moment(
                        k_prev,
                        k_curr,
                        m_prev_knm,
                        m_curr_knm,
                        s_steel_prev,
                        s_steel_curr,
                        epsy,
                    )

        if kappa_spall_m is None and s_cov_prev is not None and s_cov_curr is not None:
            if s_cov_curr <= _MPHI_PHI_U_EPS + 1e-12:
                if s_cov_prev > _MPHI_PHI_U_EPS + 1e-12:
                    kappa_spall_m, moment_spall_knm = _interp_kappa_moment(
                        k_prev,
                        k_curr,
                        m_prev_knm,
                        m_curr_knm,
                        s_cov_prev,
                        s_cov_curr,
                        _MPHI_PHI_U_EPS,
                    )
                else:
                    kappa_spall_m = k_curr
                    moment_spall_knm = m_curr_knm

        k_prev = k_curr
        m_prev_knm = m_curr_knm
        s_cov_prev = s_cov_curr
        s_steel_prev = s_steel_curr

    data_dir = str(DATA_RESULTS_DIR / f"HW2_case{case_num}_Data")
    os.makedirs(data_dir, exist_ok=True)
    out_path = os.path.join(data_dir, f"MomentCurvature_{case_label}.csv")
    kappa_arr = np.asarray(kappa_hist, dtype=float)
    kappa_per_m = kappa_arr * 1000.0
    # 纖維斷面彎矩為 kN·mm（mm × kN 基底）；kN·m = kN·mm / 1000
    moment_knmm = np.asarray(moment_hist, dtype=float)
    moment_knm = moment_knmm / 1000.0
    header = "curvature_1_per_m,curvature_1_per_mm,moment_kNm"
    np.savetxt(
        out_path,
        np.column_stack([kappa_per_m, kappa_arr, moment_knm]),
        delimiter=",",
        header=header,
        comments="",
    )
    print(f"已匯出 moment–curvature：{out_path}")
    mphi_lim = {
        "phi_y": kappa_yield_m,
        "My_knm": moment_yield_knm,
        "phi_u": kappa_spall_m,
        "Mu_knm": moment_spall_knm,
    }
    return np.column_stack([kappa_per_m, moment_knm]), out_path, mphi_lim


def _curve_moment_knm(curve):
    """curve 第二欄為彎矩（kN·m，與 CSV 之 moment_kNm 一致）。"""
    return np.asarray(curve[:, 1], dtype=float)


def plot_moment_curvature_single(curve, title, out_png):
    """curve: (N,2) 曲率(1/m)、彎矩(kN·m)；title：圖表標題（建議含 $f_c'$、$f_y$）。"""
    k = curve[:, 0]
    m_knm = _curve_moment_knm(curve)
    x_plot = k
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    ax.plot(x_plot, m_knm, "-", linewidth=1.4, color="#1f77b4")
    ax.set_xlabel(r"Curvature $\phi$ (m$^{-1}$)")
    ax.set_ylabel("Moment (kN·m)")
    ax.set_title(title, fontsize=14)
    ax.grid(True, linestyle=":", alpha=0.7)
    ax.xaxis.set_major_locator(mticker.MaxNLocator(nbins=7, prune="both"))
    ax.xaxis.set_major_formatter(mticker.FormatStrFormatter("%.3f"))
    fig.tight_layout()
    fig.savefig(out_png, dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"已匯出 moment–curvature 圖：{out_png}")


def plot_moment_curvature_combined(series, out_png):
    """series: [(label, curve_ndarray), ...]，曲率欄位為 1/m；多案合圖。"""
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f"]
    fig, ax = plt.subplots(figsize=(11, 5.8))
    for i, (label, curve) in enumerate(series):
        if curve is None or len(curve) == 0:
            continue
        k = curve[:, 0]
        m_knm = _curve_moment_knm(curve)
        x_plot = k
        c = colors[i % len(colors)]
        ax.plot(x_plot, m_knm, "-", linewidth=1.5, color=c, label=label)
    ax.set_xlabel(r"Curvature $\phi$ (m$^{-1}$)")
    ax.set_ylabel("Moment (kN·m)")
    ax.set_title("Moment–curvature (OpenSees) all cases")
    ax.grid(True, linestyle=":", alpha=0.7)
    ax.xaxis.set_major_locator(mticker.MaxNLocator(nbins=8, prune="both"))
    ax.xaxis.set_major_formatter(mticker.FormatStrFormatter("%.3f"))
    ax.legend(loc="best", fontsize=7.5)
    fig.tight_layout()
    fig.savefig(out_png, dpi=240, bbox_inches="tight")
    plt.close(fig)
    print(f"已匯出合圖：{out_png}")


def plot_moment_curvature_compare_two(curve_a, label_a, curve_b, label_b, out_png, title_line):
    """
    同一圖內兩條 moment–curvature。
    curve_*: (N,2)，曲率 1/m、彎矩 kN·m
    title_line: 圖表主標題（已含比較類型說明）
    """
    colors = ["#1f77b4", "#ff7f0e"]
    fig, ax = plt.subplots(figsize=(7, 4.8))
    for curve, lab, c in ((curve_a, label_a, colors[0]), (curve_b, label_b, colors[1])):
        if curve is None or len(curve) == 0:
            continue
        k = curve[:, 0]
        m_knm = _curve_moment_knm(curve)
        ax.plot(k, m_knm, "-", linewidth=1.6, color=c, label=lab)
    ax.set_xlabel(r"Curvature $\phi$ (m$^{-1}$)")
    ax.set_ylabel("Moment (kN·m)")
    ax.set_title(title_line)
    ax.grid(True, linestyle=":", alpha=0.7)
    ax.xaxis.set_major_locator(mticker.MaxNLocator(nbins=8, prune="both"))
    ax.xaxis.set_major_formatter(mticker.FormatStrFormatter("%.3f"))
    ax.legend(loc="best", fontsize=9)
    fig.tight_layout()
    fig.savefig(out_png, dpi=240, bbox_inches="tight")
    plt.close(fig)
    print(f"已匯出比較圖：{out_png}")


def plot_moment_curvature_compare_three(
    curve_a,
    label_a,
    curve_b,
    label_b,
    curve_c,
    label_c,
    out_png,
    title_line,
    xlim_max=None,
):
    """
    同一圖內三條 moment–curvature（例如固定 fy、比較 f'c = 28/35/42 MPa）。
    形式對齊 plot_moment_curvature_compare_two（座標、格線、字級），三種顏色區分案別。
    xlim_max：若給定（如 0.015），橫軸限制為 [0, xlim_max]（1/m），縱軸依該曲率範圍內彎矩自動放大。
    """
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c"]
    fig, ax = plt.subplots(figsize=(7, 4.8))
    for curve, lab, c in (
        (curve_a, label_a, colors[0]),
        (curve_b, label_b, colors[1]),
        (curve_c, label_c, colors[2]),
    ):
        if curve is None or len(curve) == 0:
            continue
        k = curve[:, 0]
        m_knm = _curve_moment_knm(curve)
        ax.plot(k, m_knm, "-", linewidth=1.6, color=c, label=lab)
    ax.set_xlabel(r"Curvature $\phi$ (m$^{-1}$)")
    ax.set_ylabel("Moment (kN·m)")
    ax.set_title(title_line)
    ax.grid(True, linestyle=":", alpha=0.7)
    ax.xaxis.set_major_locator(mticker.MaxNLocator(nbins=8, prune="both"))
    ax.xaxis.set_major_formatter(mticker.FormatStrFormatter("%.3f"))
    ax.legend(loc="lower right", fontsize=9)
    if xlim_max is not None:
        ax.set_xlim(0, float(xlim_max))
        ymax = 0.0
        for curve in (curve_a, curve_b, curve_c):
            if curve is None or len(curve) == 0:
                continue
            k = np.asarray(curve[:, 0], dtype=float)
            m_knm = np.asarray(_curve_moment_knm(curve), dtype=float)
            mask = (k >= 0) & (k <= float(xlim_max))
            if np.any(mask):
                ymax = max(ymax, float(np.max(m_knm[mask])))
        if ymax > 0:
            ax.set_ylim(0, ymax * 1.05)
    fig.tight_layout()
    fig.savefig(out_png, dpi=240, bbox_inches="tight")
    plt.close(fig)
    print(f"已匯出比較圖（三曲線）：{out_png}")


def export_mphi_comparison_charts(compa_dir, curves_by_case):
    """
    curves_by_case: dict，鍵為 'Case1_28MPa_SD420' 等，值為 (N,2) M–φ 曲線（彎矩欄位 kN·m）。
    產生 6 張雙曲線比較圖至 compa_dir。
    """
    compa_dir.mkdir(parents=True, exist_ok=True)
    g = curves_by_case.get

    def pair(key_a, key_b, lab_a, lab_b, fname, title_line):
        plot_moment_curvature_compare_two(
            g(key_a),
            lab_a,
            g(key_b),
            lab_b,
            compa_dir / fname,
            title_line,
        )

    t_steel = "Moment–curvature (Different steel grade comparison)"
    t_conc = "Moment–curvature (Different concrete strength)"

    # 1. (28,420) vs (28,550)
    pair(
        "Case1_28MPa_SD420",
        "Case2_28MPa_SD550",
        r"$f_c'=28$ MPa, $f_y=420$ MPa",
        r"$f_c'=28$ MPa, $f_y=550$ MPa",
        "Mphi_fc28_fy420_vs_fy550.png",
        t_steel,
    )
    # 2. (35,420) vs (35,550)
    pair(
        "Case3_35MPa_SD420",
        "Case4_35MPa_SD550",
        r"$f_c'=35$ MPa, $f_y=420$ MPa",
        r"$f_c'=35$ MPa, $f_y=550$ MPa",
        "Mphi_fc35_fy420_vs_fy550.png",
        t_steel,
    )
    # 3. (42,420) vs (42,550)
    pair(
        "Case5_42MPa_SD420",
        "Case6_42MPa_SD550",
        r"$f_c'=42$ MPa, $f_y=420$ MPa",
        r"$f_c'=42$ MPa, $f_y=550$ MPa",
        "Mphi_fc42_fy420_vs_fy550.png",
        t_steel,
    )
    # 4. (28,420) vs (35,420)
    pair(
        "Case1_28MPa_SD420",
        "Case3_35MPa_SD420",
        r"$f_c'=28$ MPa, $f_y=420$ MPa",
        r"$f_c'=35$ MPa, $f_y=420$ MPa",
        "Mphi_fy420_fc28_vs_fc35.png",
        t_conc,
    )
    # 5. (28,420) vs (42,420)
    pair(
        "Case1_28MPa_SD420",
        "Case5_42MPa_SD420",
        r"$f_c'=28$ MPa, $f_y=420$ MPa",
        r"$f_c'=42$ MPa, $f_y=420$ MPa",
        "Mphi_fy420_fc28_vs_fc42.png",
        t_conc,
    )
    # 6. (35,420) vs (42,420)
    pair(
        "Case3_35MPa_SD420",
        "Case5_42MPa_SD420",
        r"$f_c'=35$ MPa, $f_y=420$ MPa",
        r"$f_c'=42$ MPa, $f_y=420$ MPa",
        "Mphi_fy420_fc35_vs_fc42.png",
        t_conc,
    )


def export_mphi_concrete_strength_three_way(compa_dir, curves_by_case):
    """
    固定鋼筋強度、比較三種混凝土強度（28/35/42 MPa）之 M–φ。
    - fy=420：Case 1、3、5
    - fy=550：Case 2、4、6
    標題與雙曲線混凝土比較圖相同：Moment–curvature (Different concrete strength)
    """
    compa_dir.mkdir(parents=True, exist_ok=True)
    g = curves_by_case.get
    t_conc = "Moment–curvature (Different concrete strength)"
    plot_moment_curvature_compare_three(
        g("Case1_28MPa_SD420"),
        r"$f_c'=28$ MPa, $f_y=420$ MPa",
        g("Case3_35MPa_SD420"),
        r"$f_c'=35$ MPa, $f_y=420$ MPa",
        g("Case5_42MPa_SD420"),
        r"$f_c'=42$ MPa, $f_y=420$ MPa",
        compa_dir / "Mphi_fy420_fc28_fc35_fc42.png",
        t_conc,
    )
    plot_moment_curvature_compare_three(
        g("Case2_28MPa_SD550"),
        r"$f_c'=28$ MPa, $f_y=550$ MPa",
        g("Case4_35MPa_SD550"),
        r"$f_c'=35$ MPa, $f_y=550$ MPa",
        g("Case6_42MPa_SD550"),
        r"$f_c'=42$ MPa, $f_y=550$ MPa",
        compa_dir / "Mphi_fy550_fc28_fc35_fc42.png",
        t_conc,
    )
    # 同上兩圖，橫軸 0～0.015 1/m，便於判讀初始勁度（另存 compa/detail）
    detail_dir = compa_dir / "detail"
    detail_dir.mkdir(parents=True, exist_ok=True)
    _xlim_detail = 0.015
    plot_moment_curvature_compare_three(
        g("Case1_28MPa_SD420"),
        r"$f_c'=28$ MPa, $f_y=420$ MPa",
        g("Case3_35MPa_SD420"),
        r"$f_c'=35$ MPa, $f_y=420$ MPa",
        g("Case5_42MPa_SD420"),
        r"$f_c'=42$ MPa, $f_y=420$ MPa",
        detail_dir / "Mphi_fy420_fc28_fc35_fc42.png",
        t_conc,
        xlim_max=_xlim_detail,
    )
    plot_moment_curvature_compare_three(
        g("Case2_28MPa_SD550"),
        r"$f_c'=28$ MPa, $f_y=550$ MPa",
        g("Case4_35MPa_SD550"),
        r"$f_c'=35$ MPa, $f_y=550$ MPa",
        g("Case6_42MPa_SD550"),
        r"$f_c'=42$ MPa, $f_y=550$ MPa",
        detail_dir / "Mphi_fy550_fc28_fc35_fc42.png",
        t_conc,
        xlim_max=_xlim_detail,
    )


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
    fcc_28, ecc_28, fcc_35, ecc_35, fcc_42, ecc_42, fcc_70, ecc_70 = func_MatDef(BCol, HCol, DLongi, DStirUp)
        
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
    return (
        np.array(plot_data),
        fcc_28,
        ecc_28,
        fcc_35,
        ecc_35,
        fcc_42,
        ecc_42,
        fcc_70,
        ecc_70,
    )
    



# HW2 main program --------------------------------------------------------------------------------------------------------------------------------------------#


DLongi = rebar.D25
ALongi = rebar.A25
DStirUp = rebar.D13

BCol = 600 * unit_SI.mm
HCol = 600 * unit_SI.mm
LCol = 3200 * unit_SI.mm

case1 = 1
case2 = 2
case3 = 3
case4 = 4
case5 = 5
case6 = 6
case7 = 7
case8 = 8
ID28MPaUC = 4001
ID35MPaUC = 4002
ID42MPaUC = 4003
ID70MPaUC = 4004
ID28MPaCC = 4101
ID35MPaCC = 4102
ID42MPaCC = 4103
ID70MPaCC = 4104
IDSD420 = 4201
IDSD550 = 4202
IDsec1 = 5001
IDsec2 = 5002
IDsec3 = 5003
IDsec4 = 5004
IDsec5 = 5005
IDsec6 = 5006
IDsec7 = 5007
IDsec8 = 5008

# --- 執行八種組合：28／35／42／70 MPa × SD420／SD550 ---
DATA_RESULTS_DIR.mkdir(parents=True, exist_ok=True)
print(f"輸出目錄：{DATA_RESULTS_DIR}")
print("Running Column Analyses for All Cases...")
_ret1 = func_run_column_analysis(case1, IDsec1, BCol, HCol, ID28MPaUC, ID28MPaCC, IDSD420, DLongi, ALongi, DStirUp)
PDeltaCurve1 = _ret1[0]
fcc_28c1, ecc_28c1, fcc_35c1, ecc_35c1, fcc_42c1, ecc_42c1, fcc_70c1, ecc_70c1 = _ret1[1:]

_ret2 = func_run_column_analysis(case2, IDsec2, BCol, HCol, ID28MPaUC, ID28MPaCC, IDSD550, DLongi, ALongi, DStirUp)
PDeltaCurve2 = _ret2[0]
fcc_28c2, ecc_28c2, fcc_35c2, ecc_35c2, fcc_42c2, ecc_42c2, fcc_70c2, ecc_70c2 = _ret2[1:]

_ret3 = func_run_column_analysis(case3, IDsec3, BCol, HCol, ID35MPaUC, ID35MPaCC, IDSD420, DLongi, ALongi, DStirUp)
PDeltaCurve3 = _ret3[0]
fcc_28c3, ecc_28c3, fcc_35c3, ecc_35c3, fcc_42c3, ecc_42c3, fcc_70c3, ecc_70c3 = _ret3[1:]

_ret4 = func_run_column_analysis(case4, IDsec4, BCol, HCol, ID35MPaUC, ID35MPaCC, IDSD550, DLongi, ALongi, DStirUp)
PDeltaCurve4 = _ret4[0]
fcc_28c4, ecc_28c4, fcc_35c4, ecc_35c4, fcc_42c4, ecc_42c4, fcc_70c4, ecc_70c4 = _ret4[1:]

_ret5 = func_run_column_analysis(case5, IDsec5, BCol, HCol, ID42MPaUC, ID42MPaCC, IDSD420, DLongi, ALongi, DStirUp)
PDeltaCurve5 = _ret5[0]
fcc_28c5, ecc_28c5, fcc_35c5, ecc_35c5, fcc_42c5, ecc_42c5, fcc_70c5, ecc_70c5 = _ret5[1:]

_ret6 = func_run_column_analysis(case6, IDsec6, BCol, HCol, ID42MPaUC, ID42MPaCC, IDSD550, DLongi, ALongi, DStirUp)
PDeltaCurve6 = _ret6[0]
fcc_28c6, ecc_28c6, fcc_35c6, ecc_35c6, fcc_42c6, ecc_42c6, fcc_70c6, ecc_70c6 = _ret6[1:]

_ret7 = func_run_column_analysis(case7, IDsec7, BCol, HCol, ID70MPaUC, ID70MPaCC, IDSD420, DLongi, ALongi, DStirUp)
PDeltaCurve7 = _ret7[0]
fcc_28c7, ecc_28c7, fcc_35c7, ecc_35c7, fcc_42c7, ecc_42c7, fcc_70c7, ecc_70c7 = _ret7[1:]

_ret8 = func_run_column_analysis(case8, IDsec8, BCol, HCol, ID70MPaUC, ID70MPaCC, IDSD550, DLongi, ALongi, DStirUp)
PDeltaCurve8 = _ret8[0]
fcc_28c8, ecc_28c8, fcc_35c8, ecc_35c8, fcc_42c8, ecc_42c8, fcc_70c8, ecc_70c8 = _ret8[1:]

# --- Moment–curvature（zeroLengthSection 斷面分析，概念同 MomentCurvature.py；與柱模型分開 wipe）---
# 軸力取該案圍束混凝土強度之 20%×斷面積（可依設計軸力修改）
Ag = BCol * HCol

# Moment–curvature 常軸力（與下方 func_moment_curvature_export 之 axial_load_N 相同；僅終端機顯示）
_MPHI_AXIAL_CASES = (
    (1, "28MPa_SD420", fcc_28c1),
    (2, "28MPa_SD550", fcc_28c2),
    (3, "35MPa_SD420", fcc_35c3),
    (4, "35MPa_SD550", fcc_35c4),
    (5, "42MPa_SD420", fcc_42c5),
    (6, "42MPa_SD550", fcc_42c6),
    (7, "70MPa_SD420", fcc_70c7),
    (8, "70MPa_SD550", fcc_70c8),
)
print("\n" + "-" * 80)
print("Moment–curvature 常軸力（kN；壓力為負，與 ops.load(2, axial_load_N, …) 一致）")
print("-" * 80)
for _cn, _lab, _fcc in _MPHI_AXIAL_CASES:
    _P_kN = -0.2 * abs(_fcc) * Ag
    _fcc_MPa = abs(_fcc) / unit_SI.MPa
    print(
        f"  Case {_cn} ({_lab}): P = {_P_kN:.6f} kN"
        f"  |  f'cc = {_fcc_MPa:.2f} MPa  |  |P| = {abs(_P_kN):.6f} kN"
    )
print("-" * 80 + "\n")

# 單案 M–φ 圖標題與合圖圖例：與 Case1–8（fc、fy 單位 MPa）一致
# 與 HW2/Xtract/plot_hw2_mphi_xtract.py 單線圖標題形式一致，標明 (OpenSees)
_MPHI_CHART_TITLES = (
    r"Moment–curvature (OpenSees) Case 1: $f_c'=28$ MPa, $f_y=420$ MPa",
    r"Moment–curvature (OpenSees) Case 2: $f_c'=28$ MPa, $f_y=550$ MPa",
    r"Moment–curvature (OpenSees) Case 3: $f_c'=35$ MPa, $f_y=420$ MPa",
    r"Moment–curvature (OpenSees) Case 4: $f_c'=35$ MPa, $f_y=550$ MPa",
    r"Moment–curvature (OpenSees) Case 5: $f_c'=42$ MPa, $f_y=420$ MPa",
    r"Moment–curvature (OpenSees) Case 6: $f_c'=42$ MPa, $f_y=550$ MPa",
    r"Moment–curvature (OpenSees) Case 7: $f_c'=70$ MPa, $f_y=420$ MPa",
    r"Moment–curvature (OpenSees) Case 8: $f_c'=70$ MPa, $f_y=550$ MPa",
)
_MPHI_LEGEND_LABELS = (
    r"Case 1: $f_c'=28$, $f_y=420$ MPa",
    r"Case 2: $f_c'=28$, $f_y=550$ MPa",
    r"Case 3: $f_c'=35$, $f_y=420$ MPa",
    r"Case 4: $f_c'=35$, $f_y=550$ MPa",
    r"Case 5: $f_c'=42$, $f_y=420$ MPa",
    r"Case 6: $f_c'=42$, $f_y=550$ MPa",
    r"Case 7: $f_c'=70$, $f_y=420$ MPa",
    r"Case 8: $f_c'=70$, $f_y=550$ MPa",
)
_mphi_series = []

curve1, _, mphi_lim1 = func_moment_curvature_export(
    "Case1_28MPa_SD420",
    case1,
    IDsec1,
    ID28MPaUC,
    ID28MPaCC,
    IDSD420,
    BCol,
    HCol,
    DLongi,
    ALongi,
    DStirUp,
    -0.2 * abs(fcc_28c1) * Ag,
)
plot_moment_curvature_single(
    curve1,
    _MPHI_CHART_TITLES[0],
    DATA_RESULTS_DIR / "HW2_case1_Data" / "MomentCurvature_Case1_28MPa_SD420.png",
)
_mphi_series.append((_MPHI_LEGEND_LABELS[0], curve1))

curve2, _, mphi_lim2 = func_moment_curvature_export(
    "Case2_28MPa_SD550",
    case2,
    IDsec2,
    ID28MPaUC,
    ID28MPaCC,
    IDSD550,
    BCol,
    HCol,
    DLongi,
    ALongi,
    DStirUp,
    -0.2 * abs(fcc_28c2) * Ag,
)
plot_moment_curvature_single(
    curve2,
    _MPHI_CHART_TITLES[1],
    DATA_RESULTS_DIR / "HW2_case2_Data" / "MomentCurvature_Case2_28MPa_SD550.png",
)
_mphi_series.append((_MPHI_LEGEND_LABELS[1], curve2))

curve3, _, mphi_lim3 = func_moment_curvature_export(
    "Case3_35MPa_SD420",
    case3,
    IDsec3,
    ID35MPaUC,
    ID35MPaCC,
    IDSD420,
    BCol,
    HCol,
    DLongi,
    ALongi,
    DStirUp,
    -0.2 * abs(fcc_35c3) * Ag,
)
plot_moment_curvature_single(
    curve3,
    _MPHI_CHART_TITLES[2],
    DATA_RESULTS_DIR / "HW2_case3_Data" / "MomentCurvature_Case3_35MPa_SD420.png",
)
_mphi_series.append((_MPHI_LEGEND_LABELS[2], curve3))

curve4, _, mphi_lim4 = func_moment_curvature_export(
    "Case4_35MPa_SD550",
    case4,
    IDsec4,
    ID35MPaUC,
    ID35MPaCC,
    IDSD550,
    BCol,
    HCol,
    DLongi,
    ALongi,
    DStirUp,
    -0.2 * abs(fcc_35c4) * Ag,
)
plot_moment_curvature_single(
    curve4,
    _MPHI_CHART_TITLES[3],
    DATA_RESULTS_DIR / "HW2_case4_Data" / "MomentCurvature_Case4_35MPa_SD550.png",
)
_mphi_series.append((_MPHI_LEGEND_LABELS[3], curve4))

curve5, _, mphi_lim5 = func_moment_curvature_export(
    "Case5_42MPa_SD420",
    case5,
    IDsec5,
    ID42MPaUC,
    ID42MPaCC,
    IDSD420,
    BCol,
    HCol,
    DLongi,
    ALongi,
    DStirUp,
    -0.2 * abs(fcc_42c5) * Ag,
)
plot_moment_curvature_single(
    curve5,
    _MPHI_CHART_TITLES[4],
    DATA_RESULTS_DIR / "HW2_case5_Data" / "MomentCurvature_Case5_42MPa_SD420.png",
)
_mphi_series.append((_MPHI_LEGEND_LABELS[4], curve5))

curve6, _, mphi_lim6 = func_moment_curvature_export(
    "Case6_42MPa_SD550",
    case6,
    IDsec6,
    ID42MPaUC,
    ID42MPaCC,
    IDSD550,
    BCol,
    HCol,
    DLongi,
    ALongi,
    DStirUp,
    -0.2 * abs(fcc_42c6) * Ag,
)
plot_moment_curvature_single(
    curve6,
    _MPHI_CHART_TITLES[5],
    DATA_RESULTS_DIR / "HW2_case6_Data" / "MomentCurvature_Case6_42MPa_SD550.png",
)
_mphi_series.append((_MPHI_LEGEND_LABELS[5], curve6))

curve7, _, mphi_lim7 = func_moment_curvature_export(
    "Case7_70MPa_SD420",
    case7,
    IDsec7,
    ID70MPaUC,
    ID70MPaCC,
    IDSD420,
    BCol,
    HCol,
    DLongi,
    ALongi,
    DStirUp,
    -0.2 * abs(fcc_70c7) * Ag,
)
plot_moment_curvature_single(
    curve7,
    _MPHI_CHART_TITLES[6],
    DATA_RESULTS_DIR / "HW2_case7_Data" / "MomentCurvature_Case7_70MPa_SD420.png",
)
_mphi_series.append((_MPHI_LEGEND_LABELS[6], curve7))

curve8, _, mphi_lim8 = func_moment_curvature_export(
    "Case8_70MPa_SD550",
    case8,
    IDsec8,
    ID70MPaUC,
    ID70MPaCC,
    IDSD550,
    BCol,
    HCol,
    DLongi,
    ALongi,
    DStirUp,
    -0.2 * abs(fcc_70c8) * Ag,
)
plot_moment_curvature_single(
    curve8,
    _MPHI_CHART_TITLES[7],
    DATA_RESULTS_DIR / "HW2_case8_Data" / "MomentCurvature_Case8_70MPa_SD550.png",
)
_mphi_series.append((_MPHI_LEGEND_LABELS[7], curve8))

# Moment–curvature 各 Case 極值（僅終端機；曲線欄位：φ 1/m、M kN·m）
_MPHI_EXTREMA_ROWS = (
    (1, "Case1_28MPa_SD420", curve1),
    (2, "Case2_28MPa_SD550", curve2),
    (3, "Case3_35MPa_SD420", curve3),
    (4, "Case4_35MPa_SD550", curve4),
    (5, "Case5_42MPa_SD420", curve5),
    (6, "Case6_42MPa_SD550", curve6),
    (7, "Case7_70MPa_SD420", curve7),
    (8, "Case8_70MPa_SD550", curve8),
)
print("\n" + "-" * 80)
print(
    "Moment–curvature 各 Case（OpenSees）"
    "：1. 最大曲率 φ_max 發生處之 Moment"
    "  2. 曲線上 φ_max（1/m）"
)
print("-" * 80)
for _cn, _lab, _crv in _MPHI_EXTREMA_ROWS:
    if _crv is None or len(_crv) == 0:
        print(f"  Case {_cn} ({_lab}): 無資料")
        continue
    _phi = np.asarray(_crv[:, 0], dtype=float)
    _m = np.asarray(_crv[:, 1], dtype=float)
    _ix = int(np.argmax(_phi))
    _phi_max = float(_phi[_ix])
    _m_at_phi_max = float(_m[_ix])
    print(
        f"  Case {_cn} ({_lab}):"
        f" M(φ_max處) = {_m_at_phi_max:.6f} kN·m"
        f"  |  φ_max = {_phi_max:.6f} (1/m)"
    )
print("-" * 80 + "\n")

print("\n" + "-" * 80)
print(
    "延性用極限狀態（Foroughi & Yuksel 2020 式 (1) μ_φ=φ_u/φ_y；Moment–curvature 斷面；應變跨步線性內插）"
    "\n  [1] 受拉鋼筋達 fy：1-1 斷面曲率 φ_y (1/m)，1-2 對應彎矩 M_y (kN·m)  [max(ε_steel)=εy]"
    "\n  [2] 最外緣未圍束混凝土壓應變達 εcu=0.004（模型 ε=-0.004，同 Concrete01 epsU）："
    "2-1 斷面曲率 φ_u (1/m)，2-2 對應彎矩 M_u (kN·m)  [min(ε_cover)]"
    "\n  μ = φ_u / φ_y（曲率延性比；φ_y 或 φ_u 未達則顯示「—」）"
)
print("-" * 80)
_MPHI_DUCTILITY_ROWS = (
    (1, "Case1_28MPa_SD420", mphi_lim1),
    (2, "Case2_28MPa_SD550", mphi_lim2),
    (3, "Case3_35MPa_SD420", mphi_lim3),
    (4, "Case4_35MPa_SD550", mphi_lim4),
    (5, "Case5_42MPa_SD420", mphi_lim5),
    (6, "Case6_42MPa_SD550", mphi_lim6),
    (7, "Case7_70MPa_SD420", mphi_lim7),
    (8, "Case8_70MPa_SD550", mphi_lim8),
)
for _cn, _lab, _lim in _MPHI_DUCTILITY_ROWS:
    _py = _lim["phi_y"]
    _my = _lim["My_knm"]
    _pu = _lim["phi_u"]
    _mu = _lim["Mu_knm"]
    _s1 = (
        f"1-1 φ_y = {_py:.6f} (1/m)  |  1-2 M_y = {_my:.6f} kN·m"
        if _py is not None and _my is not None
        else "1-1 φ_y = —  |  1-2 M_y = —"
    )
    _s2 = (
        f"2-1 φ_u = {_pu:.6f} (1/m)  |  2-2 M_u = {_mu:.6f} kN·m"
        if _pu is not None and _mu is not None
        else "2-1 φ_u = —  |  2-2 M_u = —"
    )
    if _py is not None and _pu is not None and abs(_py) > 1e-18:
        _mu_ratio = _pu / _py
        _s3 = f"μ = φ_u/φ_y = {_mu_ratio:.6f}"
    else:
        _s3 = "μ = φ_u/φ_y = —"
    print(f"  Case {_cn} ({_lab}):")
    print(f"    {_s1}")
    print(f"    {_s2}")
    print(f"    {_s3}")
print("-" * 80)
print(
    "附註：φ_u 採 εcu=0.004 與本腳本未圍束 Concrete01 epsU 一致。"
    "若 μ<1，表示未圍束混凝土先達 εcu，而受拉鋼筋尚未達 εy（φ_u<φ_y）。"
    "若僅軸力已使 max(ε_steel)≥εy 或 min(ε_cover)≤-0.004，則對應曲率可為 0。"
)
print()

# 合圖僅含 Case1–6（不含 70MPa 之 Case7、8）
plot_moment_curvature_combined(
    _mphi_series[:6], DATA_RESULTS_DIR / "MomentCurvature_all_cases.png"
)

_mphi_curves_for_compa = {
    "Case1_28MPa_SD420": curve1,
    "Case2_28MPa_SD550": curve2,
    "Case3_35MPa_SD420": curve3,
    "Case4_35MPa_SD550": curve4,
    "Case5_42MPa_SD420": curve5,
    "Case6_42MPa_SD550": curve6,
    "Case7_70MPa_SD420": curve7,
    "Case8_70MPa_SD550": curve8,
}
export_mphi_comparison_charts(COMPA_DIR, _mphi_curves_for_compa)
export_mphi_concrete_strength_three_way(COMPA_DIR, _mphi_curves_for_compa)

results = {
    "Case1_28MPa_SD420": (
        PDeltaCurve1,
        fcc_28c1,
        ecc_28c1,
        fcc_35c1,
        ecc_35c1,
        fcc_42c1,
        ecc_42c1,
        fcc_70c1,
        ecc_70c1,
    ),
    "Case2_28MPa_SD550": (
        PDeltaCurve2,
        fcc_28c2,
        ecc_28c2,
        fcc_35c2,
        ecc_35c2,
        fcc_42c2,
        ecc_42c2,
        fcc_70c2,
        ecc_70c2,
    ),
    "Case3_35MPa_SD420": (
        PDeltaCurve3,
        fcc_28c3,
        ecc_28c3,
        fcc_35c3,
        ecc_35c3,
        fcc_42c3,
        ecc_42c3,
        fcc_70c3,
        ecc_70c3,
    ),
    "Case4_35MPa_SD550": (
        PDeltaCurve4,
        fcc_28c4,
        ecc_28c4,
        fcc_35c4,
        ecc_35c4,
        fcc_42c4,
        ecc_42c4,
        fcc_70c4,
        ecc_70c4,
    ),
    "Case5_42MPa_SD420": (
        PDeltaCurve5,
        fcc_28c5,
        ecc_28c5,
        fcc_35c5,
        ecc_35c5,
        fcc_42c5,
        ecc_42c5,
        fcc_70c5,
        ecc_70c5,
    ),
    "Case6_42MPa_SD550": (
        PDeltaCurve6,
        fcc_28c6,
        ecc_28c6,
        fcc_35c6,
        ecc_35c6,
        fcc_42c6,
        ecc_42c6,
        fcc_70c6,
        ecc_70c6,
    ),
    "Case7_70MPa_SD420": (
        PDeltaCurve7,
        fcc_28c7,
        ecc_28c7,
        fcc_35c7,
        ecc_35c7,
        fcc_42c7,
        ecc_42c7,
        fcc_70c7,
        ecc_70c7,
    ),
    "Case8_70MPa_SD550": (
        PDeltaCurve8,
        fcc_28c8,
        ecc_28c8,
        fcc_35c8,
        ecc_35c8,
        fcc_42c8,
        ecc_42c8,
        fcc_70c8,
        ecc_70c8,
    ),
}

# --- 輸出各 Case 曲線數據到 TXT 檔（PData／PDeltaData → MainScript/results_test/...）---
PDATA_EXPORT_DIR.mkdir(parents=True, exist_ok=True)
PDELTADATA_EXPORT_DIR.mkdir(parents=True, exist_ok=True)
print("正在輸出數據檔案...")
print(f"PData 目錄：{PDATA_EXPORT_DIR}")
print(f"PDeltaData 目錄：{PDELTADATA_EXPORT_DIR}")
for label, (
    data,
    fcc28,
    ecc28,
    fcc35,
    ecc35,
    fcc42,
    ecc42,
    fcc70,
    ecc70,
) in results.items():
    if len(data) > 0:

        fcc_28MPa = fcc28 / unit_SI.MPa
        fcc_35MPa = fcc35 / unit_SI.MPa
        fcc_42MPa = fcc42 / unit_SI.MPa
        fcc_70MPa = fcc70 / unit_SI.MPa

        hdr = (
            f"f'cc28MPa = {fcc_28MPa:.2f} MPa, ecc28MPa = {ecc28:.4f}\n"
            f"f'cc35MPa = {fcc_35MPa:.2f} MPa, ecc35MPa = {ecc35:.4f}\n"
            f"f'cc42MPa = {fcc_42MPa:.2f} MPa, ecc42MPa = {ecc42:.4f}\n"
            f"f'cc70MPa = {fcc_70MPa:.2f} MPa, ecc70MPa = {ecc70:.4f}\n"
            "Axial_Strain Axial_Load(kN)"
        )

        # 儲存為 txt 檔，第一欄是 Strain，第二欄是 Axial Load (kN)
        filename = str(PDELTADATA_EXPORT_DIR / f"PDeltaData_{label}.txt")
        header = hdr
        np.savetxt(filename, data, fmt='%12.6e', header=header, comments='')
        print("************************************")
        print(f"已儲存: {filename}")
        print(f"f'cc = {fcc_28MPa:.2f} MPa, ecc = {ecc28:.4f}")
        print(f"f'cc = {fcc_35MPa:.2f} MPa, ecc = {ecc35:.4f}")
        print(f"f'cc = {fcc_42MPa:.2f} MPa, ecc = {ecc42:.4f}")
        print(f"f'cc = {fcc_70MPa:.2f} MPa, ecc = {ecc70:.4f}")

        # 儲存為 txt 檔，第一欄是 Strain，第二欄是 Axial Load (kN)
        filename = str(PDATA_EXPORT_DIR / f"PData_{label}.txt")
        header = hdr
        np.savetxt(filename, data, fmt='%12.6e', header=header, comments='')
        print("************************************")
        print(f"已儲存: {filename}")
        print(f"f'cc = {fcc_28MPa:.2f} MPa, ecc = {ecc28:.4f}")
        print(f"f'cc = {fcc_35MPa:.2f} MPa, ecc = {ecc35:.4f}")
        print(f"f'cc = {fcc_42MPa:.2f} MPa, ecc = {ecc42:.4f}")
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



